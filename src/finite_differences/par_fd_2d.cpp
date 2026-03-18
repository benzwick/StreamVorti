/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors (alphabetically):
 *      George C. BOURANTAS
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

#include "StreamVorti/finite_differences/par_fd_2d.hpp"
#include "StreamVorti/finite_differences/fd.hpp"

#ifdef MFEM_USE_MPI

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace StreamVorti {

ParFiniteDiff2d::ParFiniteDiff2d(mfem::ParGridFunction &gf, int stencil_order)
    : ParFiniteDiff(gf, stencil_order),
      sh_func_dx_(nullptr),
      sh_func_dy_(nullptr),
      sh_func_dxx_(nullptr),
      sh_func_dyy_(nullptr),
      sh_func_dxy_(nullptr)
{}

ParFiniteDiff2d::~ParFiniteDiff2d()
{
    delete sh_func_dx_;
    delete sh_func_dy_;
    delete sh_func_dxx_;
    delete sh_func_dyy_;
    delete sh_func_dxy_;
}


void ParFiniteDiff2d::Update()
{
    mfem::StopWatch timer;
    timer.Start();
    if (rank_ == 0)
        std::cout << "ParFD: update 2D derivative matrices (order "
                  << stencil_order_ << ")" << std::endl;

    // -----------------------------------------------------------------
    // 1. Gather ALL node coordinates across all ranks to detect the
    //    global structured grid.  FD on structured grids is inherently
    //    a global operation (unlike DCPSE which is local).
    //    For very large grids this could be replaced by a distributed
    //    algorithm, but for the grid sizes used in practice (<10M nodes)
    //    the serial detection is adequate.
    // -----------------------------------------------------------------

    mfem::ParMesh *pmesh = pfes_->GetParMesh();
    int local_ndofs = pfes_->GetNDofs();

    // Gather local vertex coordinates
    std::vector<double> local_x(local_ndofs), local_y(local_ndofs);
    for (int i = 0; i < local_ndofs; ++i)
    {
        // Use mesh vertices for order-1 H1 spaces
        if (i < pmesh->GetNV())
        {
            const double *vtx = pmesh->GetVertex(i);
            local_x[i] = vtx[0];
            local_y[i] = vtx[1];
        }
    }

    // Gather global DOF offsets
    HYPRE_BigInt *tdof_offsets = pfes_->GetTrueDofOffsets();
    int local_true_dofs = pfes_->GetTrueVSize();
    HYPRE_BigInt global_ndofs = pfes_->GlobalTrueVSize();

    // -----------------------------------------------------------------
    // 2. Detect global grid structure.
    //    Each rank finds its local unique coordinates, then we do an
    //    Allreduce to find the global unique set.
    // -----------------------------------------------------------------

    // Find local unique x and y
    std::vector<double> local_ux, local_uy;
    FiniteDiff::FindUniqueCoordinates(local_x, local_ux);
    FiniteDiff::FindUniqueCoordinates(local_y, local_uy);

    // Gather all unique coordinates to all ranks
    int local_nux = static_cast<int>(local_ux.size());
    int local_nuy = static_cast<int>(local_uy.size());

    // Gather sizes
    std::vector<int> all_nux(nranks_), all_nuy(nranks_);
    MPI_Allgather(&local_nux, 1, MPI_INT, all_nux.data(), 1, MPI_INT, comm_);
    MPI_Allgather(&local_nuy, 1, MPI_INT, all_nuy.data(), 1, MPI_INT, comm_);

    // Gather coordinates
    std::vector<int> displs_x(nranks_), displs_y(nranks_);
    int total_x = 0, total_y = 0;
    for (int r = 0; r < nranks_; ++r)
    {
        displs_x[r] = total_x;
        displs_y[r] = total_y;
        total_x += all_nux[r];
        total_y += all_nuy[r];
    }

    std::vector<double> all_x(total_x), all_y(total_y);
    MPI_Allgatherv(local_ux.data(), local_nux, MPI_DOUBLE,
                   all_x.data(), all_nux.data(), displs_x.data(),
                   MPI_DOUBLE, comm_);
    MPI_Allgatherv(local_uy.data(), local_nuy, MPI_DOUBLE,
                   all_y.data(), all_nuy.data(), displs_y.data(),
                   MPI_DOUBLE, comm_);

    // Find global unique coordinates
    std::vector<double> ux, uy;
    FiniteDiff::FindUniqueCoordinates(all_x, ux);
    FiniteDiff::FindUniqueCoordinates(all_y, uy);

    int nx = static_cast<int>(ux.size());
    int ny = static_cast<int>(uy.size());

    if (static_cast<HYPRE_BigInt>(nx) * ny != global_ndofs)
    {
        MFEM_ABORT("ParFD: Grid is not structured rectangular. "
                   "Detected " << nx << " x " << ny << " = "
                   << static_cast<HYPRE_BigInt>(nx) * ny
                   << " grid points, but global DOF count is "
                   << global_ndofs << ".");
    }

    int min_nodes = (stencil_order_ == 4) ? 6 : 4;
    if (nx < min_nodes || ny < min_nodes)
    {
        MFEM_ABORT("ParFD: Grid too small for order-" << stencil_order_
                   << " stencils. Need at least " << min_nodes
                   << " nodes per direction, got " << nx << " x " << ny);
    }

    double hx = FiniteDiff::VerifyUniformSpacing(ux, "x");
    double hy = FiniteDiff::VerifyUniformSpacing(uy, "y");

    if (rank_ == 0)
    {
        std::cout << "ParFD: Global grid: " << nx << " x " << ny
                  << ", hx=" << hx << ", hy=" << hy << std::endl;
    }

    // -----------------------------------------------------------------
    // 3. For each local true DOF, compute its global grid indices (gix, giy)
    //    and assemble the FD stencil into a local SparseMatrix (size:
    //    local_true_dofs x global_ndofs, later split into diag/offd).
    // -----------------------------------------------------------------

    double xtol = 0.5 * hx;
    double ytol = 0.5 * hy;

    // Build global grid-index to global-DOF mapping.
    // For a structured mesh with lexicographic DOF ordering, global DOF
    // = gix + giy * nx.  However, we must verify this matches MFEM's
    // actual global ordering.  For parallel H1 order-1 on Cartesian
    // meshes, MFEM uses this ordering.
    // We'll use MFEM's own global DOF numbers via GetGlobalTDofNumber().

    // Assemble into full-width SparseMatrix, then let MFEM's
    // HypreParMatrix constructor handle the diag/offd split.
    // We use MFEM's RAP or direct assembly approach.

    // Actually, the cleanest approach: assemble a serial SparseMatrix on
    // each rank for its own rows, then construct HypreParMatrix.

    mfem::SparseMatrix local_dx(local_true_dofs, static_cast<int>(global_ndofs));
    mfem::SparseMatrix local_dy(local_true_dofs, static_cast<int>(global_ndofs));
    mfem::SparseMatrix local_dxx(local_true_dofs, static_cast<int>(global_ndofs));
    mfem::SparseMatrix local_dyy(local_true_dofs, static_cast<int>(global_ndofs));
    mfem::SparseMatrix local_dxy(local_true_dofs, static_cast<int>(global_ndofs));

    // For each local true DOF, find its global grid position
    for (int ldof = 0; ldof < pmesh->GetNV(); ++ldof)
    {
        int tdof = pfes_->GetLocalTDofNumber(ldof);
        if (tdof < 0) continue; // Not owned by this rank

        const double *vtx = pmesh->GetVertex(ldof);
        int gix = FiniteDiff::CoordinateToIndex(vtx[0], ux, xtol);
        int giy = FiniteDiff::CoordinateToIndex(vtx[1], uy, ytol);

        if (gix < 0 || giy < 0) continue;

        // Global DOF index for structured grid: gix + giy * nx
        // (lexicographic ordering)
        auto global_dof = [&](int ix, int iy) -> int {
            return ix + iy * nx;
        };

        int row = tdof;

        // --- dx: 2nd-order central/one-sided ---
        if (gix > 0 && gix < nx - 1)
        {
            local_dx.Add(row, global_dof(gix+1, giy),  1.0 / (2.0 * hx));
            local_dx.Add(row, global_dof(gix-1, giy), -1.0 / (2.0 * hx));
        }
        else if (gix == 0)
        {
            local_dx.Add(row, global_dof(0, giy), -3.0 / (2.0 * hx));
            local_dx.Add(row, global_dof(1, giy),  4.0 / (2.0 * hx));
            local_dx.Add(row, global_dof(2, giy), -1.0 / (2.0 * hx));
        }
        else
        {
            local_dx.Add(row, global_dof(nx-1, giy),  3.0 / (2.0 * hx));
            local_dx.Add(row, global_dof(nx-2, giy), -4.0 / (2.0 * hx));
            local_dx.Add(row, global_dof(nx-3, giy),  1.0 / (2.0 * hx));
        }

        // --- dy ---
        if (giy > 0 && giy < ny - 1)
        {
            local_dy.Add(row, global_dof(gix, giy+1),  1.0 / (2.0 * hy));
            local_dy.Add(row, global_dof(gix, giy-1), -1.0 / (2.0 * hy));
        }
        else if (giy == 0)
        {
            local_dy.Add(row, global_dof(gix, 0), -3.0 / (2.0 * hy));
            local_dy.Add(row, global_dof(gix, 1),  4.0 / (2.0 * hy));
            local_dy.Add(row, global_dof(gix, 2), -1.0 / (2.0 * hy));
        }
        else
        {
            local_dy.Add(row, global_dof(gix, ny-1),  3.0 / (2.0 * hy));
            local_dy.Add(row, global_dof(gix, ny-2), -4.0 / (2.0 * hy));
            local_dy.Add(row, global_dof(gix, ny-3),  1.0 / (2.0 * hy));
        }

        // --- dxx ---
        double h2x = hx * hx;
        if (gix > 0 && gix < nx - 1)
        {
            local_dxx.Add(row, global_dof(gix+1, giy),  1.0 / h2x);
            local_dxx.Add(row, global_dof(gix,   giy), -2.0 / h2x);
            local_dxx.Add(row, global_dof(gix-1, giy),  1.0 / h2x);
        }
        else if (gix == 0)
        {
            local_dxx.Add(row, global_dof(0, giy),  2.0 / h2x);
            local_dxx.Add(row, global_dof(1, giy), -5.0 / h2x);
            local_dxx.Add(row, global_dof(2, giy),  4.0 / h2x);
            local_dxx.Add(row, global_dof(3, giy), -1.0 / h2x);
        }
        else
        {
            local_dxx.Add(row, global_dof(nx-1, giy),  2.0 / h2x);
            local_dxx.Add(row, global_dof(nx-2, giy), -5.0 / h2x);
            local_dxx.Add(row, global_dof(nx-3, giy),  4.0 / h2x);
            local_dxx.Add(row, global_dof(nx-4, giy), -1.0 / h2x);
        }

        // --- dyy ---
        double h2y = hy * hy;
        if (giy > 0 && giy < ny - 1)
        {
            local_dyy.Add(row, global_dof(gix, giy+1),  1.0 / h2y);
            local_dyy.Add(row, global_dof(gix, giy  ), -2.0 / h2y);
            local_dyy.Add(row, global_dof(gix, giy-1),  1.0 / h2y);
        }
        else if (giy == 0)
        {
            local_dyy.Add(row, global_dof(gix, 0),  2.0 / h2y);
            local_dyy.Add(row, global_dof(gix, 1), -5.0 / h2y);
            local_dyy.Add(row, global_dof(gix, 2),  4.0 / h2y);
            local_dyy.Add(row, global_dof(gix, 3), -1.0 / h2y);
        }
        else
        {
            local_dyy.Add(row, global_dof(gix, ny-1),  2.0 / h2y);
            local_dyy.Add(row, global_dof(gix, ny-2), -5.0 / h2y);
            local_dyy.Add(row, global_dof(gix, ny-3),  4.0 / h2y);
            local_dyy.Add(row, global_dof(gix, ny-4), -1.0 / h2y);
        }

        // --- dxy ---
        {
            int ixm = (gix > 0) ? gix - 1 : gix;
            int ixp = (gix < nx - 1) ? gix + 1 : gix;
            int iym = (giy > 0) ? giy - 1 : giy;
            int iyp = (giy < ny - 1) ? giy + 1 : giy;
            double ds_x = (ixp - ixm) * hx;
            double ds_y = (iyp - iym) * hy;
            if (ds_x > 0 && ds_y > 0)
            {
                double s = 1.0 / (ds_x * ds_y);
                local_dxy.Add(row, global_dof(ixp, iyp),  s);
                local_dxy.Add(row, global_dof(ixp, iym), -s);
                local_dxy.Add(row, global_dof(ixm, iyp), -s);
                local_dxy.Add(row, global_dof(ixm, iym),  s);
            }
        }
    }

    local_dx.Finalize();
    local_dy.Finalize();
    local_dxx.Finalize();
    local_dyy.Finalize();
    local_dxy.Finalize();

    // -----------------------------------------------------------------
    // 4. Construct HypreParMatrix from the local rows.
    //    We use the constructor that takes row/col partitioning and a
    //    local SparseMatrix with global column indices.
    // -----------------------------------------------------------------

    auto make_par = [&](mfem::SparseMatrix &local_mat) -> mfem::HypreParMatrix*
    {
        return new mfem::HypreParMatrix(
            comm_, local_true_dofs, global_ndofs, global_ndofs,
            local_mat.GetI(), local_mat.GetJ(), local_mat.GetData(),
            tdof_offsets, tdof_offsets);
    };

    sh_func_dx_  = make_par(local_dx);
    sh_func_dy_  = make_par(local_dy);
    sh_func_dxx_ = make_par(local_dxx);
    sh_func_dyy_ = make_par(local_dyy);
    sh_func_dxy_ = make_par(local_dxy);

    if (rank_ == 0)
    {
        std::cout << "ParFD: Execution time for parallel 2D FD matrices: "
                  << timer.RealTime() << " s" << std::endl;
    }
}


void ParFiniteDiff2d::SaveDerivToFile(const std::string &deriv,
                                      const std::string &filename) const
{
    // Parallel matrix save not yet implemented
    if (rank_ == 0)
    {
        std::cout << "ParFiniteDiff2d::SaveDerivToFile: "
                  << "Parallel matrix export not yet implemented for '"
                  << deriv << "'." << std::endl;
    }
}

} // namespace StreamVorti

#endif // MFEM_USE_MPI

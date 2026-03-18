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

#include "StreamVorti/finite_differences/fd_3d.hpp"

#include <filesystem>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace StreamVorti {

void FiniteDiff3d::Update()
{
    mfem::StopWatch timer;
    timer.Start();
    std::cout << "FD: update 3D derivative matrices" << std::endl;

    const mfem::FiniteElementSpace *fes = gf_->FESpace();
    const mfem::Mesh *mesh = fes->GetMesh();
    int nnodes = fes->GetNDofs();

    // Extract node coordinates
    std::vector<double> xcoords(nnodes), ycoords(nnodes), zcoords(nnodes);
    for (int i = 0; i < nnodes; ++i)
    {
        const double *vtx = mesh->GetVertex(i);
        xcoords[i] = vtx[0];
        ycoords[i] = vtx[1];
        zcoords[i] = vtx[2];
    }

    // Find unique coordinates
    std::vector<double> ux(xcoords), uy(ycoords), uz(zcoords);
    std::sort(ux.begin(), ux.end());
    std::sort(uy.begin(), uy.end());
    std::sort(uz.begin(), uz.end());
    auto near = [](double a, double b){ return std::abs(a - b) < 1e-12; };
    ux.erase(std::unique(ux.begin(), ux.end(), near), ux.end());
    uy.erase(std::unique(uy.begin(), uy.end(), near), uy.end());
    uz.erase(std::unique(uz.begin(), uz.end(), near), uz.end());

    int nxn = static_cast<int>(ux.size());
    int nyn = static_cast<int>(uy.size());
    int nzn = static_cast<int>(uz.size());

    if (nxn * nyn * nzn != nnodes)
    {
        MFEM_ABORT("FD: Grid is not structured. "
                   "Expected " << nxn << " x " << nyn << " x " << nzn << " = "
                   << nxn * nyn * nzn << " nodes, but got " << nnodes);
    }

    double hx = (ux.back() - ux.front()) / (nxn - 1);
    double hy = (uy.back() - uy.front()) / (nyn - 1);
    double hz = (uz.back() - uz.front()) / (nzn - 1);

    std::cout << "FD: Grid dimensions: " << nxn << " x " << nyn << " x " << nzn << std::endl;
    std::cout << "FD: Grid spacing: hx = " << hx << ", hy = " << hy << ", hz = " << hz << std::endl;

    // Build 3D index map
    std::vector<int> dof_ix(nnodes), dof_iy(nnodes), dof_iz(nnodes);
    // Flatten 3D index: grid_to_dof[ix * nyn * nzn + iy * nzn + iz]
    std::vector<int> grid_to_dof(nxn * nyn * nzn, -1);

    auto grid_idx = [&](int ix, int iy, int iz) -> int {
        return ix * nyn * nzn + iy * nzn + iz;
    };

    for (int i = 0; i < nnodes; ++i)
    {
        auto itx = std::lower_bound(ux.begin(), ux.end(), xcoords[i] - 1e-12);
        auto ity = std::lower_bound(uy.begin(), uy.end(), ycoords[i] - 1e-12);
        auto itz = std::lower_bound(uz.begin(), uz.end(), zcoords[i] - 1e-12);
        int ix = std::min(static_cast<int>(itx - ux.begin()), nxn - 1);
        int iy = std::min(static_cast<int>(ity - uy.begin()), nyn - 1);
        int iz = std::min(static_cast<int>(itz - uz.begin()), nzn - 1);

        dof_ix[i] = ix;
        dof_iy[i] = iy;
        dof_iz[i] = iz;
        grid_to_dof[grid_idx(ix, iy, iz)] = i;
    }

    // Verify completeness
    for (int idx = 0; idx < nxn * nyn * nzn; ++idx)
        if (grid_to_dof[idx] < 0)
            MFEM_ABORT("FD: Missing node in 3D grid");

    // Helper to get DOF from grid indices
    auto dof = [&](int ix, int iy, int iz) -> int {
        return grid_to_dof[grid_idx(ix, iy, iz)];
    };

    // Initialize matrices
    this->sh_func_dx_  = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dy_  = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dz_  = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxx_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dzz_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxz_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyz_ = mfem::SparseMatrix(nnodes, nnodes);

    for (int node_id = 0; node_id < nnodes; ++node_id)
    {
        int ix = dof_ix[node_id];
        int iy = dof_iy[node_id];
        int iz = dof_iz[node_id];

        // --- dx ---
        if (ix > 0 && ix < nxn - 1)
        {
            sh_func_dx_.Add(node_id, dof(ix+1, iy, iz),  1.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, dof(ix-1, iy, iz), -1.0 / (2.0 * hx));
        }
        else if (ix == 0)
        {
            sh_func_dx_.Add(node_id, dof(ix,   iy, iz), -3.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, dof(ix+1, iy, iz),  4.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, dof(ix+2, iy, iz), -1.0 / (2.0 * hx));
        }
        else
        {
            sh_func_dx_.Add(node_id, dof(ix,   iy, iz),  3.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, dof(ix-1, iy, iz), -4.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, dof(ix-2, iy, iz),  1.0 / (2.0 * hx));
        }

        // --- dy ---
        if (iy > 0 && iy < nyn - 1)
        {
            sh_func_dy_.Add(node_id, dof(ix, iy+1, iz),  1.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, dof(ix, iy-1, iz), -1.0 / (2.0 * hy));
        }
        else if (iy == 0)
        {
            sh_func_dy_.Add(node_id, dof(ix, iy,   iz), -3.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, dof(ix, iy+1, iz),  4.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, dof(ix, iy+2, iz), -1.0 / (2.0 * hy));
        }
        else
        {
            sh_func_dy_.Add(node_id, dof(ix, iy,   iz),  3.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, dof(ix, iy-1, iz), -4.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, dof(ix, iy-2, iz),  1.0 / (2.0 * hy));
        }

        // --- dz ---
        if (iz > 0 && iz < nzn - 1)
        {
            sh_func_dz_.Add(node_id, dof(ix, iy, iz+1),  1.0 / (2.0 * hz));
            sh_func_dz_.Add(node_id, dof(ix, iy, iz-1), -1.0 / (2.0 * hz));
        }
        else if (iz == 0)
        {
            sh_func_dz_.Add(node_id, dof(ix, iy, iz  ), -3.0 / (2.0 * hz));
            sh_func_dz_.Add(node_id, dof(ix, iy, iz+1),  4.0 / (2.0 * hz));
            sh_func_dz_.Add(node_id, dof(ix, iy, iz+2), -1.0 / (2.0 * hz));
        }
        else
        {
            sh_func_dz_.Add(node_id, dof(ix, iy, iz  ),  3.0 / (2.0 * hz));
            sh_func_dz_.Add(node_id, dof(ix, iy, iz-1), -4.0 / (2.0 * hz));
            sh_func_dz_.Add(node_id, dof(ix, iy, iz-2),  1.0 / (2.0 * hz));
        }

        // --- dxx ---
        if (ix > 0 && ix < nxn - 1)
        {
            sh_func_dxx_.Add(node_id, dof(ix+1, iy, iz),  1.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix,   iy, iz), -2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix-1, iy, iz),  1.0 / (hx * hx));
        }
        else if (ix == 0)
        {
            sh_func_dxx_.Add(node_id, dof(ix,   iy, iz),  2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix+1, iy, iz), -5.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix+2, iy, iz),  4.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix+3, iy, iz), -1.0 / (hx * hx));
        }
        else
        {
            sh_func_dxx_.Add(node_id, dof(ix,   iy, iz),  2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix-1, iy, iz), -5.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix-2, iy, iz),  4.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, dof(ix-3, iy, iz), -1.0 / (hx * hx));
        }

        // --- dyy ---
        if (iy > 0 && iy < nyn - 1)
        {
            sh_func_dyy_.Add(node_id, dof(ix, iy+1, iz),  1.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy,   iz), -2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy-1, iz),  1.0 / (hy * hy));
        }
        else if (iy == 0)
        {
            sh_func_dyy_.Add(node_id, dof(ix, iy,   iz),  2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy+1, iz), -5.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy+2, iz),  4.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy+3, iz), -1.0 / (hy * hy));
        }
        else
        {
            sh_func_dyy_.Add(node_id, dof(ix, iy,   iz),  2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy-1, iz), -5.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy-2, iz),  4.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, dof(ix, iy-3, iz), -1.0 / (hy * hy));
        }

        // --- dzz ---
        if (iz > 0 && iz < nzn - 1)
        {
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz+1),  1.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz  ), -2.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz-1),  1.0 / (hz * hz));
        }
        else if (iz == 0)
        {
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz  ),  2.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz+1), -5.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz+2),  4.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz+3), -1.0 / (hz * hz));
        }
        else
        {
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz  ),  2.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz-1), -5.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz-2),  4.0 / (hz * hz));
            sh_func_dzz_.Add(node_id, dof(ix, iy, iz-3), -1.0 / (hz * hz));
        }

        // --- dxy ---
        {
            int ixm = (ix > 0) ? ix - 1 : ix;
            int ixp = (ix < nxn - 1) ? ix + 1 : ix;
            int iym = (iy > 0) ? iy - 1 : iy;
            int iyp = (iy < nyn - 1) ? iy + 1 : iy;
            double dx_step = (ixp - ixm) * hx;
            double dy_step = (iyp - iym) * hy;
            if (dx_step > 0 && dy_step > 0)
            {
                double s = 1.0 / (dx_step * dy_step);
                sh_func_dxy_.Add(node_id, dof(ixp, iyp, iz),  s);
                sh_func_dxy_.Add(node_id, dof(ixp, iym, iz), -s);
                sh_func_dxy_.Add(node_id, dof(ixm, iyp, iz), -s);
                sh_func_dxy_.Add(node_id, dof(ixm, iym, iz),  s);
            }
        }

        // --- dxz ---
        {
            int ixm = (ix > 0) ? ix - 1 : ix;
            int ixp = (ix < nxn - 1) ? ix + 1 : ix;
            int izm = (iz > 0) ? iz - 1 : iz;
            int izp = (iz < nzn - 1) ? iz + 1 : iz;
            double dx_step = (ixp - ixm) * hx;
            double dz_step = (izp - izm) * hz;
            if (dx_step > 0 && dz_step > 0)
            {
                double s = 1.0 / (dx_step * dz_step);
                sh_func_dxz_.Add(node_id, dof(ixp, iy, izp),  s);
                sh_func_dxz_.Add(node_id, dof(ixp, iy, izm), -s);
                sh_func_dxz_.Add(node_id, dof(ixm, iy, izp), -s);
                sh_func_dxz_.Add(node_id, dof(ixm, iy, izm),  s);
            }
        }

        // --- dyz ---
        {
            int iym = (iy > 0) ? iy - 1 : iy;
            int iyp = (iy < nyn - 1) ? iy + 1 : iy;
            int izm = (iz > 0) ? iz - 1 : iz;
            int izp = (iz < nzn - 1) ? iz + 1 : iz;
            double dy_step = (iyp - iym) * hy;
            double dz_step = (izp - izm) * hz;
            if (dy_step > 0 && dz_step > 0)
            {
                double s = 1.0 / (dy_step * dz_step);
                sh_func_dyz_.Add(node_id, dof(ix, iyp, izp),  s);
                sh_func_dyz_.Add(node_id, dof(ix, iyp, izm), -s);
                sh_func_dyz_.Add(node_id, dof(ix, iym, izp), -s);
                sh_func_dyz_.Add(node_id, dof(ix, iym, izm),  s);
            }
        }
    }

    // Finalize
    this->sh_func_dx_.Finalize();
    this->sh_func_dy_.Finalize();
    this->sh_func_dz_.Finalize();
    this->sh_func_dxx_.Finalize();
    this->sh_func_dyy_.Finalize();
    this->sh_func_dzz_.Finalize();
    this->sh_func_dxy_.Finalize();
    this->sh_func_dxz_.Finalize();
    this->sh_func_dyz_.Finalize();

    std::cout << "FD: Execution time for 3D finite difference matrices: "
              << timer.RealTime() << " s" << std::endl;
}


void FiniteDiff3d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    mfem::SparseMatrix derivative;

    if      (deriv == "dx")  { derivative = this->sh_func_dx_; }
    else if (deriv == "dy")  { derivative = this->sh_func_dy_; }
    else if (deriv == "dz")  { derivative = this->sh_func_dz_; }
    else if (deriv == "dxx") { derivative = this->sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = this->sh_func_dyy_; }
    else if (deriv == "dzz") { derivative = this->sh_func_dzz_; }
    else if (deriv == "dxy") { derivative = this->sh_func_dxy_; }
    else if (deriv == "dxz") { derivative = this->sh_func_dxz_; }
    else if (deriv == "dyz") { derivative = this->sh_func_dyz_; }

    if (derivative.Height() == 0)
    {
        MFEM_ABORT( "Could not save FD derivative. "
                    "Derivatives have not been computed." );
    }

    std::string path = "";
    std::size_t last_slash = filename.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }

    std::filesystem::path dir(path);
    if (!path.empty() && !std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }

    std::string ext = "";
    if (filename.find_last_of(".") != std::string::npos) {
        ext = filename.substr(filename.find_last_of("."));
    }

    std::string out_filename;
    if (ext != ".txt") { out_filename = filename + ".txt"; }
    else { out_filename = filename; }

    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    derivative.PrintMatlab(out);
    out.close();
}

} // namespace StreamVorti

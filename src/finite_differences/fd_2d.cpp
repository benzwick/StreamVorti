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

#include "StreamVorti/finite_differences/fd_2d.hpp"

#include <filesystem>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace StreamVorti {

void FiniteDiff2d::Update()
{
    mfem::StopWatch timer;
    timer.Start();
    std::cout << "FD: update 2D derivative matrices" << std::endl;

    const mfem::FiniteElementSpace *fes = gf_->FESpace();
    const mfem::Mesh *mesh = fes->GetMesh();
    int nnodes = fes->GetNDofs();

    // Determine grid dimensions from mesh.
    // For a Cartesian grid with nx*ny elements, there are (nx+1)*(ny+1) nodes.
    // We detect nx by finding the number of nodes along the x-axis (first row).
    //
    // Strategy: sort nodes by coordinates to build a lexicographic index map.
    // The DOF ordering in MFEM's H1 space on Cartesian meshes is already
    // lexicographic for order 1, but we don't rely on that assumption.

    // Extract node coordinates
    mfem::GridFunction nodes(const_cast<mfem::FiniteElementSpace*>(fes));
    mesh->GetNodes(nodes);

    // If GetNodes returned empty (order-1 mesh without explicit nodes),
    // build coordinates from the mesh vertices.
    std::vector<double> xcoords(nnodes), ycoords(nnodes);
    if (nodes.Size() == 0)
    {
        // For order-1 H1, DOFs correspond to mesh vertices
        for (int i = 0; i < nnodes; ++i)
        {
            const double *vtx = mesh->GetVertex(i);
            xcoords[i] = vtx[0];
            ycoords[i] = vtx[1];
        }
    }
    else
    {
        for (int i = 0; i < nnodes; ++i)
        {
            xcoords[i] = nodes(fes->DofToVDof(i, 0));
            ycoords[i] = nodes(fes->DofToVDof(i, 1));
        }
    }

    // Find unique x and y coordinates to determine grid structure
    std::vector<double> ux(xcoords), uy(ycoords);
    std::sort(ux.begin(), ux.end());
    std::sort(uy.begin(), uy.end());
    ux.erase(std::unique(ux.begin(), ux.end(),
        [](double a, double b){ return std::abs(a - b) < 1e-12; }), ux.end());
    uy.erase(std::unique(uy.begin(), uy.end(),
        [](double a, double b){ return std::abs(a - b) < 1e-12; }), uy.end());

    int nx_nodes = static_cast<int>(ux.size()); // number of nodes in x
    int ny_nodes = static_cast<int>(uy.size()); // number of nodes in y

    if (nx_nodes * ny_nodes != nnodes)
    {
        MFEM_ABORT("FD: Grid is not structured. "
                   "Expected " << nx_nodes << " x " << ny_nodes << " = "
                   << nx_nodes * ny_nodes << " nodes, but got " << nnodes);
    }

    double hx = (ux.back() - ux.front()) / (nx_nodes - 1);
    double hy = (uy.back() - uy.front()) / (ny_nodes - 1);

    std::cout << "FD: Grid dimensions: " << nx_nodes << " x " << ny_nodes << std::endl;
    std::cout << "FD: Grid spacing: hx = " << hx << ", hy = " << hy << std::endl;

    // Build a DOF index map: given (ix, iy) grid indices, find the DOF id.
    // Also build inverse map: given DOF id, find (ix, iy).
    std::vector<int> dof_ix(nnodes), dof_iy(nnodes);
    std::vector<std::vector<int>> grid_to_dof(nx_nodes, std::vector<int>(ny_nodes, -1));

    for (int i = 0; i < nnodes; ++i)
    {
        // Find grid indices by matching coordinates
        auto itx = std::lower_bound(ux.begin(), ux.end(), xcoords[i] - 1e-12);
        auto ity = std::lower_bound(uy.begin(), uy.end(), ycoords[i] - 1e-12);
        int ix = static_cast<int>(itx - ux.begin());
        int iy = static_cast<int>(ity - uy.begin());

        // Clamp in case of floating point issues
        ix = std::min(ix, nx_nodes - 1);
        iy = std::min(iy, ny_nodes - 1);

        dof_ix[i] = ix;
        dof_iy[i] = iy;
        grid_to_dof[ix][iy] = i;
    }

    // Verify all grid positions are filled
    for (int ix = 0; ix < nx_nodes; ++ix)
        for (int iy = 0; iy < ny_nodes; ++iy)
            if (grid_to_dof[ix][iy] < 0)
                MFEM_ABORT("FD: Missing node at grid position (" << ix << ", " << iy << ")");

    // Initialize derivative matrices
    this->sh_func_dx_  = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dy_  = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxx_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dyy_ = mfem::SparseMatrix(nnodes, nnodes);
    this->sh_func_dxy_ = mfem::SparseMatrix(nnodes, nnodes);

    for (int node_id = 0; node_id < nnodes; ++node_id)
    {
        int ix = dof_ix[node_id];
        int iy = dof_iy[node_id];

        // --- dx: d/dx ---
        if (ix > 0 && ix < nx_nodes - 1)
        {
            // Central difference: (f(i+1,j) - f(i-1,j)) / (2*hx)
            sh_func_dx_.Add(node_id, grid_to_dof[ix+1][iy],  1.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, grid_to_dof[ix-1][iy], -1.0 / (2.0 * hx));
        }
        else if (ix == 0)
        {
            // Forward difference: (-3f(i,j) + 4f(i+1,j) - f(i+2,j)) / (2*hx)
            sh_func_dx_.Add(node_id, grid_to_dof[ix  ][iy], -3.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, grid_to_dof[ix+1][iy],  4.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, grid_to_dof[ix+2][iy], -1.0 / (2.0 * hx));
        }
        else // ix == nx_nodes - 1
        {
            // Backward difference: (3f(i,j) - 4f(i-1,j) + f(i-2,j)) / (2*hx)
            sh_func_dx_.Add(node_id, grid_to_dof[ix  ][iy],  3.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, grid_to_dof[ix-1][iy], -4.0 / (2.0 * hx));
            sh_func_dx_.Add(node_id, grid_to_dof[ix-2][iy],  1.0 / (2.0 * hx));
        }

        // --- dy: d/dy ---
        if (iy > 0 && iy < ny_nodes - 1)
        {
            // Central difference: (f(i,j+1) - f(i,j-1)) / (2*hy)
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy+1],  1.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy-1], -1.0 / (2.0 * hy));
        }
        else if (iy == 0)
        {
            // Forward difference
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy  ], -3.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy+1],  4.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy+2], -1.0 / (2.0 * hy));
        }
        else // iy == ny_nodes - 1
        {
            // Backward difference
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy  ],  3.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy-1], -4.0 / (2.0 * hy));
            sh_func_dy_.Add(node_id, grid_to_dof[ix][iy-2],  1.0 / (2.0 * hy));
        }

        // --- dxx: d²/dx² ---
        if (ix > 0 && ix < nx_nodes - 1)
        {
            // Central: (f(i+1,j) - 2f(i,j) + f(i-1,j)) / hx²
            sh_func_dxx_.Add(node_id, grid_to_dof[ix+1][iy],  1.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix  ][iy], -2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix-1][iy],  1.0 / (hx * hx));
        }
        else if (ix == 0)
        {
            // Forward: (2f(i,j) - 5f(i+1,j) + 4f(i+2,j) - f(i+3,j)) / hx²
            sh_func_dxx_.Add(node_id, grid_to_dof[ix  ][iy],  2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix+1][iy], -5.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix+2][iy],  4.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix+3][iy], -1.0 / (hx * hx));
        }
        else // ix == nx_nodes - 1
        {
            // Backward: (2f(i,j) - 5f(i-1,j) + 4f(i-2,j) - f(i-3,j)) / hx²
            sh_func_dxx_.Add(node_id, grid_to_dof[ix  ][iy],  2.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix-1][iy], -5.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix-2][iy],  4.0 / (hx * hx));
            sh_func_dxx_.Add(node_id, grid_to_dof[ix-3][iy], -1.0 / (hx * hx));
        }

        // --- dyy: d²/dy² ---
        if (iy > 0 && iy < ny_nodes - 1)
        {
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy+1],  1.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy  ], -2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy-1],  1.0 / (hy * hy));
        }
        else if (iy == 0)
        {
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy  ],  2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy+1], -5.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy+2],  4.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy+3], -1.0 / (hy * hy));
        }
        else // iy == ny_nodes - 1
        {
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy  ],  2.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy-1], -5.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy-2],  4.0 / (hy * hy));
            sh_func_dyy_.Add(node_id, grid_to_dof[ix][iy-3], -1.0 / (hy * hy));
        }

        // --- dxy: d²/dxdy ---
        // Use central differences in both directions where possible.
        // (f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1)) / (4*hx*hy)
        // At boundaries, fall back to one-sided differences.
        int ixm = (ix > 0) ? ix - 1 : ix;
        int ixp = (ix < nx_nodes - 1) ? ix + 1 : ix;
        int iym = (iy > 0) ? iy - 1 : iy;
        int iyp = (iy < ny_nodes - 1) ? iy + 1 : iy;

        double dx_step = (ixp - ixm) * hx;
        double dy_step = (iyp - iym) * hy;

        if (dx_step > 0 && dy_step > 0)
        {
            double scale = 1.0 / (dx_step * dy_step);
            sh_func_dxy_.Add(node_id, grid_to_dof[ixp][iyp],  scale);
            sh_func_dxy_.Add(node_id, grid_to_dof[ixp][iym], -scale);
            sh_func_dxy_.Add(node_id, grid_to_dof[ixm][iyp], -scale);
            sh_func_dxy_.Add(node_id, grid_to_dof[ixm][iym],  scale);
        }
    }

    // Finalize sparse matrices
    this->sh_func_dx_.Finalize();
    this->sh_func_dy_.Finalize();
    this->sh_func_dxx_.Finalize();
    this->sh_func_dyy_.Finalize();
    this->sh_func_dxy_.Finalize();

    std::cout << "FD: Execution time for 2D finite difference matrices: "
              << timer.RealTime() << " s" << std::endl;
}


void FiniteDiff2d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    mfem::SparseMatrix derivative;

    if      (deriv == "dx")  { derivative = this->sh_func_dx_; }
    else if (deriv == "dy")  { derivative = this->sh_func_dy_; }
    else if (deriv == "dxx") { derivative = this->sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = this->sh_func_dyy_; }
    else if (deriv == "dxy") { derivative = this->sh_func_dxy_; }

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

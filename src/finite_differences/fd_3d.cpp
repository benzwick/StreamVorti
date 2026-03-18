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
#include <vector>

namespace StreamVorti {

// Helper: apply 1D 2nd-order first derivative stencil
static void Apply1D_D1_O2(mfem::SparseMatrix &mat, int row, int col_m, int col_0, int col_p,
                           int idx, int n, double h)
{
    if (idx > 0 && idx < n - 1)
    {
        mat.Add(row, col_p,  1.0 / (2.0 * h));
        mat.Add(row, col_m, -1.0 / (2.0 * h));
    }
    else if (idx == 0)
    {
        // Forward: (-3, 4, -1) / 2h — need col_0, col_p, col_p2
        // Caller must supply col_p2 separately; this helper is simplified.
        // We handle it inline in Update() for 3D.
        mat.Add(row, col_0, -3.0 / (2.0 * h));
        mat.Add(row, col_p,  4.0 / (2.0 * h));
        // col_p+1 handled by caller
    }
    // Boundary cases handled inline in Update()
}

void FiniteDiff3d::Update()
{
    mfem::StopWatch timer;
    timer.Start();
    std::cout << "FD: update 3D derivative matrices (order " << stencil_order_ << ")" << std::endl;

    if (stencil_order_ != 2)
    {
        MFEM_ABORT("FD: 4th-order stencils for 3D are not yet implemented. "
                   "Use stencil_order=2 for 3D finite differences.");
    }

    // --- Extract node coordinates from GridFunction ---
    std::vector<std::vector<double>> coords;
    this->ExtractNodeCoordinates(coords);
    const std::vector<double> &xc = coords[0];
    const std::vector<double> &yc = coords[1];
    const std::vector<double> &zc = coords[2];

    // --- Find unique coordinates and verify structured grid ---
    std::vector<double> ux, uy, uz;
    FindUniqueCoordinates(xc, ux);
    FindUniqueCoordinates(yc, uy);
    FindUniqueCoordinates(zc, uz);

    int nxn = static_cast<int>(ux.size());
    int nyn = static_cast<int>(uy.size());
    int nzn = static_cast<int>(uz.size());

    if (nxn * nyn * nzn != nnodes_)
    {
        MFEM_ABORT("FD: Grid is not structured rectangular. "
                   "Detected " << nxn << " x " << nyn << " x " << nzn << " = "
                   << nxn * nyn * nzn << " grid points, but the mesh has "
                   << nnodes_ << " nodes. "
                   "Finite differences require a structured Cartesian grid.");
    }

    if (nxn < 4 || nyn < 4 || nzn < 4)
    {
        MFEM_ABORT("FD: Grid too small for 2nd-order stencils. Need at least "
                   "4 nodes in each direction, got "
                   << nxn << " x " << nyn << " x " << nzn << ".");
    }

    // --- Verify uniform spacing ---
    double hx = VerifyUniformSpacing(ux, "x");
    double hy = VerifyUniformSpacing(uy, "y");
    double hz = VerifyUniformSpacing(uz, "z");

    std::cout << "FD: Grid: " << nxn << " x " << nyn << " x " << nzn
              << ", hx=" << hx << ", hy=" << hy << ", hz=" << hz << std::endl;

    // --- Build 3D index map ---
    double xtol = 0.5 * hx, ytol = 0.5 * hy, ztol = 0.5 * hz;

    std::vector<int> dof_ix(nnodes_), dof_iy(nnodes_), dof_iz(nnodes_);
    std::vector<int> g2d(nxn * nyn * nzn, -1);

    auto grid_idx = [&](int ix, int iy, int iz) -> int {
        return ix * nyn * nzn + iy * nzn + iz;
    };

    for (int i = 0; i < nnodes_; ++i)
    {
        int gix = CoordinateToIndex(xc[i], ux, xtol);
        int giy = CoordinateToIndex(yc[i], uy, ytol);
        int giz = CoordinateToIndex(zc[i], uz, ztol);

        if (gix < 0 || giy < 0 || giz < 0)
        {
            MFEM_ABORT("FD: Node " << i << " at (" << xc[i] << ", " << yc[i]
                       << ", " << zc[i] << ") could not be mapped to the "
                       "structured grid.");
        }

        dof_ix[i] = gix;
        dof_iy[i] = giy;
        dof_iz[i] = giz;

        int flat = grid_idx(gix, giy, giz);
        if (g2d[flat] >= 0)
        {
            MFEM_ABORT("FD: Duplicate node at grid position ("
                       << gix << ", " << giy << ", " << giz << ").");
        }
        g2d[flat] = i;
    }

    for (int idx = 0; idx < nxn * nyn * nzn; ++idx)
        if (g2d[idx] < 0)
            MFEM_ABORT("FD: Missing node in 3D grid at flat index " << idx);

    auto dof = [&](int ix, int iy, int iz) -> int {
        return g2d[grid_idx(ix, iy, iz)];
    };

    // --- Initialize matrices ---
    sh_func_dx_  = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dy_  = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dz_  = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dxx_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dyy_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dzz_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dxy_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dxz_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dyz_ = mfem::SparseMatrix(nnodes_, nnodes_);

    // --- Macro for 2nd-order 1D stencils along arbitrary axis ---
    // D1: first derivative, D2: second derivative
    #define FD_D1_O2(mat, nid, idx, n, h, DOF_FUNC) \
    do { \
        if ((idx) > 0 && (idx) < (n) - 1) { \
            (mat).Add((nid), DOF_FUNC((idx)+1),  1.0 / (2.0 * (h))); \
            (mat).Add((nid), DOF_FUNC((idx)-1), -1.0 / (2.0 * (h))); \
        } else if ((idx) == 0) { \
            (mat).Add((nid), DOF_FUNC(0), -3.0 / (2.0 * (h))); \
            (mat).Add((nid), DOF_FUNC(1),  4.0 / (2.0 * (h))); \
            (mat).Add((nid), DOF_FUNC(2), -1.0 / (2.0 * (h))); \
        } else { \
            (mat).Add((nid), DOF_FUNC((n)-1),  3.0 / (2.0 * (h))); \
            (mat).Add((nid), DOF_FUNC((n)-2), -4.0 / (2.0 * (h))); \
            (mat).Add((nid), DOF_FUNC((n)-3),  1.0 / (2.0 * (h))); \
        } \
    } while(0)

    #define FD_D2_O2(mat, nid, idx, n, h, DOF_FUNC) \
    do { \
        double _h2 = (h) * (h); \
        if ((idx) > 0 && (idx) < (n) - 1) { \
            (mat).Add((nid), DOF_FUNC((idx)+1),  1.0 / _h2); \
            (mat).Add((nid), DOF_FUNC((idx)  ), -2.0 / _h2); \
            (mat).Add((nid), DOF_FUNC((idx)-1),  1.0 / _h2); \
        } else if ((idx) == 0) { \
            (mat).Add((nid), DOF_FUNC(0),  2.0 / _h2); \
            (mat).Add((nid), DOF_FUNC(1), -5.0 / _h2); \
            (mat).Add((nid), DOF_FUNC(2),  4.0 / _h2); \
            (mat).Add((nid), DOF_FUNC(3), -1.0 / _h2); \
        } else { \
            (mat).Add((nid), DOF_FUNC((n)-1),  2.0 / _h2); \
            (mat).Add((nid), DOF_FUNC((n)-2), -5.0 / _h2); \
            (mat).Add((nid), DOF_FUNC((n)-3),  4.0 / _h2); \
            (mat).Add((nid), DOF_FUNC((n)-4), -1.0 / _h2); \
        } \
    } while(0)

    for (int nid = 0; nid < nnodes_; ++nid)
    {
        int ix = dof_ix[nid];
        int iy = dof_iy[nid];
        int iz = dof_iz[nid];

        // dx, dxx
        auto dof_x = [&](int i) { return dof(i, iy, iz); };
        FD_D1_O2(sh_func_dx_, nid, ix, nxn, hx, dof_x);
        FD_D2_O2(sh_func_dxx_, nid, ix, nxn, hx, dof_x);

        // dy, dyy
        auto dof_y = [&](int j) { return dof(ix, j, iz); };
        FD_D1_O2(sh_func_dy_, nid, iy, nyn, hy, dof_y);
        FD_D2_O2(sh_func_dyy_, nid, iy, nyn, hy, dof_y);

        // dz, dzz
        auto dof_z = [&](int k) { return dof(ix, iy, k); };
        FD_D1_O2(sh_func_dz_, nid, iz, nzn, hz, dof_z);
        FD_D2_O2(sh_func_dzz_, nid, iz, nzn, hz, dof_z);

        // Mixed derivatives: dxy, dxz, dyz
        // Use product of central 1st-derivative stencils with one-sided fallback
        auto add_mixed = [&](mfem::SparseMatrix &mat,
                             int i1, int n1, double h1,
                             int i2, int n2, double h2,
                             auto dof_func)
        {
            int i1m = (i1 > 0) ? i1 - 1 : i1;
            int i1p = (i1 < n1 - 1) ? i1 + 1 : i1;
            int i2m = (i2 > 0) ? i2 - 1 : i2;
            int i2p = (i2 < n2 - 1) ? i2 + 1 : i2;
            double ds1 = (i1p - i1m) * h1;
            double ds2 = (i2p - i2m) * h2;
            if (ds1 > 0 && ds2 > 0)
            {
                double s = 1.0 / (ds1 * ds2);
                mat.Add(nid, dof_func(i1p, i2p),  s);
                mat.Add(nid, dof_func(i1p, i2m), -s);
                mat.Add(nid, dof_func(i1m, i2p), -s);
                mat.Add(nid, dof_func(i1m, i2m),  s);
            }
        };

        add_mixed(sh_func_dxy_, ix, nxn, hx, iy, nyn, hy,
                  [&](int a, int b){ return dof(a, b, iz); });
        add_mixed(sh_func_dxz_, ix, nxn, hx, iz, nzn, hz,
                  [&](int a, int c){ return dof(a, iy, c); });
        add_mixed(sh_func_dyz_, iy, nyn, hy, iz, nzn, hz,
                  [&](int b, int c){ return dof(ix, b, c); });
    }

    #undef FD_D1_O2
    #undef FD_D2_O2

    // Finalize
    sh_func_dx_.Finalize();
    sh_func_dy_.Finalize();
    sh_func_dz_.Finalize();
    sh_func_dxx_.Finalize();
    sh_func_dyy_.Finalize();
    sh_func_dzz_.Finalize();
    sh_func_dxy_.Finalize();
    sh_func_dxz_.Finalize();
    sh_func_dyz_.Finalize();

    std::cout << "FD: Execution time for 3D FD matrices: "
              << timer.RealTime() << " s" << std::endl;
}


void FiniteDiff3d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    const mfem::SparseMatrix *derivative = nullptr;

    if      (deriv == "dx")  { derivative = &sh_func_dx_; }
    else if (deriv == "dy")  { derivative = &sh_func_dy_; }
    else if (deriv == "dz")  { derivative = &sh_func_dz_; }
    else if (deriv == "dxx") { derivative = &sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = &sh_func_dyy_; }
    else if (deriv == "dzz") { derivative = &sh_func_dzz_; }
    else if (deriv == "dxy") { derivative = &sh_func_dxy_; }
    else if (deriv == "dxz") { derivative = &sh_func_dxz_; }
    else if (deriv == "dyz") { derivative = &sh_func_dyz_; }
    else
    {
        MFEM_ABORT("FiniteDiff3d::SaveDerivToFile: Unknown derivative '"
                   << deriv << "'.");
    }

    if (derivative->Height() == 0)
    {
        MFEM_ABORT("Could not save FD derivative '" << deriv
                   << "'. Derivatives have not been computed.");
    }

    std::string path = "";
    std::size_t last_slash = filename.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }
    if (!path.empty()) {
        std::filesystem::create_directories(path);
    }

    std::string out_filename = filename;
    if (filename.find_last_of(".") == std::string::npos ||
        filename.substr(filename.find_last_of(".")) != ".txt") {
        out_filename = filename + ".txt";
    }

    std::ofstream out(out_filename, std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
        MFEM_ABORT("Could not open file for writing: " << out_filename);
    }
    derivative->PrintMatlab(out);
    out.close();
}

} // namespace StreamVorti

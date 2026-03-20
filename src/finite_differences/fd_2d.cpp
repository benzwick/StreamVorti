/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
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
#include <vector>

namespace StreamVorti {

// ====================================================================
// 2nd-order stencil helpers
// ====================================================================

// d/dx, central: (-1, 0, 1) / 2h
// d/dx, forward: (-3, 4, -1) / 2h
// d/dx, backward: (3, -4, 1) / 2h
static void AddDx_O2(mfem::SparseMatrix &mat, int row,
                      const std::vector<std::vector<int>> &g2d,
                      int ix, int iy, int nx, double hx)
{
    if (ix > 0 && ix < nx - 1)
    {
        mat.Add(row, g2d[ix+1][iy],  1.0 / (2.0 * hx));
        mat.Add(row, g2d[ix-1][iy], -1.0 / (2.0 * hx));
    }
    else if (ix == 0)
    {
        mat.Add(row, g2d[ix  ][iy], -3.0 / (2.0 * hx));
        mat.Add(row, g2d[ix+1][iy],  4.0 / (2.0 * hx));
        mat.Add(row, g2d[ix+2][iy], -1.0 / (2.0 * hx));
    }
    else
    {
        mat.Add(row, g2d[ix  ][iy],  3.0 / (2.0 * hx));
        mat.Add(row, g2d[ix-1][iy], -4.0 / (2.0 * hx));
        mat.Add(row, g2d[ix-2][iy],  1.0 / (2.0 * hx));
    }
}

// d²/dx², central: (1, -2, 1) / h²
// d²/dx², forward: (2, -5, 4, -1) / h²
// d²/dx², backward: (2, -5, 4, -1) / h²  (mirror)
static void AddDxx_O2(mfem::SparseMatrix &mat, int row,
                       const std::vector<std::vector<int>> &g2d,
                       int ix, int iy, int nx, double hx)
{
    double h2 = hx * hx;
    if (ix > 0 && ix < nx - 1)
    {
        mat.Add(row, g2d[ix+1][iy],  1.0 / h2);
        mat.Add(row, g2d[ix  ][iy], -2.0 / h2);
        mat.Add(row, g2d[ix-1][iy],  1.0 / h2);
    }
    else if (ix == 0)
    {
        mat.Add(row, g2d[ix  ][iy],  2.0 / h2);
        mat.Add(row, g2d[ix+1][iy], -5.0 / h2);
        mat.Add(row, g2d[ix+2][iy],  4.0 / h2);
        mat.Add(row, g2d[ix+3][iy], -1.0 / h2);
    }
    else
    {
        mat.Add(row, g2d[ix  ][iy],  2.0 / h2);
        mat.Add(row, g2d[ix-1][iy], -5.0 / h2);
        mat.Add(row, g2d[ix-2][iy],  4.0 / h2);
        mat.Add(row, g2d[ix-3][iy], -1.0 / h2);
    }
}

// ====================================================================
// 4th-order stencil helpers
// ====================================================================

// d/dx, central 5-point: (1, -8, 0, 8, -1) / 12h
// d/dx, boundary: use progressively one-sided stencils
static void AddDx_O4(mfem::SparseMatrix &mat, int row,
                      const std::vector<std::vector<int>> &g2d,
                      int ix, int iy, int nx, double hx)
{
    double c = 1.0 / (12.0 * hx);
    if (ix >= 2 && ix <= nx - 3)
    {
        // Central 5-point
        mat.Add(row, g2d[ix-2][iy],  1.0 * c);
        mat.Add(row, g2d[ix-1][iy], -8.0 * c);
        mat.Add(row, g2d[ix+1][iy],  8.0 * c);
        mat.Add(row, g2d[ix+2][iy], -1.0 * c);
    }
    else if (ix == 0)
    {
        // Forward 5-point: (-25, 48, -36, 16, -3) / 12h
        mat.Add(row, g2d[0][iy], -25.0 * c);
        mat.Add(row, g2d[1][iy],  48.0 * c);
        mat.Add(row, g2d[2][iy], -36.0 * c);
        mat.Add(row, g2d[3][iy],  16.0 * c);
        mat.Add(row, g2d[4][iy],  -3.0 * c);
    }
    else if (ix == 1)
    {
        // Near-boundary: use forward-biased 5-point
        // (-3, -10, 18, -6, 1) / 12h
        mat.Add(row, g2d[0][iy], -3.0 * c);
        mat.Add(row, g2d[1][iy],-10.0 * c);
        mat.Add(row, g2d[2][iy], 18.0 * c);
        mat.Add(row, g2d[3][iy], -6.0 * c);
        mat.Add(row, g2d[4][iy],  1.0 * c);
    }
    else if (ix == nx - 1)
    {
        // Backward 5-point: (25, -48, 36, -16, 3) / 12h
        mat.Add(row, g2d[nx-1][iy],  25.0 * c);
        mat.Add(row, g2d[nx-2][iy], -48.0 * c);
        mat.Add(row, g2d[nx-3][iy],  36.0 * c);
        mat.Add(row, g2d[nx-4][iy], -16.0 * c);
        mat.Add(row, g2d[nx-5][iy],   3.0 * c);
    }
    else // ix == nx - 2
    {
        // Near-boundary backward-biased
        mat.Add(row, g2d[nx-1][iy],  3.0 * c);
        mat.Add(row, g2d[nx-2][iy], 10.0 * c);
        mat.Add(row, g2d[nx-3][iy],-18.0 * c);
        mat.Add(row, g2d[nx-4][iy],  6.0 * c);
        mat.Add(row, g2d[nx-5][iy], -1.0 * c);
    }
}

// d²/dx², central 5-point: (-1, 16, -30, 16, -1) / 12h²
static void AddDxx_O4(mfem::SparseMatrix &mat, int row,
                       const std::vector<std::vector<int>> &g2d,
                       int ix, int iy, int nx, double hx)
{
    double c = 1.0 / (12.0 * hx * hx);
    if (ix >= 2 && ix <= nx - 3)
    {
        mat.Add(row, g2d[ix-2][iy],  -1.0 * c);
        mat.Add(row, g2d[ix-1][iy],  16.0 * c);
        mat.Add(row, g2d[ix  ][iy], -30.0 * c);
        mat.Add(row, g2d[ix+1][iy],  16.0 * c);
        mat.Add(row, g2d[ix+2][iy],  -1.0 * c);
    }
    else if (ix == 0)
    {
        // Forward 6-point: (45, -154, 214, -156, 61, -10) / 12h²
        mat.Add(row, g2d[0][iy],  45.0 * c);
        mat.Add(row, g2d[1][iy],-154.0 * c);
        mat.Add(row, g2d[2][iy], 214.0 * c);
        mat.Add(row, g2d[3][iy],-156.0 * c);
        mat.Add(row, g2d[4][iy],  61.0 * c);
        mat.Add(row, g2d[5][iy], -10.0 * c);
    }
    else if (ix == 1)
    {
        // Near-boundary: (10, -15, -4, 14, -6, 1) / 12h²
        mat.Add(row, g2d[0][iy],  10.0 * c);
        mat.Add(row, g2d[1][iy], -15.0 * c);
        mat.Add(row, g2d[2][iy],  -4.0 * c);
        mat.Add(row, g2d[3][iy],  14.0 * c);
        mat.Add(row, g2d[4][iy],  -6.0 * c);
        mat.Add(row, g2d[5][iy],   1.0 * c);
    }
    else if (ix == nx - 1)
    {
        mat.Add(row, g2d[nx-1][iy],  45.0 * c);
        mat.Add(row, g2d[nx-2][iy],-154.0 * c);
        mat.Add(row, g2d[nx-3][iy], 214.0 * c);
        mat.Add(row, g2d[nx-4][iy],-156.0 * c);
        mat.Add(row, g2d[nx-5][iy],  61.0 * c);
        mat.Add(row, g2d[nx-6][iy], -10.0 * c);
    }
    else // ix == nx - 2
    {
        mat.Add(row, g2d[nx-1][iy],  10.0 * c);
        mat.Add(row, g2d[nx-2][iy], -15.0 * c);
        mat.Add(row, g2d[nx-3][iy],  -4.0 * c);
        mat.Add(row, g2d[nx-4][iy],  14.0 * c);
        mat.Add(row, g2d[nx-5][iy],  -6.0 * c);
        mat.Add(row, g2d[nx-6][iy],   1.0 * c);
    }
}

// ====================================================================
// Update
// ====================================================================

void FiniteDiff2d::Update()
{
    mfem::StopWatch timer;
    timer.Start();
    std::cout << "FD: update 2D derivative matrices (order " << stencil_order_ << ")" << std::endl;

    // --- Extract node coordinates from GridFunction (not mesh vertices) ---
    std::vector<std::vector<double>> coords;
    this->ExtractNodeCoordinates(coords);
    const std::vector<double> &xcoords = coords[0];
    const std::vector<double> &ycoords = coords[1];

    // --- Find unique coordinates and verify structured grid ---
    std::vector<double> ux, uy;
    FindUniqueCoordinates(xcoords, ux);
    FindUniqueCoordinates(ycoords, uy);

    int nx = static_cast<int>(ux.size());
    int ny = static_cast<int>(uy.size());

    if (nx * ny != nnodes_)
    {
        MFEM_ABORT("FD: Grid is not structured rectangular. "
                   "Detected " << nx << " x " << ny << " = "
                   << nx * ny << " grid points, but the mesh has "
                   << nnodes_ << " nodes. "
                   "Finite differences require a structured Cartesian grid. "
                   "Use DCPSE for unstructured meshes.");
    }

    // Minimum grid size for stencil order
    int min_nodes = (stencil_order_ == 4) ? 6 : 4;
    if (nx < min_nodes || ny < min_nodes)
    {
        MFEM_ABORT("FD: Grid too small for order-" << stencil_order_
                   << " stencils. Need at least " << min_nodes
                   << " nodes in each direction, got " << nx << " x " << ny << ".");
    }

    // --- Verify uniform spacing ---
    double hx = VerifyUniformSpacing(ux, "x");
    double hy = VerifyUniformSpacing(uy, "y");

    std::cout << "FD: Grid: " << nx << " x " << ny
              << ", hx = " << hx << ", hy = " << hy << std::endl;

    // --- Build DOF index map ---
    double xtol = 0.5 * hx;
    double ytol = 0.5 * hy;

    std::vector<int> dof_ix(nnodes_), dof_iy(nnodes_);
    std::vector<std::vector<int>> g2d(nx, std::vector<int>(ny, -1));

    for (int i = 0; i < nnodes_; ++i)
    {
        int gix = CoordinateToIndex(xcoords[i], ux, xtol);
        int giy = CoordinateToIndex(ycoords[i], uy, ytol);

        if (gix < 0 || giy < 0)
        {
            MFEM_ABORT("FD: Node " << i << " at (" << xcoords[i] << ", "
                       << ycoords[i] << ") could not be mapped to the "
                       "structured grid. This indicates the mesh is not "
                       "a proper Cartesian grid.");
        }

        dof_ix[i] = gix;
        dof_iy[i] = giy;

        if (g2d[gix][giy] >= 0)
        {
            MFEM_ABORT("FD: Duplicate node detected at grid position ("
                       << gix << ", " << giy << "). Nodes " << g2d[gix][giy]
                       << " and " << i << " map to the same grid point.");
        }
        g2d[gix][giy] = i;
    }

    // Verify all grid positions are filled
    for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy)
            if (g2d[ix][iy] < 0)
                MFEM_ABORT("FD: Missing node at grid position ("
                           << ix << ", " << iy << "). The mesh does not "
                           "cover the full rectangular grid.");

    // --- Initialize derivative matrices ---
    sh_func_dx_  = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dy_  = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dxx_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dyy_ = mfem::SparseMatrix(nnodes_, nnodes_);
    sh_func_dxy_ = mfem::SparseMatrix(nnodes_, nnodes_);

    // --- Assemble stencils ---
    for (int nid = 0; nid < nnodes_; ++nid)
    {
        int ix = dof_ix[nid];
        int iy = dof_iy[nid];

        // Transpose g2d for y-direction stencils: swap ix/iy roles
        // We create a lambda to access g2d with swapped indices
        // g2d_yx[iy][ix] == g2d[ix][iy]

        if (stencil_order_ == 2)
        {
            AddDx_O2(sh_func_dx_, nid, g2d, ix, iy, nx, hx);
            AddDxx_O2(sh_func_dxx_, nid, g2d, ix, iy, nx, hx);
        }
        else
        {
            AddDx_O4(sh_func_dx_, nid, g2d, ix, iy, nx, hx);
            AddDxx_O4(sh_func_dxx_, nid, g2d, ix, iy, nx, hx);
        }

        // dy: same stencils but along y-axis.
        // Build a transposed view for the y-direction helpers.
        // Since the helpers take g2d[index_along_axis][fixed_index],
        // we need g2d_y where g2d_y[iy_idx][ix_idx] == g2d[ix_idx][iy_idx].
        // Instead of building a transposed array, apply the stencils directly.
        if (stencil_order_ == 2)
        {
            // dy
            if (iy > 0 && iy < ny - 1)
            {
                sh_func_dy_.Add(nid, g2d[ix][iy+1],  1.0 / (2.0 * hy));
                sh_func_dy_.Add(nid, g2d[ix][iy-1], -1.0 / (2.0 * hy));
            }
            else if (iy == 0)
            {
                sh_func_dy_.Add(nid, g2d[ix][iy  ], -3.0 / (2.0 * hy));
                sh_func_dy_.Add(nid, g2d[ix][iy+1],  4.0 / (2.0 * hy));
                sh_func_dy_.Add(nid, g2d[ix][iy+2], -1.0 / (2.0 * hy));
            }
            else
            {
                sh_func_dy_.Add(nid, g2d[ix][iy  ],  3.0 / (2.0 * hy));
                sh_func_dy_.Add(nid, g2d[ix][iy-1], -4.0 / (2.0 * hy));
                sh_func_dy_.Add(nid, g2d[ix][iy-2],  1.0 / (2.0 * hy));
            }

            // dyy
            double h2y = hy * hy;
            if (iy > 0 && iy < ny - 1)
            {
                sh_func_dyy_.Add(nid, g2d[ix][iy+1],  1.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy  ], -2.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy-1],  1.0 / h2y);
            }
            else if (iy == 0)
            {
                sh_func_dyy_.Add(nid, g2d[ix][iy  ],  2.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy+1], -5.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy+2],  4.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy+3], -1.0 / h2y);
            }
            else
            {
                sh_func_dyy_.Add(nid, g2d[ix][iy  ],  2.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy-1], -5.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy-2],  4.0 / h2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy-3], -1.0 / h2y);
            }
        }
        else // 4th order y-stencils
        {
            double cy = 1.0 / (12.0 * hy);
            // dy
            if (iy >= 2 && iy <= ny - 3)
            {
                sh_func_dy_.Add(nid, g2d[ix][iy-2],  1.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][iy-1], -8.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][iy+1],  8.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][iy+2], -1.0 * cy);
            }
            else if (iy == 0)
            {
                sh_func_dy_.Add(nid, g2d[ix][0], -25.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][1],  48.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][2], -36.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][3],  16.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][4],  -3.0 * cy);
            }
            else if (iy == 1)
            {
                sh_func_dy_.Add(nid, g2d[ix][0], -3.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][1],-10.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][2], 18.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][3], -6.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][4],  1.0 * cy);
            }
            else if (iy == ny - 1)
            {
                sh_func_dy_.Add(nid, g2d[ix][ny-1],  25.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-2], -48.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-3],  36.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-4], -16.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-5],   3.0 * cy);
            }
            else // iy == ny - 2
            {
                sh_func_dy_.Add(nid, g2d[ix][ny-1],  3.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-2], 10.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-3],-18.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-4],  6.0 * cy);
                sh_func_dy_.Add(nid, g2d[ix][ny-5], -1.0 * cy);
            }

            // dyy
            double c2y = 1.0 / (12.0 * hy * hy);
            if (iy >= 2 && iy <= ny - 3)
            {
                sh_func_dyy_.Add(nid, g2d[ix][iy-2],  -1.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy-1],  16.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy  ], -30.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy+1],  16.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][iy+2],  -1.0 * c2y);
            }
            else if (iy == 0)
            {
                sh_func_dyy_.Add(nid, g2d[ix][0],  45.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][1],-154.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][2], 214.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][3],-156.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][4],  61.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][5], -10.0 * c2y);
            }
            else if (iy == 1)
            {
                sh_func_dyy_.Add(nid, g2d[ix][0],  10.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][1], -15.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][2],  -4.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][3],  14.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][4],  -6.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][5],   1.0 * c2y);
            }
            else if (iy == ny - 1)
            {
                sh_func_dyy_.Add(nid, g2d[ix][ny-1],  45.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-2],-154.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-3], 214.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-4],-156.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-5],  61.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-6], -10.0 * c2y);
            }
            else // iy == ny - 2
            {
                sh_func_dyy_.Add(nid, g2d[ix][ny-1],  10.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-2], -15.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-3],  -4.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-4],  14.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-5],  -6.0 * c2y);
                sh_func_dyy_.Add(nid, g2d[ix][ny-6],   1.0 * c2y);
            }
        }

        // --- dxy: d²/dxdy ---
        // Use the product of 1D first-derivative stencils.
        // For 2nd order: central difference in both directions at interior,
        // one-sided at boundaries.
        {
            int ixm = (ix > 0) ? ix - 1 : ix;
            int ixp = (ix < nx - 1) ? ix + 1 : ix;
            int iym = (iy > 0) ? iy - 1 : iy;
            int iyp = (iy < ny - 1) ? iy + 1 : iy;

            double dx_step = (ixp - ixm) * hx;
            double dy_step = (iyp - iym) * hy;

            if (dx_step > 0 && dy_step > 0)
            {
                double s = 1.0 / (dx_step * dy_step);
                sh_func_dxy_.Add(nid, g2d[ixp][iyp],  s);
                sh_func_dxy_.Add(nid, g2d[ixp][iym], -s);
                sh_func_dxy_.Add(nid, g2d[ixm][iyp], -s);
                sh_func_dxy_.Add(nid, g2d[ixm][iym],  s);
            }
        }
    }

    // --- Finalize ---
    sh_func_dx_.Finalize();
    sh_func_dy_.Finalize();
    sh_func_dxx_.Finalize();
    sh_func_dyy_.Finalize();
    sh_func_dxy_.Finalize();

    std::cout << "FD: Execution time for 2D FD matrices: "
              << timer.RealTime() << " s" << std::endl;
}


void FiniteDiff2d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    const mfem::SparseMatrix *derivative = nullptr;

    if      (deriv == "dx")  { derivative = &sh_func_dx_; }
    else if (deriv == "dy")  { derivative = &sh_func_dy_; }
    else if (deriv == "dxx") { derivative = &sh_func_dxx_; }
    else if (deriv == "dyy") { derivative = &sh_func_dyy_; }
    else if (deriv == "dxy") { derivative = &sh_func_dxy_; }
    else
    {
        MFEM_ABORT("FiniteDiff2d::SaveDerivToFile: Unknown derivative '"
                   << deriv << "'. Use dx, dy, dxx, dyy, or dxy.");
    }

    if (derivative->Height() == 0)
    {
        MFEM_ABORT("Could not save FD derivative '" << deriv
                   << "'. Derivatives have not been computed (call Update first).");
    }

    // Create output directory if needed
    std::string path = "";
    std::size_t last_slash = filename.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        path = filename.substr(0, last_slash);
    }
    if (!path.empty()) {
        std::filesystem::create_directories(path);
    }

    // Ensure .txt extension
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

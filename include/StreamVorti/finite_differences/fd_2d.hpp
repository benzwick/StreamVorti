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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_FD_2D_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_FD_2D_HPP_

#include "StreamVorti/finite_differences/fd.hpp"

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class FiniteDiff2d
 * \brief Finite difference derivatives on 2D structured grids.
 *
 * Builds sparse derivative matrices for a 2D Cartesian grid.
 *
 * Supported stencil orders (set via constructor):
 *
 * Order 2 (default):
 *   1st deriv interior: central 3-point     (-1, 0, 1) / 2h
 *   1st deriv boundary: one-sided 3-point   (-3, 4, -1) / 2h
 *   2nd deriv interior: central 3-point     (1, -2, 1) / h²
 *   2nd deriv boundary: one-sided 4-point   (2, -5, 4, -1) / h²
 *
 * Order 4:
 *   1st deriv interior: central 5-point     (1, -8, 0, 8, -1) / 12h
 *   1st deriv boundary: one-sided 5-point   (-25, 48, -36, 16, -3) / 12h
 *   2nd deriv interior: central 5-point     (-1, 16, -30, 16, -1) / 12h²
 *   2nd deriv boundary: one-sided 6-point   (45, -154, 214, -156, 61, -10) / 12h²
 */
class FiniteDiff2d: public FiniteDiff
{
public:
    /*!
     * \brief Construct 2D FD operator.
     * \param gf GridFunction on an H1 FE space over a structured mesh.
     * \param stencil_order Accuracy order: 2 or 4 (default 2).
     */
    FiniteDiff2d(mfem::GridFunction &gf, int stencil_order = 2)
        : FiniteDiff(gf, stencil_order) {}

    void Update();

    void SaveDerivToFile(const std::string &deriv, const std::string &filename) const;

    const mfem::SparseMatrix & D(int i) const {
        switch (i)
        {
        case 0: return this->sh_func_dx_;
        case 1: return this->sh_func_dy_;
        }
        MFEM_ABORT( "Index " << i << " out of bounds." );
    }

    inline const mfem::SparseMatrix & ShapeFunctionDx() const { return this->sh_func_dx_; }
    inline const mfem::SparseMatrix & ShapeFunctionDy() const { return this->sh_func_dy_; }
    inline const mfem::SparseMatrix & ShapeFunctionDxx() const { return this->sh_func_dxx_; }
    inline const mfem::SparseMatrix & ShapeFunctionDyy() const { return this->sh_func_dyy_; }
    inline const mfem::SparseMatrix & ShapeFunctionDxy() const { return this->sh_func_dxy_; }

    const mfem::SparseMatrix & Dxx() const override { return this->sh_func_dxx_; }
    const mfem::SparseMatrix & Dyy() const override { return this->sh_func_dyy_; }
    const mfem::SparseMatrix & Dxy() const override { return this->sh_func_dxy_; }

private:
    mfem::SparseMatrix sh_func_dx_;  /*!< 1st x derivative matrix */
    mfem::SparseMatrix sh_func_dy_;  /*!< 1st y derivative matrix */
    mfem::SparseMatrix sh_func_dxx_; /*!< 2nd xx derivative matrix */
    mfem::SparseMatrix sh_func_dyy_; /*!< 2nd yy derivative matrix */
    mfem::SparseMatrix sh_func_dxy_; /*!< 2nd xy derivative matrix */
};

} // namespace StreamVorti

#endif // STREAMVORTI_FINITE_DIFFERENCES_FD_2D_HPP_

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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_FD_3D_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_FD_3D_HPP_

#include "StreamVorti/finite_differences/fd.hpp"

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class FiniteDiff3d
 * \brief Finite difference derivatives on 3D structured grids.
 *
 * Same stencil options as FiniteDiff2d, extended to three dimensions.
 */
class FiniteDiff3d: public FiniteDiff
{
public:
    FiniteDiff3d(mfem::GridFunction &gf, int stencil_order = 2)
        : FiniteDiff(gf, stencil_order) {}

    void Update();

    void SaveDerivToFile(const std::string &deriv, const std::string &filename) const;

    const mfem::SparseMatrix & D(int i) const {
        switch (i)
        {
        case 0: return this->sh_func_dx_;
        case 1: return this->sh_func_dy_;
        case 2: return this->sh_func_dz_;
        }
        MFEM_ABORT( "Index " << i << " out of bounds." );
    }

    inline const mfem::SparseMatrix & ShapeFunctionDx() const { return this->sh_func_dx_; }
    inline const mfem::SparseMatrix & ShapeFunctionDy() const { return this->sh_func_dy_; }
    inline const mfem::SparseMatrix & ShapeFunctionDz() const { return this->sh_func_dz_; }
    inline const mfem::SparseMatrix & ShapeFunctionDxx() const { return this->sh_func_dxx_; }
    inline const mfem::SparseMatrix & ShapeFunctionDyy() const { return this->sh_func_dyy_; }
    inline const mfem::SparseMatrix & ShapeFunctionDzz() const { return this->sh_func_dzz_; }
    inline const mfem::SparseMatrix & ShapeFunctionDxy() const { return this->sh_func_dxy_; }
    inline const mfem::SparseMatrix & ShapeFunctionDxz() const { return this->sh_func_dxz_; }
    inline const mfem::SparseMatrix & ShapeFunctionDyz() const { return this->sh_func_dyz_; }

    const mfem::SparseMatrix & Dxx() const override { return this->sh_func_dxx_; }
    const mfem::SparseMatrix & Dyy() const override { return this->sh_func_dyy_; }
    const mfem::SparseMatrix & Dxy() const override { return this->sh_func_dxy_; }

private:
    mfem::SparseMatrix sh_func_dx_;
    mfem::SparseMatrix sh_func_dy_;
    mfem::SparseMatrix sh_func_dz_;
    mfem::SparseMatrix sh_func_dxx_;
    mfem::SparseMatrix sh_func_dyy_;
    mfem::SparseMatrix sh_func_dzz_;
    mfem::SparseMatrix sh_func_dxy_;
    mfem::SparseMatrix sh_func_dxz_;
    mfem::SparseMatrix sh_func_dyz_;
};

} // namespace StreamVorti

#endif // STREAMVORTI_FINITE_DIFFERENCES_FD_3D_HPP_

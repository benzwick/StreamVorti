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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_2D_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_2D_HPP_

#include "StreamVorti/finite_differences/par_fd.hpp"

#ifdef MFEM_USE_MPI

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class ParFiniteDiff2d
 * \brief Parallel 2D finite difference derivatives.
 *
 * Algorithm:
 * 1. Each rank detects the global structured grid from the parallel mesh.
 * 2. Each rank identifies which global grid nodes it owns (true DOFs).
 * 3. FD stencil coefficients are computed analytically.
 * 4. Entries referencing local DOFs go into the diagonal SparseMatrix block;
 *    entries referencing DOFs on other ranks go into the off-diagonal block.
 * 5. HypreParMatrix is constructed from the two blocks.
 */
class ParFiniteDiff2d : public ParFiniteDiff
{
public:
    ParFiniteDiff2d(mfem::ParGridFunction &gf, int stencil_order = 2);
    virtual ~ParFiniteDiff2d();

    void Update() override;

    void SaveDerivToFile(const std::string &deriv,
                         const std::string &filename) const override;

    const mfem::HypreParMatrix & D(int i) const override {
        switch (i)
        {
        case 0: return *sh_func_dx_;
        case 1: return *sh_func_dy_;
        }
        MFEM_ABORT("ParFiniteDiff2d::D: Index " << i << " out of bounds.");
    }

    inline const mfem::HypreParMatrix & ShapeFunctionDx() const { return *sh_func_dx_; }
    inline const mfem::HypreParMatrix & ShapeFunctionDy() const { return *sh_func_dy_; }
    inline const mfem::HypreParMatrix & ShapeFunctionDxx() const { return *sh_func_dxx_; }
    inline const mfem::HypreParMatrix & ShapeFunctionDyy() const { return *sh_func_dyy_; }
    inline const mfem::HypreParMatrix & ShapeFunctionDxy() const { return *sh_func_dxy_; }

    const mfem::HypreParMatrix & Dxx() const override { return *sh_func_dxx_; }
    const mfem::HypreParMatrix & Dyy() const override { return *sh_func_dyy_; }
    const mfem::HypreParMatrix & Dxy() const override { return *sh_func_dxy_; }

private:
    mfem::HypreParMatrix *sh_func_dx_;
    mfem::HypreParMatrix *sh_func_dy_;
    mfem::HypreParMatrix *sh_func_dxx_;
    mfem::HypreParMatrix *sh_func_dyy_;
    mfem::HypreParMatrix *sh_func_dxy_;
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_2D_HPP_

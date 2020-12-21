/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
 */

#ifndef STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_
#define STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_

#include "StreamVorti/approximants/dcpse.hpp"

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class Dcpse2d
 * \brief DC PSE derivatives in 2D.
 */
class Dcpse2d: public Dcpse
{
public:
    /*!
     * \brief Dcpse2d constructor to match nodes of an MFEM H1 GridFunction.
     */
    Dcpse2d(mfem::GridFunction &gf, int NumNeighbors)
        : Dcpse(gf, NumNeighbors) {}

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

private:
    mfem::SparseMatrix sh_func_dx_; /*!< The shape function 1st x derivative matrix. */
    mfem::SparseMatrix sh_func_dy_; /*!< The shape function 1st y derivative matrix. */
    mfem::SparseMatrix sh_func_dxx_; /*!< The shape function 2nd xx derivative matrix. */
    mfem::SparseMatrix sh_func_dyy_; /*!< The shape function 2nd yy derivative matrix. */
    mfem::SparseMatrix sh_func_dxy_; /*!< The shape function 2nd xy derivative matrix. */
};

} // namespace StreamVorti

#endif // STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_

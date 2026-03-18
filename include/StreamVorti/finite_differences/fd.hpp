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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class FiniteDiff
 * \brief Finite difference derivatives base class for structured grids.
 *
 * Constructs sparse derivative matrices using standard central finite
 * difference stencils on structured (Cartesian) MFEM meshes. For nodes at
 * domain boundaries, one-sided (forward/backward) stencils are used for
 * first derivatives, while second derivatives use the same 3-point stencil
 * (which naturally reaches one node inward).
 *
 * Usage:
 *
 *     GridFunction u_gf(fespace_h1);
 *     StreamVorti::FiniteDiff2d fd(u_gf);
 *     fd.Update();
 *     const SparseMatrix &dx = fd.ShapeFunctionDx();
 *     GridFunction dudx(fespace_h1);
 *     dx.Mult(u_gf, dudx);
 */
class FiniteDiff
{
public:
    /*!
     * \brief Constructor for structured MFEM grids.
     *
     * Extracts grid dimensions and spacing from the mesh associated with
     * the given GridFunction. The mesh must be a Cartesian structured grid.
     */
    FiniteDiff(mfem::GridFunction &gf);

    virtual ~FiniteDiff();

    virtual void Update() = 0;

    virtual void SaveDerivToFile(const std::string &deriv, const std::string &filename) const = 0;

    virtual const mfem::SparseMatrix & D(int i) const = 0;

    inline int NumNodes() const { return nnodes_; }

protected:
    mfem::GridFunction *gf_;
    int nnodes_;
    int dim_;
};

} // namespace StreamVorti

#endif // STREAMVORTI_FINITE_DIFFERENCES_FD_HPP_

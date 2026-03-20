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

#ifndef STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_HPP_
#define STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_HPP_

#include "StreamVorti/derivative_operator.hpp"

#ifdef MFEM_USE_MPI

#include "mfem.hpp"

#include <string>

namespace StreamVorti {

/*!
 * \class ParFiniteDiff
 * \brief Parallel finite difference derivatives base class.
 *
 * Returns HypreParMatrix objects for use with parallel MFEM solvers.
 * The structured grid is detected globally, then each rank assembles
 * its local rows using the standard FD stencils.  Off-diagonal entries
 * (referencing DOFs owned by other ranks) are placed in the
 * off-diagonal block of the HypreParMatrix.
 */
class ParFiniteDiff : public ParDerivativeOperator
{
public:
    ParFiniteDiff(mfem::ParGridFunction &gf, int stencil_order = 2);

    virtual ~ParFiniteDiff();

    virtual void Update() = 0;

    virtual void SaveDerivToFile(const std::string &deriv,
                                 const std::string &filename) const = 0;

    virtual const mfem::HypreParMatrix & D(int i) const = 0;

    inline int StencilOrder() const { return stencil_order_; }

protected:
    mfem::ParGridFunction *gf_;
    mfem::ParFiniteElementSpace *pfes_;
    MPI_Comm comm_;
    int rank_;
    int nranks_;
    int stencil_order_;
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_FINITE_DIFFERENCES_PAR_FD_HPP_

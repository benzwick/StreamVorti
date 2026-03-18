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

#ifndef STREAMVORTI_DERIVATIVE_OPERATOR_HPP_
#define STREAMVORTI_DERIVATIVE_OPERATOR_HPP_

#include "mfem.hpp"

#include <string>

namespace StreamVorti {

/*!
 * \class DerivativeOperator
 * \brief Abstract interface for serial spatial derivative operators.
 *
 * Both DCPSE (meshless) and finite differences (structured grid) implement
 * this interface, allowing the solver to use either method through a single
 * pointer without branching on the discretization type.
 *
 * Usage:
 *
 *     DerivativeOperator *op = CreateDerivativeOperator(...);
 *     op->Update();
 *     const SparseMatrix &dx = op->D(0);
 *     dx.Mult(u, dudx);
 */
class DerivativeOperator
{
public:
    virtual ~DerivativeOperator() = default;

    /*! \brief Compute or recompute derivative matrices. */
    virtual void Update() = 0;

    /*! \brief Save a named derivative matrix to file.
     *  \param deriv  Name: "dx", "dy", "dxx", "dyy", "dxy", etc.
     *  \param filename Output path.
     */
    virtual void SaveDerivToFile(const std::string &deriv,
                                 const std::string &filename) const = 0;

    /*! \brief Get first-derivative matrix by index.
     *  \param i  0 → d/dx, 1 → d/dy [, 2 → d/dz]
     */
    virtual const mfem::SparseMatrix & D(int i) const = 0;
};

} // namespace StreamVorti

// ======================================================================
// Parallel version
// ======================================================================

#ifdef MFEM_USE_MPI

namespace StreamVorti {

/*!
 * \class ParDerivativeOperator
 * \brief Abstract interface for parallel spatial derivative operators.
 *
 * Parallel counterpart of DerivativeOperator.  Returns HypreParMatrix
 * instead of SparseMatrix.
 */
class ParDerivativeOperator
{
public:
    virtual ~ParDerivativeOperator() = default;

    virtual void Update() = 0;

    virtual void SaveDerivToFile(const std::string &deriv,
                                 const std::string &filename) const = 0;

    virtual const mfem::HypreParMatrix & D(int i) const = 0;
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_DERIVATIVE_OPERATOR_HPP_

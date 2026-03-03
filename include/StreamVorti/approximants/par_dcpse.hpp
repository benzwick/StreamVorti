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

#ifndef STREAMVORTI_APPROXIMANTS_PAR_DCPSE_HPP_
#define STREAMVORTI_APPROXIMANTS_PAR_DCPSE_HPP_

#include "StreamVorti/support_domain/par_support_domain.hpp"

#ifdef MFEM_USE_MPI

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class ParDcpse
 * \brief Parallel DC-PSE derivatives base class.
 *
 * This is the parallel equivalent of the Dcpse base class. It inherits from
 * ParSupportDomain to access ghost nodes needed for computing DCPSE stencils
 * across partition boundaries.
 *
 * Key differences from serial Dcpse:
 * - Uses ParGridFunction instead of GridFunction
 * - Returns HypreParMatrix instead of SparseMatrix
 * - Requires ghost exchange before computing derivatives
 *
 * Usage:
 *     ParGridFunction u_gf(par_fespace);
 *     StreamVorti::ParDcpse2d derivs(u_gf, 30);
 *     derivs.Update();  // Exchanges ghosts and computes matrices
 *     const HypreParMatrix& dx = derivs.D(0);
 *     ParGridFunction dudx_dcpse(par_fespace);
 *     dx.Mult(u_gf, dudx_dcpse);
 */
class ParDcpse : public ParSupportDomain
{
public:
    /*!
     * \brief ParDcpse constructor to match nodes of an MFEM H1 ParGridFunction.
     * \param gf Parallel grid function containing nodal coordinates
     * \param num_neighbors Number of neighbors for k-NN search
     */
    ParDcpse(mfem::ParGridFunction &gf, int num_neighbors);

    /*!
     * \brief ParDcpse destructor.
     */
    virtual ~ParDcpse();

    /*!
     * \brief Update derivative matrices after mesh or field changes.
     *
     * This method must be called before using the derivative matrices.
     * It triggers:
     * 1. Ghost node exchange (ParSupportDomain::Update())
     * 2. DCPSE stencil computation with extended k-NN
     * 3. HypreParMatrix assembly
     */
    virtual void Update() = 0;

    /*!
     * \brief Save derivative matrix to file for debugging.
     * \param deriv Derivative name ("dx", "dy", "dxx", etc.)
     * \param filename Output file path
     */
    virtual void SaveDerivToFile(const std::string &deriv, const std::string &filename) const = 0;

    /*!
     * \brief Get derivative matrix by index.
     * \param i Derivative index (0=dx, 1=dy for 2D; 0=dx, 1=dy, 2=dz for 3D)
     * \return Reference to HypreParMatrix for the requested derivative
     */
    virtual const mfem::HypreParMatrix & D(int i) const = 0;

    // Limits on condition number of A matrix (same as serial)
    const double cond_A_limit_warn  = 1e6;
    const double cond_A_limit_abort = 1e30;

protected:
    MPI_Comm comm_;  /*!< MPI communicator from ParGridFunction */
    int rank_;       /*!< MPI rank */
    int nranks_;     /*!< Number of MPI ranks */
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_APPROXIMANTS_PAR_DCPSE_HPP_

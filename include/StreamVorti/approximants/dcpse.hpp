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

#ifndef STREAMVORTI_APPROXIMANTS_DCPSE_HPP_
#define STREAMVORTI_APPROXIMANTS_DCPSE_HPP_

#include "mfem.hpp"

#include "StreamVorti/support_domain/support_domain.hpp"

namespace StreamVorti {

/*!
 * \class Dcpse
 * \brief DC PSE derivatives base class.
 */
class Dcpse: public SupportDomain
{
public:
    /*!
     * \brief Dcpse constructor to match nodes of an MFEM H1 GridFunction.
     *
     * Use this in user code:
     *
     *     GridFunction u_gf(fespace_h1);                  // Nodal values
     *     StreamVorti::Dcpse[2d/3d] derivs(u_gf, 30, 5);  // DC PSE derivatives object
     *     derivs.Update();                                // Update/compute matrices
     *     SparseMatrix dx(derivs.ShapeFunctionDx());      // DC PSE derivatives matrix
     *     SparseMatrix dx(derivs.D(0));                   // <--------- (or like this)
     *     GridFunction dudx_dcpse(fespace_h1);            // Derivatives of u wrt x
     *     dx.Mult(u_gf, dudx_dcpse);
     */
    Dcpse(mfem::GridFunction &gf, int NumNeighbors);

    virtual ~Dcpse();

    // Limits on condition number of A matrix
    const double cond_A_limit_warn  = 1e6;
    const double cond_A_limit_abort = 1e6;

    virtual void Update() = 0;

    virtual void SaveDerivToFile(const std::string &deriv, const std::string &filename) const = 0;

    virtual const mfem::SparseMatrix & D(int i) const = 0;

    // TODO: Create these and let user output to file themselves
    // mfem::GridFunction num_support_nodes; // is this already in support domain? Best to put it there instead
    // mfem::GridFunction condition_number;
};

} // namespace StreamVorti

#endif // STREAMVORTI_APPROXIMANTS_DCPSE_HPP_

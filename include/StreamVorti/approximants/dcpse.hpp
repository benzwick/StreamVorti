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

/*!
   \file dcpse.hpp
   \brief Dcpse base class header file.
   \author Benjamin F. Zwick
   \date 06/12/2020
*/

#ifndef STREAMVORTI_APPROXIMANTS_DCPSE_HPP_
#define STREAMVORTI_APPROXIMANTS_DCPSE_HPP_

#include "mfem.hpp"

namespace StreamVorti {

/*!
 *  \addtogroup Approximants
 *  @{
 */

/*!
 * \class Dcpse
 * \brief Base class implemmenting modified DCPSE function derivatives approximants.
 */

class Dcpse
{

public:
    /*!
     * \brief Dcpse constructor to match nodes of an MFEM H1 GridFunction.
     */
    // TODO: move to base class
    // Dcpse(mfem::GridFunction &gf,
    //         int CutoffRadAtNeighbor = 30,
    //         int SupportRadAtNeighbor = 5);

    /*!
     * \brief Dcpse destructor.
     */
    virtual ~Dcpse() = 0;

    inline const mfem::SparseMatrix & D(int i);
};

inline Dcpse::~Dcpse() {};

/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_APPROXIMANTS_DCPSE_HPP_

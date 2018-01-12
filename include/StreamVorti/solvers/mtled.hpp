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



#ifndef STREAMVORTI_SOLVERS_MTLED_HPP_
#define STREAMVORTI_SOLVERS_MTLED_HPP_

/*!
   \file mtled.hpp
   \brief Mtled class header file.
   \author Konstantinos A. Mountris
   \date 27/09/2017
*/


#include "StreamVorti/models/models.hpp"
#include "StreamVorti/approximants/dcpse_2d.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>

#include <stdexcept>
#include <exception>

#include <thread>
#include<mutex>


namespace StreamVorti {

/*!
 *  \addtogroup Solvers
 *  @{
 */


/*!
 * \class Mtled
 * \brief Class implemmenting the Meshfree Total Lagrangian Explicit Dynamics (MTLED) pde solver.
 *
 */

class Mtled {
public:
    /*!
     * \brief Mtled constructor.
     */
    Mtled();


    /*!
     * \brief Mtled destructor.
     */
    virtual ~Mtled();



};


/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_SOLVERS_MTLED_HPP_

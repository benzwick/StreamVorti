/*
 * ExplicitSim - Software for solving PDEs using explicit methods.
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
 *      Grand R. JOLDES
 *      Konstantinos A. MOUNTRIS
 */



#ifndef EXPLICITSIM_CONDITIONS_DIRICHLET_HPP_
#define EXPLICITSIM_CONDITIONS_DIRICHLET_HPP_

/*!
   \file dirichlet.hpp
   \brief Dirichlet class header file.
   \author Konstantinos A. Mountris
   \date 30/09/2017
*/

#include <Eigen/Dense>

#include <string>
#include <vector>

#include <stdexcept>
#include <exception>

namespace ExplicitSim {

/*!
 *  \addtogroup Conditions
 *  @{
 */


/*!
 * \class Dirichlet
 * \brief Class implemmenting the Dirichlet boundary condition.
 *
 */

class Dirichlet {
public:


    /*!
     * \brief The Dirichlet initialization constructor.
     * \param [in] x The conditional stating if the dirichlet condition is applied on the x axis.
     * \param [in] y The conditional stating if the dirichlet condition is applied on the y axis.
     * \param [in] z The conditional stating if the dirichlet condition is applied on the z axis.
     * \param [in] boundary_name The name of the boundary that the dirichlet condition will be assigned.
     */
    Dirichlet(const bool &x, const bool &y, const bool &z, const std::string &boundary_name);


    /*!
     * \brief The Dirichlet destructor.
     */
    virtual ~Dirichlet();


    /*!
     * \brief Edit the container of the nodes indices where the dirichlet condition is applied.
     * \return [std::vector<int>] The container of the nodes indices where the dirichlet condition is applied with editing permission.
     */
    inline std::vector<int> & EditNodesIds() { return this->nodes_ids_; }


    /*!
     * \brief Get the direction at which the dirichlet condition is applied.
     * \return [Eigen::RowVector3i] The direction at which the dirichlet condition is applied.
     */
    inline const Eigen::RowVector3i & Direction() const { return this->direction_; }


    /*!
     * \brief Get the container of the nodes indices where the dirichlet condition is applied.
     * \return [std::vector<int>] The container of the nodes indices where the dirichlet condition is applied.
     */
    inline const std::vector<int> & NodesIds() const { return this->nodes_ids_; }


    /*!
     * \brief Get the index of the boundary that must have the nodes, where the dirichlet condition will be applied.
     * \return [int] The index of the boundary that must have the nodes, where the dirichlet condition will be applied.
     */
    inline const int & BoundaryId() const { return this->boundary_id_; }


    /*!
     * \brief Get the name of the boundary that must have the nodes, where the dirichlet condition will be applied.
     * \return [std::string] The name of the boundary that must have the nodes, where the dirichlet condition will be applied.
     */
    inline const std::string & BoundaryName() const { return this->boundary_name_; }



private:

    /*!
     * \brief The Dirichlet default constructor.
     *
     * Is made private to disable uninitialized declaration of dirichlet conditions.
     *
     */
    Dirichlet();

    Eigen::RowVector3i direction_;         /*!< The direction at which the dirichlet condition is applied. */

    std::vector<int> nodes_ids_;           /*!< The container of the nodes indices where the dirichlet condition is applied. */

    int boundary_id_;                      /*!< The index of the boundary that must have the nodes, where the dirichlet condition will be applied. */

    std::string boundary_name_;            /*!< The name of the boundary that must have the nodes, where the dirichlet condition will be applied. */



};




/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_CONDITIONS_DIRICHLET_HPP_


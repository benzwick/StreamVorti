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



#ifndef EXPLICITSIM_CONDITIONS_LOADING_HPP_
#define EXPLICITSIM_CONDITIONS_LOADING_HPP_

/*!
   \file loading.hpp
   \brief Loading class header file.
   \author Konstantinos A. Mountris
   \date 05/10/2017
*/


#include <ExplicitSim/conditions/load_curve.hpp>

#include <Eigen/Dense>

#include <cmath>
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
 * \class Loading
 * \brief Class implemmenting the loading condition.
 */

class Loading {
public:

    /*!
     * \brief The Loading initialization constructor.
     * \param [in] curve The load curve of the loading condition.
     * \param [in] x The conditional stating if the loading condition is applied on the x axis.
     * \param [in] y The conditional stating if the loading condition is applied on the y axis.
     * \param [in] z The conditional stating if the loading condition is applied on the z axis.
     * \param [in] boundary_name The name of the boundary that the loading condition will be assigned.
     */
    Loading(const LoadCurve &curve, const bool &x, const bool &y, const bool &z, const std::string &boundary_name);


    /*!
     * \brief The Loading default destructor.
     */
    virtual ~Loading();


    /*!
     * \brief Edit the container of the nodes indices where the loading condition is applied.
     * \return [std::vector<int>] The container of the nodes indices where the loading condition is applied with editing permission.
     */
    inline std::vector<int> & EditNodesIds() { return this->nodes_ids_; }


    /*!
     * \brief Get the loading curve of the loading condition.
     * \return [ExplicitSim::LoadCurve] The loading curve of the loading condition.
     */
    inline const LoadCurve & Curve() const { return this->curve_; }


    /*!
     * \brief Get the direction at which the loading condition is applied.
     * \return [Eigen::RowVector3i] The direction at which the loading condition is applied.
     */
    inline const Eigen::RowVector3i & Direction() const { return this->direction_; }


    /*!
     * \brief Get the container of the nodes indices where the loading condition is applied.
     * \return [std::vector<int>] The container of the nodes indices where the loading condition is applied.
     */
    inline const std::vector<int> & NodesIds() const { return this->nodes_ids_; }


    /*!
     * \brief Get the index of the boundary that must have the nodes, where the loading condition will be applied.
     * \return [int] The index of the boundary that must have the nodes, where the loading condition will be applied.
     */
    inline const int & BoundaryId() const { return this->boundary_id_; }


    /*!
     * \brief Get the name of the boundary that must have the nodes, where the loading condition will be applied.
     * \return [std::string] The name of the boundary that must have the nodes, where the loading condition will be applied.
     */
    inline const std::string & BoundaryName() const { return this->boundary_name_; }



private:

    /*!
     * \brief The Loading default constructor.
     *
     * Is made private to disable uninitialized declaration of loading conditions.
     *
     */
    Loading();

    LoadCurve curve_;                      /*!< The loading curve of the loading condition. */

    Eigen::RowVector3i direction_;         /*!< The direction at which the loading condition is applied. */

    std::vector<int> nodes_ids_;           /*!< The container of the nodes indices where the loading condition is applied. */

    int boundary_id_;                      /*!< The index of the boundary that must have the nodes, where the loading condition will be applied. */

    std::string boundary_name_;            /*!< The name of the boundary that must have the nodes, where the loading condition will be applied. */


};




/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_CONDITIONS_LOADING_HPP_

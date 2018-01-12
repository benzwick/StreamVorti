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



#ifndef EXPLICITSIM_CONDITIONS_CONDITIONS_HANDLER_HPP_
#define EXPLICITSIM_CONDITIONS_CONDITIONS_HANDLER_HPP_

/*!
   \file conditions_handler.hpp
   \brief ConditionsHandler class header file.
   \author Konstantinos A. Mountris
   \date 05/10/2017
*/


#include "ExplicitSim/elements/elements.hpp"
#include "ExplicitSim/sets/node_set.hpp"
#include "ExplicitSim/conditions/loading.hpp"
#include "ExplicitSim/conditions/dirichlet.hpp"
#include "ExplicitSim/conditions/ebciem.hpp"
#include "ExplicitSim/conditions/load_curve.hpp"
#include "ExplicitSim/models/weak_model_3d.hpp"
#include "ExplicitSim/support_domain/support_domain.hpp"

#include <Eigen/Dense>

#include <algorithm>
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
 * \class ConditionsHandler
 * \brief Class handling the imposition of the boundary conditions.
 *
 * \note Conditions are handled by reading the boundary id of the mesh.
 * Some nodes have more than one conditions and should see how to process correctly
 * under general conditions for extraction of the nodes ids.
 *
 */

class ConditionsHandler {
public:

    /*!
     * \brief The ConditionsHandler default constructor.
     */
    ConditionsHandler();


    /*!
     * \brief The ConditionsHandler default destructor.
     */
    virtual ~ConditionsHandler();


    /*!
     * \brief Add a loading condition in the imposition handler of boundary conditions.
     * \param [in] load_curve The load curve of the added loading condition for imposition.
     * \param [in] x The conditional stating if the added loading condition is applied on the x axis.
     * \param [in] y The conditional stating if the added loading condition is applied on the y axis.
     * \param [in] z The conditional stating if the added loading condition is applied on the z axis.
     * \param [in] load_boundary_name The name of the boundary that the loading condition will be assigned.
     * \return [void]
     */
    void AddLoading(const LoadCurve &load_curve, const bool &x, const bool &y, const bool &z,
                    const std::string &load_boundary_name);


    /*!
     * \brief Add a dirichlet condition in the imposition handler of boundary conditions.
     * \param [in] x The conditional stating if the added dirichlet condition is applied on the x axis.
     * \param [in] y The conditional stating if the added dirichlet condition is applied on the y axis.
     * \param [in] z The conditional stating if the added dirichlet condition is applied on the z axis.
     * \param [in] dirichlet_boundary_name The name of the boundary that the dirichlet condition will be assigned.
     * \return [void]
     */
    void AddDirichlet(const bool &x, const bool &y, const bool &z, const std::string &dirichlet_boundary_name);


    void AddEbciem(const WeakModel3D &model, const SupportDomain &support,
                   const std::string &mmls_base_func_type, const bool &use_mmls_exact_derivs, const bool &is_simplified);


    /*!
     * \brief Extract the indices of the nodes that boundary conditions will be applied from the given nodesets.
     * \param [in] node_sets The nodesets from where node indices for boundary conditions application will be extracted.
     * \return [void]
     */
    void ExtractBoundaryNodeIds(const std::vector<NodeSet> &node_sets);


    /*!
     * \brief Apply the loading conditions on the displacements matrix.
     * \param [in] time_step The timestep of which the corresponding loading conditions will be applied on the displacements matrix.
     * \param [out] displacements The displacements matrix to be processed for loading conditions application.
     * \return [void]
     * \note Maybe when load is applied in more than one directions should be divided homogeneously. Currently the same load applies to all directions.
     */
    void ApplyLoadingConditions(const int &time_step, Eigen::MatrixXd &displacements) const;


    /*!
     * \brief Reset to zero the forces matrix rows that correspond to nodes where loading conditions are applied.
     * \param [out] forces The forces matrix where the rows corresponding to nodes that loading conditions are applied will be set to zero.
     * \return [void]
     */
    void ResetLoadingConditionsForces(Eigen::MatrixXd &forces) const;


    /*!
     * \brief Apply the dirichlet conditions on the displacements matrix.
     * \param [out] displacements The displacements matrix to be processed for dirichlet conditions application.
     * \return [void]
     */
    void ApplyDirichletConditions(Eigen::MatrixXd &displacements) const;


    /*!
     * \brief Apply the EBCIEM imposition correction on the displacements matrix.
     * \param [out] displacements The displacements matrix to be processed for EBCIEM imposition correction.
     * \return [void]
     * \note Maybe when load is applied in more than one directions should be divided homogeneously. Currently the same load applies to all directions.
     */
    void ApplyEbciem(const int &load_time_step, Eigen::MatrixXd &displacements) const;


    /*!
     * \brief Get the loading conditions of the conditions handler.
     * \return [std::vector<ExplicitSim::Loading>] The loading conditions of the conditions handler.
     */
    inline const std::vector<Loading> & LoadingConds() const { return this->loading_conds_; }


    /*!
     * \brief Get the dirichlet conditions of the conditions handler.
     * \return [std::vector<ExplicitSim::Dirichlet>] The dirichlet conditions of the conditions handler.
     */
    inline const std::vector<Dirichlet> & DirichletConds() const { return this->dirichlet_conds_; }


    /*!
     * \brief Get the boolean declaring the status of boundary nodes indices extraction.
     * \return [bool] The boolean declaring the status of boundary nodes indices extraction.
     */
    inline const bool & BoundaryNodeIdsAreExtracted() const { return this->bound_node_ids_are_extracted_; }



private:
    std::vector<Loading> loading_conds_;         /*!< The container of the loading conditions to be imposed. */

    std::vector<Dirichlet> dirichlet_conds_;     /*!< The container of the dirichlet conditions to be imposed. */

    bool bound_node_ids_are_extracted_;          /*!< Boolean declaring if the nodes indices where boundary conditions are applying have been extracted. */

    Ebciem ebciem_;                              /*!< The EBCIEM correction for essential conditions imposition. */

};




/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_CONDITIONS_CONDITIONS_HANDLER_HPP_

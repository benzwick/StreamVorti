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


#include "ExplicitSim/conditions/conditions_handler.hpp"


namespace ExplicitSim {

ConditionsHandler::ConditionsHandler() : bound_node_ids_are_extracted_(false)
{}


ConditionsHandler::~ConditionsHandler()
{}


void ConditionsHandler::AddLoading(const LoadCurve &load_curve, const bool &x, const bool &y, const bool &z, const std::string &load_boundary_name)
{
    // Convert the boundary name to low case.
    std::string low_case_boundary_name = load_boundary_name;
    std::transform(low_case_boundary_name.begin(), low_case_boundary_name.end(), low_case_boundary_name.begin(), ::tolower);

    // Initialize a new loading condition with the input variables.
    Loading new_load(load_curve, x, y, z, low_case_boundary_name);

    // Add the new loading condition in the loadings container.
    this->loading_conds_.push_back(new_load);

}


void ConditionsHandler::AddDirichlet(const bool &x, const bool &y, const bool &z, const std::string &dirichlet_boundary_name)
{
    // Convert the boundary name to low case.
    std::string low_case_boundary_name = dirichlet_boundary_name;
    std::transform(low_case_boundary_name.begin(), low_case_boundary_name.end(), low_case_boundary_name.begin(), ::tolower);

    // Initialize a new dirichlet condition with the input variables.
    Dirichlet new_dirichlet(x, y, z, low_case_boundary_name);

    // Add the new dirichlet condition in the dirichlet conditions container.
    this->dirichlet_conds_.push_back(new_dirichlet);

}


void ConditionsHandler::AddEbciem(const WeakModel3D &model, const SupportDomain &support,
                                  const std::string &mmls_base_func_type, const bool &use_mmls_exact_derivs, const bool &is_simplified)
{
    // Check if dirichlet and loading boundary nodes indices have been set before adding the EBCIEM correction.
    if (!this->bound_node_ids_are_extracted_) {
        throw std::invalid_argument(Logger::Error("Could not apply EBCIEM correction. "
                                                  "Node indices for dirichlet or loading conditions are not extracted.").c_str());
    }

    // Set the properties of the model's mmls approximants for the EBCIEM.
    this->ebciem_.SetNodalMmlsProperties(mmls_base_func_type, use_mmls_exact_derivs);

    // Compute the EBCIEM correction matrix.
    this->ebciem_.ComputeCorrectionMatrices(model, support, this->dirichlet_conds_, this->loading_conds_, is_simplified);
}


void ConditionsHandler::ExtractBoundaryNodeIds(const std::vector<NodeSet> &node_sets)
{
    // Check if any conditions have been added in the conditions handler.
    if ((this->loading_conds_.size() == 0) && (this->dirichlet_conds_.size() == 0) ) {
        std::runtime_error(Logger::Error("Cannot extract node indices for conditions application. "
                                         "No conditions available in the conditions handler.").c_str());
    }

    // Clear the loading conditions nodes indices.
    if (this->loading_conds_.size() != 0) {
        for (auto &loading : this->loading_conds_) {
            loading.EditNodesIds().clear();
        }
    }

    // Clear the dirichlet conditions nodes indices.
    if (this->dirichlet_conds_.size() != 0) {
        for (auto &dirichlet : this->dirichlet_conds_) {
            dirichlet.EditNodesIds().clear();
        }
    }

    // Iterate over node sets to extract indices of the nodes where conditions will be applied.
    for (auto &nset : node_sets) {

        // Extract nodes indices for the loading conditions.
        if (this->loading_conds_.size() != 0) {
            // Iterate over the loading conditions.
            for (auto &loading : this->loading_conds_) {
                // Check if nset's name complies with the loading condition's boundary name.
                if (loading.BoundaryName() ==  nset.Name()) {
                    // Assign the nset's indices to the loading condition's nodes indices.
                    loading.EditNodesIds() = nset.NodeIds();
                }
            }
        } // End of Extract nodes indices for the loading conditions.

        // Extract nodes indices for the dirichlet conditions.
        if (this->dirichlet_conds_.size() != 0) {
            // Iterate over the dirichlet conditions.
            for (auto &dirichlet : this->dirichlet_conds_) {
                // Check if nset's name complies with the dirichlet condition's boundary name.
                if (dirichlet.BoundaryName() == nset.Name()) {
                    // Assign the nset's indices to the dirichlet condition's nodes indices.
                    dirichlet.EditNodesIds() = nset.NodeIds();
                }

            }
        } // End of Extract nodes indices for the dirichlet conditions.

    } // End of Iteration over nodes.


    // Check if condition nodes indices were extracted for all the loading conditions.
    for (auto &loading : this->loading_conds_) {
        // The loading condition's index.
        auto loading_id = &loading - &this->loading_conds_[0];

        if (loading.NodesIds().size() == 0) {
            std::cout << Logger::Warning("No condition nodes indices were found for the loading condition [")
                      << loading_id << "]. Remove the condition if is not used, otherwise check the given mesh file.\n";
        }
    }

    // Check if condition nodes indices were extracted for all the dirichlet conditions.
    for (auto &dirichlet : this->dirichlet_conds_) {
        // The dirichlet condition's index.
        auto dirichlet_id = &dirichlet - &this->dirichlet_conds_[0];

        if (dirichlet.NodesIds().size() == 0) {
            std::cout << Logger::Warning("No condition nodes indices were found for the dirichlet condition [")
                      << dirichlet_id << "]. Remove the condition if is not used, otherwise check the given mesh file.\n";
        }
    }

    // Set the boolean declaring the nodes indices extraction to true.
    this->bound_node_ids_are_extracted_ = true;


}


void ConditionsHandler::ApplyLoadingConditions(const int &time_step, Eigen::MatrixXd &displacements) const
{
    // Check if loading conditions have been added in the conditions handler.
    if (this->loading_conds_.size() == 0) {
        std::runtime_error(Logger::Error("Cannot apply loading conditions. No loading "
                                         "conditions available in the conditions handler.").c_str());
    }

    // Iterate over loading conditions.
    for (auto &loading : this->loading_conds_) {
        // Iterate over loading condition's nodes indices for loading application in the displacement matrix.
        for (auto &node_id : loading.NodesIds()) {
            // Application on the X direction.
            if (loading.Direction().coeff(0) != 0) { displacements.coeffRef(node_id, 0) = loading.Curve().LoadDispAt(time_step); }

            // Application on the Y direction.
            if (loading.Direction().coeff(1) != 0) { displacements.coeffRef(node_id, 1) = loading.Curve().LoadDispAt(time_step); }

            // Application on the Z direction.
            if (loading.Direction().coeff(2) != 0) { displacements.coeffRef(node_id, 2) = loading.Curve().LoadDispAt(time_step); }
        }

    }

}


void ConditionsHandler::ResetLoadingConditionsForces(Eigen::MatrixXd &forces) const
{
    // Check if loading conditions have been added in the conditions handler.
    if (this->loading_conds_.size() == 0) {
        std::runtime_error(Logger::Error("Cannot apply loading conditions. No loading "
                                         "conditions available in the conditions handler.").c_str());
    }

    // Iterate over loading conditions.
    for (auto &loading : this->loading_conds_) {
        // Iterate over loading condition's nodes indices to reset force to zero.
        for (auto &node_id : loading.NodesIds()) {
            // Application on the X direction.
            if (loading.Direction().coeff(0) != 0) { forces.coeffRef(node_id, 0) = 0.; }

            // Application on the Y direction.
            if (loading.Direction().coeff(1) != 0) { forces.coeffRef(node_id, 1) = 0.; }

            // Application on the Z direction.
            if (loading.Direction().coeff(2) != 0) { forces.coeffRef(node_id, 2) = 0.; }
        }
    } // End Iterate over loading conditions.

}


void ConditionsHandler::ApplyDirichletConditions(Eigen::MatrixXd &displacements) const
{
    // Check if dirichlet conditions have been added in the conditions handler.
    if (this->dirichlet_conds_.size() == 0) {
        std::runtime_error(Logger::Error("Cannot apply dirichlet conditions. No dirichlet "
                                         "conditions available in the conditions handler.").c_str());
    }

    // Iterate over dirichlet conditions.
    for (auto &dirichlet : this->dirichlet_conds_) {
        // Iterate over dirichlet condition's nodes indices for dirichlet application.
        for (auto &node_id : dirichlet.NodesIds()) {
            // Application on the X direction.
            if (dirichlet.Direction().coeff(0) != 0) { displacements.coeffRef(node_id, 0) = 0.; }

            // Application on the Y direction.
            if (dirichlet.Direction().coeff(1) != 0) { displacements.coeffRef(node_id, 1) = 0.; }

            // Application on the Z direction.
            if (dirichlet.Direction().coeff(2) != 0) { displacements.coeffRef(node_id, 2) = 0.; }
        }
    }

}


void ConditionsHandler::ApplyEbciem(const int &load_time_step, Eigen::MatrixXd &displacements) const
{
    // Check if dirichlet and loading conditions have been added in the conditions handler.
    if (this->dirichlet_conds_.size() == 0 || this->loading_conds_.size() == 0) {
        std::runtime_error(Logger::Error("Cannot apply dirichlet conditions. No dirichlet of loading "
                                         "conditions available in the conditions handler.").c_str());
    }

    // Create modified displacement matrix.
    Eigen::MatrixXd mod_disp = this->ebciem_.NodalMmls().ShapeFunction().transpose()*displacements;

    // Initialize total imposition matrix for displacements correction.
    Eigen::MatrixXd total_imposed = Eigen::MatrixXd::Zero(displacements.rows(), displacements.cols());

    // Iterate over dirichlet conditions.
    for (auto &dirichlet : this->dirichlet_conds_) {
        // The condition's index in the container.
        auto cond_id = &dirichlet - &this->dirichlet_conds_[0];

        // Initialize dirichlet imposition matrix for displacements correction.
        Eigen::MatrixXd dirichlet_imposed = Eigen::MatrixXd::Zero(dirichlet.NodesIds().size(), displacements.cols());

        // Iterate over dirichlet condition's nodes indices for dirichlet application.
        for (auto &node_id : dirichlet.NodesIds()) {
            // The row index of the imposed_disp matrix.
            auto row_id = &node_id - &dirichlet.NodesIds()[0];

            // Application on the X direction.
            if (dirichlet.Direction().coeff(0) != 0) { dirichlet_imposed.coeffRef(row_id, 0) = -1.*mod_disp.coeff(node_id, 0); }

            // Application on the Y direction.
            if (dirichlet.Direction().coeff(1) != 0) { dirichlet_imposed.coeffRef(row_id, 1) = -1.*mod_disp.coeff(node_id, 1); }

            // Application on the Z direction.
            if (dirichlet.Direction().coeff(2) != 0) { dirichlet_imposed.coeffRef(row_id, 2) = -1.*mod_disp.coeff(node_id, 2); }
        }

        // Update total imposition matrix (add dirichlet correction matrix * dirichlet imposition matrix).
        total_imposed += this->ebciem_.DirichletCorrMats()[cond_id]*dirichlet_imposed;

    } // End Iterate over dirichlet conditions.

    // Iterate over loading conditions.
    for (auto &loading : this->loading_conds_) {
        // The condition's index in the container.
        auto cond_id = &loading - &this->loading_conds_[0];

        // Initialize loading imposition matrix for displacements correction.
        Eigen::MatrixXd loading_imposed = Eigen::MatrixXd::Zero(loading.NodesIds().size(), displacements.cols());

        // Iterate over loading condition's nodes indices for loading application.
        for (auto &node_id : loading.NodesIds()) {
            // The row index of the imposed_disp matrix.
            auto row_id = &node_id - &loading.NodesIds()[0];

            // Application on the X direction.
            if (loading.Direction().coeff(0) != 0) {
                loading_imposed.coeffRef(row_id, 0) = loading.Curve().LoadDispAt(load_time_step) - mod_disp.coeff(node_id, 0);
            }

            // Application on the Y direction.
            if (loading.Direction().coeff(1) != 0) {
                loading_imposed.coeffRef(row_id, 1) = loading.Curve().LoadDispAt(load_time_step) - mod_disp.coeff(node_id, 1);
            }

            // Application on the Z direction.
            if (loading.Direction().coeff(2) != 0) {
                loading_imposed.coeffRef(row_id, 2) = loading.Curve().LoadDispAt(load_time_step) - mod_disp.coeff(node_id, 2);
            }

        } // End Iterate over loading condition's nodes indices for loading application.

        // Update total imposition matrix (add loading correction matrix * loading imposition matrix).
        total_imposed += this->ebciem_.LoadingCorrMats()[cond_id]*loading_imposed;

    } // End Iterate over loading conditions.

    // Apply imposition correction to displacements.
    displacements += total_imposed;

}


} //end of namespace ExplicitSim


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


#include "ExplicitSim/solvers/mtled.hpp"
#include<ctime>
namespace ExplicitSim {

Mtled::Mtled() : min_step_(0.), max_step_(0.), stable_step_(0.), total_time_steps_num(0), save_progress_steps_(1)
{}


Mtled::~Mtled()
{}


void Mtled::ComputeTimeSteps(const std::vector<double> &wave_speed,
                             const std::vector< std::vector<int> > &neighbors_ids,
                             const Mmls3d &model_approximant)
{
    // Check that size of containers is consistent.
    if ( (wave_speed.size() != neighbors_ids.size()) ||
         (static_cast<int>(wave_speed.size()) != model_approximant.ShapeFunctionDx().cols()) ||
         (static_cast<int>(wave_speed.size()) != model_approximant.ShapeFunctionDy().cols()) ||
         (static_cast<int>(wave_speed.size()) != model_approximant.ShapeFunctionDz().cols()) ) {

        throw std::invalid_argument(Logger::Error("Could not compute time steps. "
                                                  "Check the size consistency of the input variables.").c_str());
    }

    // Clear the solver's time steps container.
    this->time_steps_.clear();

    // Compute time step for each evaluation point.
    double step = 0.;
    for (auto &point_speed : wave_speed) {

        // The index of the ith evaluation point.
        auto i = &point_speed - &wave_speed[0];

        // Gather x, y, z derivatives in single matrix.
        Eigen::MatrixXd derivs_mat(neighbors_ids[i].size(), 3);

        int row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDx(),i); it; ++it) {
            derivs_mat(row_id, 0) = it.value();
            row_id++;
        }

        row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDy(),i); it; ++it) {
            derivs_mat(row_id, 1) = it.value();
            row_id++;
        }

        row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDz(),i); it; ++it) {
            derivs_mat(row_id, 2) = it.value();
            row_id++;
        }

        // Square the elements of the derivatives matrix.
        derivs_mat = derivs_mat.cwiseProduct(derivs_mat);

        // Compute the time step for the ith evaluation point.
        step = point_speed * std::sqrt(neighbors_ids[i].size() * derivs_mat.sum());
        step = 2. / step;

        // Store the time step for the ith evaluation point.
        this->time_steps_.emplace_back(step);

    }

    // Get the value of the minimum time step.
    this->min_step_ = *std::min_element(this->time_steps_.begin(), this->time_steps_.end());

    // Get the value of the maximum time step.
    this->max_step_ = *std::max_element(this->time_steps_.begin(), this->time_steps_.end());

}


void Mtled::ComputeStableStep(const std::vector<double> &mass, const bool &is_mass_scaled, double safety_factor)
{
    // Check if mass has been computed.
    if (mass.size() == 0) {
        std::string error = "[ExplicitSim ERROR] Cannot compute stable step. The given mass container is empty.";
        throw std::invalid_argument(error.c_str());
    }

    // Check if time steps have been already computed.
    if (this->time_steps_.size() == 0) {
        std::string error = "[ExplicitSim ERROR] Cannot compute stable step. Time steps have not been computed.";
        throw std::invalid_argument(error.c_str());
    }

    // Set the stable time step according to the mass scaling.
    if (is_mass_scaled) { this->stable_step_ = this->max_step_ / safety_factor; }
    else { this->stable_step_ = this->min_step_ / safety_factor; }
}


void Mtled::ComputeTotalTimeStepsNum(const int &load_steps_num, const int &equilibrium_steps_num)
{
    // Check if the number of loading application time steps has been initialized.
    if (load_steps_num <= 0) {
        std::string error = "[ExplicitSim ERROR] Cannot compute the total time steps number for the explicit solution. "
                            "The number of loading application time steps is not initialized.";
        throw std::invalid_argument(error.c_str());
    }

    // Check if the number of equilibrium application time steps has been initialized.
    if (equilibrium_steps_num <= 0) {
        std::string error = "[ExplicitSim ERROR] Cannot compute the total time steps number for the explicit solution. "
                            "The number of dynamic relaxation equilibrium time steps is not initialized.";
        throw std::invalid_argument(error.c_str());
    }

    this->total_time_steps_num = load_steps_num + equilibrium_steps_num;

}


void Mtled::Solve(const WeakModel3D &weak_model_3d, const std::vector<std::vector<int> > &neighbor_ids, const ConditionsHandler &cond_handler,
                  const Mmls3d &model_approximant, const NeoHookean &material, const DynRelaxProp &dyn_relax_prop, const bool &use_ebciem)
{

    // Check if the dynamic relaxation properties have been initialized.
    if (!dyn_relax_prop.IsInitialized()) {
        std::string error = "[ExplicitSim ERROR] Cannot generate the explicit dynamics solution. "
                            "One or more dynamic relaxation properties have not been initialized.";
        throw std::invalid_argument(error.c_str());
    }

    // Displacements and forces matrices initialization.
    Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);
    Eigen::MatrixXd disp_new = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);
    Eigen::MatrixXd disp_old = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);
    Eigen::MatrixXd disp_saved = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);
    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);
    Eigen::MatrixXd forces_saved = Eigen::MatrixXd::Zero(weak_model_3d.Grid().NodesNum(), 3);

    // Clear saved disps and forces.
    this->saved_disps_.clear();
    this->saved_forces_.clear();

    // Reserve memory for the states to be recorded
    if (this->save_progress_steps_ != 0) {
        this->saved_disps_.reserve(this->total_time_steps_num/this->save_progress_steps_);
        this->saved_forces_.reserve(this->total_time_steps_num/this->save_progress_steps_);
    }
    else {  // only for initial and final state.
        this->saved_disps_.reserve(2);
        this->saved_forces_.reserve(2);
    }

    // Store initial disps and forces.
    this->saved_disps_.emplace_back(disp);
    this->saved_forces_.emplace_back(forces);

    // Retrieve the load steps number from the total steps and the equilibrium steps difference.
    int step_num_load = this->total_time_steps_num - dyn_relax_prop.EquilibriumStepsNum();

    // Dynamic Relaxation variables.
    bool stabilized_conv_rate = false;
    bool conv_disp_updated = false;

    double conv_rate = dyn_relax_prop.LoadConvRate();
    double old_conv_rate = dyn_relax_prop.LoadConvRate();

    int termination_count = 0;
    int samples_num = 0;
    int no_update_steps_num = 0;

    // Apply load condition at first time step (0) on new displacements.
    //cond_handler.ApplyLoadingConditions(0, disp_new);

    // Iterate over integration points to collect xyz derivatives matrices.
    std::vector<Eigen::MatrixXd> deriv_mats;
    deriv_mats.reserve(weak_model_3d.IntegrationPoints().PointsNum());
    for (int ip = 0; ip != weak_model_3d.IntegrationPoints().PointsNum(); ++ip) {

        // Gather x, y, z derivatives in single matrix.
        Eigen::MatrixXd xyz_derivs(neighbor_ids[ip].size(), 3);

        int row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDx(),ip); it; ++it) {
            xyz_derivs(row_id, 0) = it.value();
            row_id++;
        }

        row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDy(),ip); it; ++it) {
            xyz_derivs(row_id, 1) = it.value();
            row_id++;
        }

        row_id = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(model_approximant.ShapeFunctionDz(),ip); it; ++it) {
            xyz_derivs(row_id, 2) = it.value();
            row_id++;
        }

        // Store the derivatives matrix for each integration point in the container.
        deriv_mats.emplace_back(xyz_derivs);

    }

    // Iterate over the total number of time steps.
    int steps_counter = 0;
    for (auto step = 0; step != this->total_time_steps_num; ++step) {
        // Increase steps_counter to count the performing steps.
        steps_counter++;

        // Update displacements and make force zero (note the order).
        disp_old = disp;
        disp = disp_new;
        forces.setZero(weak_model_3d.Grid().NodesNum(), 3);

        // Compute the forces at each time step.
        this->ComputeForces(weak_model_3d, neighbor_ids, deriv_mats, material, disp, forces);

        // Update saved displacements, forces and convergence rates.
        if (steps_counter == step_num_load) {
            disp_saved = disp;
            forces_saved = forces;
            conv_rate = dyn_relax_prop.AfterLoadConvRate();
            old_conv_rate = dyn_relax_prop.AfterLoadConvRate();
        }

        // Use forces to update displacements using explicit integration and
        // mass proportional damping (Dynamic Relaxation)
        double f8x = (conv_rate + 1.) * (this->stable_step_/2.);

        // Compute new displacements.
        disp_new = -f8x*f8x*(forces.array() / weak_model_3d.MassMatrix().array()).matrix() -
                    conv_rate*conv_rate*disp_old + (1. + conv_rate*conv_rate)*disp;

        // Apply boundary conditions.
        if (use_ebciem) { // EBCIEM imposition of boundary conditions.

            // Apply EBCIEM for suitable loading step.
            if (steps_counter < step_num_load) { cond_handler.ApplyEbciem(steps_counter, disp_new); }
            else { cond_handler.ApplyEbciem(step_num_load-1, disp_new); }

        }
        else { // Direct imposition of boundary conditions.

            // Apply load conditions on the new displacements for suitable loading step.
            if (steps_counter < step_num_load) { cond_handler.ApplyLoadingConditions(steps_counter, disp_new); }
            else { cond_handler.ApplyLoadingConditions(step_num_load-1, disp_new); }

            // Apply dirichlet conditions on the new displacements.
            cond_handler.ApplyDirichletConditions(disp_new);

        } // End Apply boundary conditions.

        // Check for divergence.
        double max_current_disp = disp_new.cwiseAbs().col(0).maxCoeff();
        double max_load_disp = 1.5 * cond_handler.LoadingConds()[0].Curve().MaxDisplacement();


        // Use squared values to avoid using absolute.
        if (max_current_disp*max_current_disp > max_load_disp*max_load_disp) {
            std::cout << "[ExplicitSim WARNING] MTLED solution has become unbounded at step " +
                                std::to_string(step) + ". Reduce time step.";
            // Stop solution if become unstable and return.
            return;
        }

        // Termination criteria and convergence rate.
        if (steps_counter > step_num_load) {
            // Estimate the lower oscilation frequency.
            Eigen::MatrixXd disp_diff = std::move(disp - disp_saved);
            Eigen::MatrixXd force_diff = std::move(forces - forces_saved);

            // Apply loading conditions to force difference matrix.
            cond_handler.ResetLoadingConditionsForces(force_diff);
            cond_handler.ApplyDirichletConditions(force_diff);

            double k_sum = (disp_diff.array() * force_diff.array()).sum();
            double m_sum = (disp_diff.array() * disp_diff.array() * weak_model_3d.MassMatrix().array()).sum();

            // Update convergence rate adaptively.
            if (steps_counter < (step_num_load + dyn_relax_prop.StopUpdateConvRateStepsNum()) ) {
                // Ensure k_sum is possitive.
                k_sum = std::abs(k_sum);

                // Updates in convergence rates.
                if ((m_sum > 1.e-13) && (k_sum > 1.e-13)) {
                    // Reset no update steps number.
                    no_update_steps_num = 0;

                    // Minimum frequency.
                    double min_freq = std::sqrt(k_sum / m_sum);

                    // Kmat condition number (square root).
                    double k_cond_number_root = 2. / (min_freq * this->stable_step_);

                    double temp_conv_rate = (k_cond_number_root - 1.) / (k_cond_number_root + 1.);

                    if (std::abs(temp_conv_rate - old_conv_rate) < dyn_relax_prop.ConvRateDeviation()) {
                        samples_num++;
                        if (samples_num >= dyn_relax_prop.StableConvRateStepsNum()) {
                            if (conv_disp_updated) {
                                conv_disp_updated = false;

                                // Update stabilized state if convergence rate deviation criterion is satisfied.
                                if (std::abs(temp_conv_rate - conv_rate) < dyn_relax_prop.ConvRateStopDeviation()) { stabilized_conv_rate = true; }

                                // Update convergence rate.
                                conv_rate = temp_conv_rate;
                            }
                        }
                    }
                    else {
                        // Reset number of samples to zero.
                        samples_num = 0;
                    }

                    // Update old convergence rate.
                    old_conv_rate = temp_conv_rate;

                } // End of Updates in convergence rates.

            }
            else {
                no_update_steps_num++;
                if (no_update_steps_num >= 10*dyn_relax_prop.StableConvRateStepsNum()) { stabilized_conv_rate = true; }
            }


            // Check termination criteria.
            if (stabilized_conv_rate) {

                // Compute maximum displacement variation.
                double max_disp_var = (disp_new - disp).cwiseAbs().maxCoeff();

                // Adjust convergence rate value.
                double estim_conv_rate = conv_rate + dyn_relax_prop.StopConvRateError()*(1. - conv_rate);
                double estim_error = max_disp_var * estim_conv_rate / (1. - estim_conv_rate);

                if (estim_error < dyn_relax_prop.StopAbsError()) {
                    termination_count++;
                    if (termination_count >= dyn_relax_prop.StopStepsNum()) {
                        std::cout << "[ExplicitSim] MTLED solution tolerance has been satisfied at step: " << step+1 << "\n";
                        break;
                    }
                }
                else { termination_count = 0; }

            } //End of check termination criteria.

            // Update saved forces and displacements.
            if ( (steps_counter - step_num_load) % dyn_relax_prop.ForceDispUpdateStepsNum() == 0 ) {
                disp_saved = disp;
                forces_saved = forces;
                conv_disp_updated = true;
            }

        } // End of Termination criteria and convergence rate.

        // Output MTLED progress.
        if (this->save_progress_steps_ != 0) {
            if(steps_counter % this->save_progress_steps_ == 0) {
                this->saved_disps_.emplace_back(disp);
                this->saved_forces_.emplace_back(forces);
                std::cout << "[ExplicitSim] MTLED solver completed " << steps_counter
                          << " / " << this->total_time_steps_num << " steps\n";
            }
        }

    } //End of time steps iteration.

    // Check for convergence satisfaction.
    if (steps_counter == this->total_time_steps_num) {
        std::string error = "[ExplicitSim ERROR] MTLED solution convergence rate at the end of the simulation: "
                + std::to_string(conv_rate) + ".\n Solution tolerance has not been satisfied. Reduce the size of the stable time step.";
        throw std::runtime_error(error.c_str());
    }
    else {
        std::cout << Logger::Message("Convergence rate of MTLED solution at termination: ") << conv_rate << std::endl;
    }


    // Store final disps and forces if they haven't been stored during progress storing.
    if (this->save_progress_steps_ != 0) {
        if (steps_counter % this->save_progress_steps_ != 0) {
            this->saved_disps_.emplace_back(disp);
            this->saved_forces_.emplace_back(forces);
        }
    }
    else {  // Store final state
        this->saved_disps_.emplace_back(disp);
        this->saved_forces_.emplace_back(forces);
    }

    // Check size consistency between displacements and forces containers.
    if (this->saved_disps_.size() != this->saved_forces_.size()) {
        std::string error = "[ExplicitSim ERROR] The saved displacements and saved forces do not correspond to the same number of steps.";
        throw std::runtime_error(error.c_str());
    }


}


void Mtled::ApplyShapeFuncToDisplacements(const WeakModel3D &weak_model_3d, const Mmls3d &nodal_approximant)
{
    // Create matrix with initial nodal positions.
    Eigen::MatrixXd nodal_pos_init(weak_model_3d.TetrahedralMesh().NodesNum(), 3);

    // Assemble initial nodal positions.
    for (const auto &node : weak_model_3d.TetrahedralMesh().Nodes()) {
        auto row_id = &node - &weak_model_3d.TetrahedralMesh().Nodes()[0];

        nodal_pos_init.coeffRef(row_id, 0) = node.Coordinates().X();
        nodal_pos_init.coeffRef(row_id, 1) = node.Coordinates().Y();
        nodal_pos_init.coeffRef(row_id, 2) = node.Coordinates().Z();

    }

    // Iterate over saved displacements. Skip the first one (zero displacements)
    for (std::size_t disp_id = 1; disp_id != this->saved_disps_.size(); ++disp_id) {

        // Compute the displaced nodal positions.
        Eigen::MatrixXd nodal_pos_disp = nodal_pos_init + this->saved_disps_[disp_id];

        // Iterate over model nodes.
        for (const auto &node : weak_model_3d.TetrahedralMesh().Nodes()) {
            // The node index.
            auto node_id = &node - &weak_model_3d.TetrahedralMesh().Nodes()[0];

            // Initialize the final nodal position.
            Vec3<double> final_pos(0., 0., 0.);

            // Add the values of the displaced nodal positions * the shape function value for the neighbor nodes of the nodal node.
            for (Eigen::SparseMatrix<double>::InnerIterator sf(nodal_approximant.ShapeFunction(),node_id); sf; ++sf) {
                final_pos.SetX(final_pos.X() + nodal_pos_disp.coeff(sf.row(), 0)*sf.value());
                final_pos.SetY(final_pos.Y() + nodal_pos_disp.coeff(sf.row(), 1)*sf.value());
                final_pos.SetZ(final_pos.Z() + nodal_pos_disp.coeff(sf.row(), 2)*sf.value());
            }

            // Update the nodal values of the saved displacements with the final nodal positions.
            this->saved_disps_[disp_id].coeffRef(node_id, 0) = final_pos.X();
            this->saved_disps_[disp_id].coeffRef(node_id, 1) = final_pos.Y();
            this->saved_disps_[disp_id].coeffRef(node_id, 2) = final_pos.Z();

        }

        // Isolate only the displacement values.
        this->saved_disps_[disp_id] -= nodal_pos_init;

    } // End Iterate over saved displacements.


}


void Mtled::ComputeForces(const WeakModel3D &weak_model_3d, const std::vector<std::vector<int> > &neighbor_ids,
                          const std::vector<Eigen::MatrixXd> &deriv_mats, const NeoHookean &material,
                          const Eigen::MatrixXd &displacements, Eigen::MatrixXd &forces)
{
    // Iterate over integration points for force generation.
    for (auto &int_point_weight : weak_model_3d.IntegrationPoints().Weights()) {

        // Initialize deformation gradient and 2nd Piola-Kirchhoff stress tensor.
        Eigen::Matrix3d FT = Eigen::Matrix3d::Zero(3, 3);
        Eigen::Matrix3d spk_stress = Eigen::Matrix3d::Zero(3, 3);

        // The integration point's index.
        auto int_point_id = &int_point_weight - &weak_model_3d.IntegrationPoints().Weights()[0];

        // Local displacements and forces at integration point's support domain.
        Eigen::MatrixXd disp_local = Eigen::MatrixXd::Zero(neighbor_ids[int_point_id].size(), 3);
        Eigen::MatrixXd forces_local = Eigen::MatrixXd::Zero(neighbor_ids[int_point_id].size(), 3);

        // Iterate over neighbor nodes indices.
        for (auto &neigh_id : neighbor_ids[int_point_id]) {
            // The neighbors index in the container.
            auto id = &neigh_id - &neighbor_ids[int_point_id][0];

            // Populate disp_local.
            disp_local.row(id) = displacements.row(neigh_id);

        }

        // Compute deformation gradient.
        FT.noalias() = deriv_mats[int_point_id].transpose() * disp_local;
        // Add identity matrix contribution to diagonal elements of the deformation gradient tensor.
        FT.coeffRef(0,0) += 1.; FT.coeffRef(1,1) += 1.; FT.coeffRef(2,2) += 1.;

        // Compute the 2nd Piola-Kirchhoff stress tensor.
        spk_stress.noalias() = std::move(material.SpkStress(FT, int_point_id));

        // Compute the force contribution of the current integration point.
        forces_local = deriv_mats[int_point_id] * spk_stress.transpose() * FT * int_point_weight;

        // Update total force.
        for (auto &neigh_id : neighbor_ids[int_point_id]) {
            // The neighbors index in the container.
            auto id = &neigh_id - &neighbor_ids[int_point_id][0];

            // Add integration point's contribution in force matrix.
            forces.row(neigh_id) += forces_local.row(id);
        }

    } // End iteration over integration points.
}


} //end of namespace ExplicitSim

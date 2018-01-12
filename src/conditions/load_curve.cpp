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


#include "ExplicitSim/conditions/load_curve.hpp"


namespace ExplicitSim {

LoadCurve::LoadCurve() : load_time_(0.), max_displacement_(0.), load_steps_num_(0)
{}


LoadCurve::~LoadCurve()
{}


void LoadCurve::ComputeLoadStepsNum(const double &solver_time_step)
{
    // Check if loading curve application time is initialized to a "reasonable" value.
    if (this->load_time_ < 0.0000001) {
        std::string error = "[ExplicitSim ERROR] Cannot compute number of load curve application steps. "
                            "The load curve application time has not been set or is too small.";
        throw std::runtime_error(error.c_str());
    }

    // Compute the number of time steps for load curve application.
    this->load_steps_num_ = static_cast<int>(std::floor((this->load_time_ / solver_time_step) + 0.5));

}


void LoadCurve::ComputeLoadStepDisplacements(const double &solver_time_step)
{

    // Check if the maximum displacement that is applied by the load curve is initialized.
    if (std::abs(this->max_displacement_) < 0.0000001) {
        std::string error = "[ExplicitSim ERROR] Cannnot compute the displacement variation in each loading time step. "
                            "The maximum displacement to be applied by the load curve is not initialized.";
        throw std::runtime_error(error.c_str());
    }

    // Check if the number of time steps for load curve application has been computed.
    if (this->load_steps_num_ <= 0) {
        std::string error = "[ExplicitSim ERROR] Cannnot compute the displacement variation in each loading time step. "
                            "The number of loading steps has not been computed.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the container of the displacement variation in each loading time step.
    this->load_step_displacements_.clear();

    // The total loading application time with respect to the solver's time step.
    double total_load_time = this->load_steps_num_ * solver_time_step;

    // Compute the displacement variation in each loading step.
    double load_step_displacement = 0.;
    for (auto load_step = 0; load_step != this->load_steps_num_; ++load_step) {

        // The normalized loading time step.
        double norm_tstep = ((load_step+1) * solver_time_step) / total_load_time;

        // Calculate loading step's displacement.
        load_step_displacement = this->max_displacement_ * (10.*norm_tstep*norm_tstep*norm_tstep -
                                                            15.*norm_tstep*norm_tstep*norm_tstep*norm_tstep +
                                                            6.*norm_tstep*norm_tstep*norm_tstep*norm_tstep*norm_tstep);

        // Store the loading step's displacement in the container.
        this->load_step_displacements_.emplace_back(load_step_displacement);

    }

}


const double & LoadCurve::LoadDispAt(const size_t &step) const
{
    // Check if displacement exists for the required time step.
    if (step > this->load_step_displacements_.size()) {
        std::string error = "[ExplicitSim ERROR] Cannot retreive displacement variation for the requested time step. "
                            "The requested step exceeds the size of the variation displacement container.";
        throw std::invalid_argument(error.c_str());
    }

    return this->load_step_displacements_[step];
}


bool LoadCurve::operator == (const LoadCurve &load_curve) const
{
    // Return true if the variables have equal values in both load curves.
    return ((this->load_time_ == load_curve.load_time_) &&
            (this->max_displacement_ == load_curve.max_displacement_) &&
            (this->load_steps_num_ == load_curve.load_steps_num_) &&
            (this->load_step_displacements_ == load_curve.load_step_displacements_) );
}


bool LoadCurve::operator != (const LoadCurve &load_curve) const
{
    // Return true if the variables have not equal values in both load curves.
    return !(*this == load_curve);
}


LoadCurve & LoadCurve::operator = (const LoadCurve &load_curve) {

    // Assign the variables of the given load curve.
    if (this != &load_curve) {
        this->load_time_ = load_curve.load_time_;
        this->max_displacement_ = load_curve.max_displacement_;
        this->load_steps_num_ = load_curve.load_steps_num_;
        this->load_step_displacements_ = load_curve.load_step_displacements_;
    }
    return *this;
}




} //end of namespace ExplicitSim


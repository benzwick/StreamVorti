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


#include "ExplicitSim/materials/neo_hookean.hpp"


namespace ExplicitSim {


NeoHookean::NeoHookean() : points_number_(-1)
{}


NeoHookean::~NeoHookean()
{}


void NeoHookean::SetPointsNumber(const int &points_number)
{
    //Set the number of points of the material.
    this->points_number_ = points_number;
}


void NeoHookean::SetDensity(const double &d_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the density container.
    this->density_.clear();

    // Assign the given density value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->density_.push_back(d_value);
    }
}


void NeoHookean::SetYoungModulus(const double &ym_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the Young modulus container.
    this->young_modulus_.clear();

    // Assign the given Young modulus value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->young_modulus_.push_back(ym_value);
    }
}


void NeoHookean::SetPoissonRatio(const double &pr_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the Poisson's ratio container.
    this->poisson_ratio_.clear();

    // Assign the given Poisson's ratio value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->poisson_ratio_.push_back(pr_value);
    }
}


void NeoHookean::SetBulkModulus(const double &bulk_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the Bulk modulus container.
    this->bulk_modulus_.clear();

    // Assign the given Bulk modulus value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->bulk_modulus_.push_back(bulk_value);
    }
}


void NeoHookean::SetLameLambda(const double &l_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the Lame lambda constant container.
    this->lambda_.clear();

    // Assign the given Lame lambda value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->lambda_.push_back(l_value);
    }

}


void NeoHookean::SetLameMu(const double &mu_value)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the Lame mu (shear modulus) constant container.
    this->mu_.clear();

    // Assign the given Lame mu (shear modulus) value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->mu_.push_back(mu_value);
    }

}


void NeoHookean::SetWaveSpeed(const double &wv_speed)
{
    // Check if number of material points has been assigned.
    if (this->points_number_ == -1) {
        std::string error = "ERROR: No number of points has been assigned to the material.";
        throw std::runtime_error(error.c_str());
    }

    // Clear the wave speed constant container.
    this->wave_speed_.clear();

    // Assign the given wave speed value to all the material points.
    for (int i = 0; i != this->points_number_; ++i) {
        this->wave_speed_.push_back(wv_speed);
    }

}


void NeoHookean::ComputeLameLambdaMu()
{
    // Check if Young modulus and Poisson's ratio are initialized and are equal.
    if ((this->young_modulus_.size() == 0) &&
            (this->young_modulus_.size() != this->poisson_ratio_.size()) ) {
        std::string error = "ERROR: Young modulus and Poisson's ratio containers are not consistently initialized.";
        throw std::runtime_error(error.c_str());
    }

    // Clear Lame constants containers.
    this->lambda_.clear();
    this->mu_.clear();

    // Calculate the lame constants.
    double lame_l = 0.; double lame_m = 0.;
    for (std::vector<double>::size_type i = 0; i != this->young_modulus_.size(); ++i) {
        // Lame lambda parameter.
        lame_l = (this->young_modulus_.at(i) * this->poisson_ratio_.at(i)) /
                ( (1. + this->poisson_ratio_.at(i)) * (1. - 2.*this->poisson_ratio_.at(i)) );

        this->lambda_.push_back(lame_l);

        // Lame mu parameter.
        lame_m = this->young_modulus_.at(i) / (2. * (1. + this->poisson_ratio_.at(i)));

        this->mu_.push_back(lame_m);
    }

}


void NeoHookean::ComputeBulkModulus()
{
    // Check if Young modulus and Poisson's ratio are initialized and are equal.
    if ((this->young_modulus_.size() == 0) &&
            (this->young_modulus_.size() != this->poisson_ratio_.size()) ) {
        std::string error = "ERROR: Young modulus and Poisson's ratio containers are not consistently initialized.";
        throw std::runtime_error(error.c_str());
    }

    // Clear Bulk modulus container.
    this->bulk_modulus_.clear();

    // Calculate the Bulk modulus.
    double bulk = 0.;
    for (std::vector<double>::size_type i = 0; i != this->young_modulus_.size(); ++i) {
        bulk = this->young_modulus_.at(i) / (3. * (1. - 2.*this->poisson_ratio_.at(i) ) );

        this->bulk_modulus_.push_back(bulk);

    }

}


void NeoHookean::ComputeWaveSpeed()
{
    // Check if Lame constants and density are initialized and are equal.
    if ((this->lambda_.size() == 0) ||
            (this->lambda_.size() != this->mu_.size()) ||
            (this->lambda_.size() != this->density_.size())) {
        std::string error = "[ExplicitSim ERROR] Lame constants (lambda, mu) and density containers are not consistently initialized.";
        throw std::runtime_error(error.c_str());
    }

    // Clear wave speed container.
    this->wave_speed_.clear();

    // Calculate the wave speed.
    double speed = 0.;
    for (std::vector<double>::size_type i = 0; i != this->young_modulus_.size(); ++i) {
        speed = std::sqrt((this->lambda_.at(i) + 2.*this->mu_.at(i)) / this->density_.at(i) );

        this->wave_speed_.push_back(speed);

    }

}


Eigen::Matrix3d NeoHookean::SpkStress(const Eigen::Matrix3d &FT, const int &integ_point_id) const
{

    // Determinant of deformation gradient.
    double det = std::abs(FT.determinant());

    // 3x3 Identity matrix.
    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity(3, 3);

    // Right Cauchy Green deformation tensor (Bathe P506).
    Eigen::Matrix3d C = Eigen::Matrix3d::Zero(3, 3);
    C.noalias() = FT * FT.transpose();

    // Inverse of the right Cauchy Green deformation tensor.
    Eigen::Matrix3d Cinv = C.inverse();

    return this->mu_[integ_point_id] * std::pow(det, -(2./3.)) * (identity - C.trace()/3.*Cinv) +
           this->bulk_modulus_[integ_point_id] * det * (det - 1.) * Cinv;

}


} //end of namespace ExplicitSim

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


#ifndef EXPLICITSIM_MATERIALS_NEO_HOOKEAN_HPP_
#define EXPLICITSIM_MATERIALS_NEO_HOOKEAN_HPP_

/*!
   \file neo_hookean.hpp
   \brief NeoHookean class header file.
   \author Konstantinos A. Mountris
   \date 20/05/2017
*/

#include <Eigen/Dense>

#include <vector>
#include <string>

#include <stdexcept>
#include <exception>

#include <cmath>


namespace ExplicitSim {

/*!
 *  \addtogroup Materials
 *  @{
 */


/*!
 * \class NeoHookean
 * \brief Class implemmenting a neo-hookean material for meshless models [strong-form/weak-form].
 */

class NeoHookean{
public:

    /*!
     * \brief NeoHookean constructor.
     */
    NeoHookean();


    /*!
     * \brief NeoHookean destructor.
     */
    virtual ~NeoHookean();


    /*!
     * \brief Set the number of points represented by the material.
     * \param [in] points_number The number of points.
     * \return [void]
     */
    void SetPointsNumber(const int &point_number);


    /*!
     * \brief Set the density of the material points.
     * \param [in] d_value The density value.
     * \return [void]
     */
    void SetDensity(const double &d_value);


    /*!
     * \brief Set the Young modulus of the material points.
     * \param [in] ym_value The Young modulus value.
     * \return [void]
     */
    void SetYoungModulus(const double &ym_value);


    /*!
     * \brief Set the Poisson's ratio of the material points.
     * \param [in] pr_value The Poisson's ratio value.
     * \return [void]
     */
    void SetPoissonRatio(const double &pr_value);


    /*!
     * \brief Set the Bulk modulus of the material points.
     * \param [in] bulk_value The Bulk modulus value.
     * \return [void]
     */
    void SetBulkModulus(const double &bulk_value);


    /*!
     * \brief Set the Lame lambda of the material points.
     * \param [in] l_value The Lame lambda value.
     * \return [void]
     */
    void SetLameLambda(const double &l_value);


    /*!
     * \brief Set the Lame mu (shear modulus) of the material points.
     * \param [in] mu_value The Lame mu (shear modulus) value.
     * \return [void]
     */
    void SetLameMu(const double &mu_value);


    /*!
     * \brief Set the wave speed of the material points.
     * \param [in] wv_speed The wave speed value.
     * \return [void]
     */
    void SetWaveSpeed(const double &wv_speed);


    /*!
     * \brief Compute the Lame lambda and mu (shear modulus) constants of the material points.
     * \return [void]
     */
    void ComputeLameLambdaMu();


    /*!
     * \brief Compute the Bulk modulus of the material points.
     * \return [void]
     */
    void ComputeBulkModulus();


    /*!
     * \brief Compute the wave speed of the material points.
     * \return [void]
     */
    void ComputeWaveSpeed();


    Eigen::Matrix3d SpkStress(const Eigen::Matrix3d &FT, const int &integ_point_id) const;


    /*!
     * \brief Get the number of points associated with the material.
     * \return [int] The number of points associated with the material.
     */
    inline const int & PointsNumber() { return this->points_number_; }


    /*!
     * \brief Get the density of the material points.
     * \return [std::vector<double>] The density of the material points.
     */
    inline const std::vector<double> Density() { return this->density_; }


    /*!
     * \brief Get the Young modulus of the material points.
     * \return [std::vector<double>] The Young modulus of the material points.
     */
    inline const std::vector<double> YoungModulus() { return this->young_modulus_; }


    /*!
     * \brief Get the Poisson's ratio of the material points.
     * \return [std::vector<double>] The Poisson's ratio of the material points.
     */
    inline const std::vector<double> PoissonRatio() { return this->poisson_ratio_; }


    /*!
     * \brief Get the Bulk modulus of the material points.
     * \return [std::vector<double>] The Bulk modulus of the material points.
     */
    inline const std::vector<double> BulkModulus() { return this->bulk_modulus_; }


    /*!
     * \brief Get the Lame lambda of the material points.
     * \return [std::vector<double>] The Lame lambda of the material points.
     */
    inline const std::vector<double> LameLambda() { return this->lambda_; }


    /*!
     * \brief Get the Lame mu of the material points.
     * \return [std::vector<double>] The Lame mu of the material points.
     */
    inline const std::vector<double> LameMu() { return this->mu_; }


    /*!
     * \brief Get the wave speed of the material points.
     * \return [std::vector<double>] The wave speed of the material points.
     */
    inline const std::vector<double> WaveSpeed() { return this->wave_speed_; }


private:
    int points_number_;                     /*!< The number of the material points. */

    std::vector<double> density_;           /*!< The density values of the material points. */

    std::vector<double> young_modulus_;     /*!< The Young modulus value of the material points. */

    std::vector<double> poisson_ratio_;     /*!< The Poisson's ratio value of the material points. */

    std::vector<double> bulk_modulus_;      /*!< The bulk modulus of the material points. */

    std::vector<double> lambda_;            /*!< The Lame lambda constant of the material points. */

    std::vector<double> mu_;                /*!< The shear modulus of the material points. */

    std::vector<double> wave_speed_;        /*!< The wave speed of the material points. */

};


/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_MATERIALS_NEO_HOOKEAN_HPP_

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



#ifndef EXPLICITSIM_CONDITIONS_LOAD_CURVE_HPP_
#define EXPLICITSIM_CONDITIONS_LOAD_CURVE_HPP_

/*!
   \file load_curve.hpp
   \brief LoadCurve class header file.
   \author Konstantinos A. Mountris
   \date 30/09/2017
*/


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
 * \class LoadCurve
 * \brief Class implemmenting a load curve to be used in load condition.
 *
 * The load curve applies the loading progressively, in load time steps, to maintain stability of the solution.
 *
 */

class LoadCurve {
public:

    /*!
     * \brief The LoadCurve constructor.
     *
     * Provides zero initialization for variables: load_time_, max_displacement_, and load_steps_num_.
     *
     */
    LoadCurve();


    /*!
     * \brief The LoadCurve destructor.
     */
    virtual ~LoadCurve();


    /*!
     * \brief Set the load curve's application time.
     * \param [in] load_time The load curve's application time.
     * \return [void]
     */
    inline void SetLoadTime(const double &load_time) { this->load_time_ = load_time; }


    /*!
     * \brief Set the maximum displacement to be applied by the load curve.
     * \param [in] max_displacement The maximum displacement to be applied by the load curve.
     * \return [void]
     */
    inline void SetMaxDisplacement(const double &max_displacement) { this->max_displacement_ = max_displacement; }


    /*!
     * \brief Compute the number of time steps required for the load curve's application.
     * \param [in] solver_time_step The time step of the explicit solver for stable solution.
     * \return [void]
     */
    void ComputeLoadStepsNum(const double &solver_time_step);


    /*!
     * \brief Compute the displacement variation in each loading time step.
     * \param [in] solver_time_step The time step of the explicit solver for stable solution.
     * \return [void]
     */
    void ComputeLoadStepDisplacements(const double &solver_time_step);


    /*!
     * \brief Get the load curve's application time.
     * \return [double] The load curve's application time.
     */
    inline const double & LoadTime() const { return this->load_time_; }


    /*!
     * \brief Get the maximum displacement that is applied by the load curve.
     * \return [double] The maximum displacement that is applied by the load curve.
     */
    inline const double & MaxDisplacement() const { return this->max_displacement_; }


    /*!
     * \brief Get the number of time steps required for the load curve's application.
     * \return [int] The number of time steps required for the load curve's application.
     */
    inline const int & LoadStepsNum() const { return this->load_steps_num_; }


    /*!
     * \brief Get the displacement variation in each loading time step.
     * \return [std::vector<double>] The displacement variation in each loading time step.
     */
    inline const std::vector<double> & LoadStepDisplacements() const { return this->load_step_displacements_; }


    /*!
     * \brief Get the displacement variation at a specific time step.
     * \param [in] step The step at which the displacement variation is requested.
     * \return [double] The displacement variation at a specific time step.
     */
    const double & LoadDispAt(const size_t &step) const;


    /*!
     * \brief Equality operator.
     * \param [in] load_curve The load curve to check for equality
     * \return [bool] True if the load curves are equal.
     */
    bool operator == (const LoadCurve &load_curve) const;


    /*!
     * \brief Non-equality operator.
     * \param [in] load_curve The load curve to check for non-equality
     * \return [bool] True if the load curves are not equal.
     */
    bool operator != (const LoadCurve &load_curve) const;


    /*!
     * \brief Assignment operator.
     * \param [in] load_curve the load curve to be assigned.
     * \return [ExplicitSim::LoadCurve] Returns the load curve for multiple assignment.
     */
    LoadCurve & operator = (const LoadCurve &load_curve);



private:
    double load_time_;                                  /*!< The load curve's application time. */

    double max_displacement_;                           /*!< The maximum displacement that is applied by the load curve. */

    int load_steps_num_;                                /*!< The number of time steps required for the load curve's application. */

    std::vector<double> load_step_displacements_;       /*!< The displacement variation in each loading time step. */

};




/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_CONDITIONS_LOAD_CURVE_HPP_

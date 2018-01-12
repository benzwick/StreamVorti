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



#ifndef EXPLICITSIM_CONDITIONS_EBCIEM_HPP_
#define EXPLICITSIM_CONDITIONS_EBCIEM_HPP_

/*!
   \file ebciem.hpp
   \brief Ebciem class header file.
   \author Konstantinos A. Mountris
   \date 30/09/2017
*/

#include "ExplicitSim/approximants/mmls_3d.hpp"

#include "ExplicitSim/models/weak_model_3d.hpp"

#include "ExplicitSim/support_domain/support_domain.hpp"

#include "ExplicitSim/conditions/dirichlet.hpp"
#include "ExplicitSim/conditions/loading.hpp"

#include "ExplicitSim/utilities/logger.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

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
 * \class Ebciem
 * \brief Class implemmenting the EBCIEM (Essential Boundary Condition Imposition in Explicit Meshless methods)
 *        for exact imposition of essential boundary conditions.
 *
 * \note Applicable when the used function approximants do not possess the Kronecker delta property (e.g. MMLS).
 *
 */

class Ebciem {

public:

    /*!
     * \brief The Ebciem default constructor.
     */
    Ebciem();


    /*!
     * \brief The Ebciem default destructor.
     */
    virtual ~Ebciem();


    /*!
     * \brief Set the properties of the nodal mmls approximants.
     * \param [in] basis_func_type The type of the basi function of the approximants.
     * \param [in] use_exact_derivatives Boolean to declare if exact derivatives should be used.
     * \return [void]
     */
    void SetNodalMmlsProperties(const std::string &basis_func_type, const bool &use_exact_derivatives);


    /*!
     * \brief Compute the correction matrices for fixed and displaced boundaries.
     * \param [in] model The model where the boundaries belong to.
     * \param [in] support The support domain for each node of the model.
     * \param [in] dirichlet_conds The fixed boundary conditions container.
     * \param [in] loading_conds The displaced boundary conditions container.
     * \param [in] is_simplified Boolean to declare is correction matrices should be computed with simplified EBCIEM.
     * \return [void]
     */
    void ComputeCorrectionMatrices(const WeakModel3D &model, const SupportDomain &support,
                                   const std::vector<Dirichlet> &dirichlet_conds,
                                   const std::vector<Loading> &loading_conds, bool is_simplified);


    /*!
     * \brief Get the nodal mmsl approximant.
     * \return [ExplicitSim::Mmls3d] The nodal mmsl approximant.
     */
    inline const Mmls3d & NodalMmls() const { return this->nodal_mmls_; }


    /*!
     * \brief Get the correction matrices of the fixed boundaries.
     * \return [std::vector<Eigen::SparseMatrix<double> >] The correction matrices of the fixed boundaries.
     */
    inline const std::vector<Eigen::SparseMatrix<double> > & DirichletCorrMats() const { return this->dirichlet_corr_mats_; }


    /*!
     * \brief Get the correction matrices of the displaced boundaries.
     * \return [std::vector<Eigen::SparseMatrix<double> >] The correction matrices of the displaced boundaries.
     */
    inline const std::vector<Eigen::SparseMatrix<double> > & LoadingCorrMats() const { return this->loading_corr_mats_; }



protected:

    /*!
     * \brief Create and get the correction matrix for a specific boundary using (s)implified EBCIEM.
     * \param [in] model_nodes_num The number of the model's geometry nodes.
     * \param [in] bound_nodes_ids The indices of the nodes belonging to the specific boundary.
     * \param [in] model_inv_mass_mat The inverse of the model's nodal mass matrix.
     * \return [void]
     */
    Eigen::SparseMatrix<double> SebciemCorrMat(const int &model_nodes_num, const std::vector<int> &bound_nodes_ids,
                                               const Eigen::SparseMatrix<double> &model_inv_mass_mat) const;



private:

    Mmls3d nodal_mmls_;                                             /*!< The mmls approximants of the nodes of the specified model. */

    std::vector<Eigen::SparseMatrix<double> > dirichlet_corr_mats_;     /*!< The correction matrices for all the fixed boundaries. */

    std::vector<Eigen::SparseMatrix<double> > loading_corr_mats_;      /*!< The correction matrices for all the displaced boundaries. */

};




/*! @} End of Doxygen Groups*/
} //end of namespace ExplicitSim

#endif //EXPLICITSIM_CONDITIONS_EBCIEM_HPP_


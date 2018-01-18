/*
 * StreamVorti - Software for solving PDEs using explicit methods.
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
 *      Konstantinos A. MOUNTRIS
 */

/*!
   \file dcpse_2d.hpp
   \brief Dcpse2d class header file.
   \author Konstantinos A. Mountris
   \date 12/01/2018
*/

#ifndef STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_
#define STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_


#include "StreamVorti/vectors/vectors.hpp"
#include "StreamVorti/elements/node.hpp"
#include "StreamVorti/utilities/logger.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <exception>

namespace StreamVorti {

/*!
 *  \addtogroup Approximants
 *  @{
 */


/*!
 * \class Dcpse2d
 * \brief Class implemmenting modified DCPSE function derivatives approximants in 2 dimensions.
 */

class Dcpse2d
{

public:
    /*!
     * \brief Dcpse2d constructor.
     */
    Dcpse2d();


    /*!
     * \brief Dcpse2d destructor.
     */
    virtual ~Dcpse2d();


    void ComputeDerivs(const std::vector<Node> &geom_nodes,
                       const std::vector<std::vector<int> > &support_nodes_ids,
                       const std::vector<double> &support_radiuses);


    void SaveDerivToFile(const std::string &deriv, const std::string &filename) const;


    inline const Eigen::SparseMatrix<double> & ShapeFunctionDx() const { return this->sh_func_dx_; }


    inline const Eigen::SparseMatrix<double> & ShapeFunctionDy() const { return this->sh_func_dy_; }


    inline const Eigen::SparseMatrix<double> & ShapeFunctionDxx() const { return this->sh_func_dxx_; }


    inline const Eigen::SparseMatrix<double> & ShapeFunctionDyy() const { return this->sh_func_dyy_; }


    inline const Eigen::SparseMatrix<double> & ShapeFunctionDxy() const { return this->sh_func_dxy_; }



private:

    Eigen::SparseMatrix<double> sh_func_dx_;     /*!< The shape function 1st x derivative matrix. */

    Eigen::SparseMatrix<double> sh_func_dy_;     /*!< The shape function 1st y derivative matrix. */

    Eigen::SparseMatrix<double> sh_func_dxx_;     /*!< The shape function 2nd x derivative matrix. */

    Eigen::SparseMatrix<double> sh_func_dyy_;     /*!< The shape function 2nd y derivative matrix. */

    Eigen::SparseMatrix<double> sh_func_dxy_;     /*!< The shape function 2nd xy derivative matrix. */
};



/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_APPROXIMANTS_DCPSE_2D_HPP_

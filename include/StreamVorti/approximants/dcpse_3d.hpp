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
   \file dcpse_3d.hpp
   \brief Dcpse3d class header file.
   \author Konstantinos A. Mountris
   \date 12/01/2018
*/

#ifndef STREAMVORTI_APPROXIMANTS_DCPSE_3D_HPP_
#define STREAMVORTI_APPROXIMANTS_DCPSE_3D_HPP_


#include "StreamVorti/approximants/dcpse.hpp"
#include "StreamVorti/vectors/vectors.hpp"
#include "StreamVorti/elements/node.hpp"
#include "StreamVorti/utilities/logger.hpp"

#include "mfem.hpp"

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
 * \class Dcpse3d
 * \brief Class implemmenting modified DCPSE function derivatives approximants in 2 dimensions.
 */

class Dcpse3d: public Dcpse
{

public:
    /*!
     * \brief Dcpse3d constructor.
     */
    Dcpse3d();


    /*!
     * \brief Dcpse3d constructor to match nodes of an MFEM H1 GridFunction.
     */
    Dcpse3d(mfem::GridFunction &gf,
            int CutoffRadAtNeighbor = 30,
            int SupportRadAtNeighbor = 5);


    /*!
     * \brief Dcpse3d destructor.
     */
    virtual ~Dcpse3d();


    void ComputeDerivs(const mfem::GridFunction &geom_nodes,
                       const std::vector<std::vector<int> > &support_nodes_ids,
                       const std::vector<double> &support_radiuses);


    void SaveDerivToFile(const std::string &deriv, const std::string &filename) const;


    const mfem::SparseMatrix & D(int i) const {
        switch (i)
        {
        case 0: return this->sh_func_dx_;
        case 1: return this->sh_func_dy_;
        case 2: return this->sh_func_dz_;
        }
        MFEM_ABORT( "Index " << i << " out of bounds." );
    }


    inline const mfem::SparseMatrix & ShapeFunctionDx() const { return this->sh_func_dx_; }


    inline const mfem::SparseMatrix & ShapeFunctionDy() const { return this->sh_func_dy_; }


    inline const mfem::SparseMatrix & ShapeFunctionDz() const { return this->sh_func_dz_; }


    inline const mfem::SparseMatrix & ShapeFunctionDxx() const { return this->sh_func_dxx_; }


    inline const mfem::SparseMatrix & ShapeFunctionDyy() const { return this->sh_func_dyy_; }


    inline const mfem::SparseMatrix & ShapeFunctionDzz() const { return this->sh_func_dzz_; }


    inline const mfem::SparseMatrix & ShapeFunctionDxy() const { return this->sh_func_dxy_; }


    inline const mfem::SparseMatrix & ShapeFunctionDxz() const { return this->sh_func_dxz_; }


    inline const mfem::SparseMatrix & ShapeFunctionDyz() const { return this->sh_func_dyz_; }


private:

    mfem::SparseMatrix sh_func_dx_;     /*!< The shape function 1st x derivative matrix. */

    mfem::SparseMatrix sh_func_dy_;     /*!< The shape function 1st y derivative matrix. */

    mfem::SparseMatrix sh_func_dz_;     /*!< The shape function 1st z derivative matrix. */

    mfem::SparseMatrix sh_func_dxx_;     /*!< The shape function 2nd x derivative matrix. */

    mfem::SparseMatrix sh_func_dyy_;     /*!< The shape function 2nd y derivative matrix. */

    mfem::SparseMatrix sh_func_dzz_;     /*!< The shape function 2nd z derivative matrix. */

    mfem::SparseMatrix sh_func_dxy_;     /*!< The shape function 2nd xy derivative matrix. */

    mfem::SparseMatrix sh_func_dxz_;     /*!< The shape function 2nd xz derivative matrix. */

    mfem::SparseMatrix sh_func_dyz_;     /*!< The shape function 2nd yz derivative matrix. */
};



/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_APPROXIMANTS_DCPSE_3D_HPP_

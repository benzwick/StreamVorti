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


#ifndef STREAMVORTI_MODELS_WEAK_MODEL_3D_HPP_
#define STREAMVORTI_MODELS_WEAK_MODEL_3D_HPP_

/*!
   \file strong_model_2d.hpp
   \brief StrongModel2d class header file.
   \author Konstantinos A. Mountris
   \date 12/01/2018
*/


#include "StreamVorti/support_domain/support_domain.hpp"
#include "StreamVorti/utilities/logger.hpp"
#include "StreamVorti/grid/grids.hpp"
#include "StreamVorti/mesh/mesh.hpp"

#include <Eigen/Dense>

#include <string>
#include <algorithm>

#include <stdexcept>
#include <exception>


namespace StreamVorti {

/*!
 *  \addtogroup Models
 *  @{
 */


/*!
 * \class StrongModel2d
 * \brief Class implemmenting a geometrical 2D model for strong-form meshless analysis.
 *
 */

class StrongModel2d {
public:
    /*!
     * \brief StrongModel2d constructor.
     */
    StrongModel2d();


    /*!
     * \brief StrongModel2d destructor.
     */
    virtual ~StrongModel2d();


    /*!
     * \brief Load the mesh representation of the model.
     * \param [in] mesh_filename The filename of the tetrahedral mesh to load.
     * \return [void]
     */
    void LoadMeshRepresentation(const std::string &mesh_filename);


    /*!
     * \brief Create the grid representation of the model.
     *
     * In order to create the grid representation of the model,
     * the mesh representation must has been already loaded.
     *
     * \return [void]
     */
    void CreateGridRepresentation();



    /*!
     * \brief Get the 3D grid representation of the model.
     * \return [ExplitSim::Grid3D] The model's 3D grid representation.
     */
    inline const Grid3D & Grid() const { return this->grid_; }




private:

    Grid3D grid_;                           /*!< The 3D grid representation of the model. */

};




/*! @} End of Doxygen Groups*/
} //end of namespace StreamVorti

#endif //STREAMVORTI_MODELS_WEAK_MODEL_3D_HPP_



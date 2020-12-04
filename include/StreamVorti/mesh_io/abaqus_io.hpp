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
   \file abaqus_io.hpp
   \brief Abaqus input/output class header file.
   \author Konstantinos A. Mountris
   \date 15/01/2018
*/

#ifndef STREAMVORTI_MESH_IO_ABAQUS_IO_HPP_
#define STREAMVORTI_MESH_IO_ABAQUS_IO_HPP_


#include "StreamVorti/vectors/vectors.hpp"
#include "StreamVorti/elements/node.hpp"
#include "StreamVorti/utilities/logger.hpp"

#include <vector>
#include <map>

#include <utility>
#include <cctype>

#include <exception>
#include <stdexcept>

#include <sstream>
#include <iostream>
#include <fstream>

#include <algorithm>
#include <string>

namespace StreamVorti {

/*!
 *  \addtogroup MeshIO
 *  @{
 */

/*!
 * \class AbaqusIO
 * \brief Class implemmenting input/output functionality for mesh in Abaqus format (.inp).
 *
 */
class AbaqusIO{
public:
    /*!
     * \brief AbaqusIO constructor.
     */
    AbaqusIO();


    /*!
     * \brief AbaqusIO destructor.
     */
    virtual ~AbaqusIO();


    /*!
     * \brief Load an abaqus mesh.
     *
     * The mesh to be loaded should be in abaqus format (.inp).
     *
     * \param [in] mesh_filename The filename (full path) of the mesh to be loaded.
     * \return [void]
     */
    void LoadMeshFrom(const std::string &mesh_filename);


    /*!
     * \brief Load the vertices of the readed mesh in the given vertices container.
     * \param [out] vertices The vertices container to load the vertices of the mesh.
     * \return [void] The vertices of the readed mesh.
     */
    void LoadNodesIn(std::vector<Node> &nodes);


private:
    std::vector<std::string> input_mesh_;                    /*!< The parsed mesh with an index per line for easy access. */

    int nodes_startline_;                                /*!< The index of the starting line of the nodes set in the mesh file. */

    std::vector<std::pair<int,int> > offsetted_nodes_;       /*!< The indices of the offsetted nodes from the storage order and their offsets. */
};


/*! @} End of Doxygen Groups*/

} //end of namespace StreamVorti

#endif //STREAMVORTI_MESH_IO_ABAQUS_IO_HPP_

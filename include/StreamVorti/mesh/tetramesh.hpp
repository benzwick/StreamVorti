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


/*!
   \file tetramesh.hpp
   \brief Tetramesh class header file.
   \author Konstantinos A. Mountris
   \date 10/05/2017
*/


#ifndef STREAMVORTI_MESH_TETRAMESH_HPP_
#define STREAMVORTI_MESH_TETRAMESH_HPP_

#include "StreamVorti/mesh_io/mesh_io.hpp"
#include "StreamVorti/mesh/mesh_properties.hpp"
#include "StreamVorti/vectors/vectors.hpp"
#include "StreamVorti/elements/elements.hpp"
#include "StreamVorti/sets/node_set.hpp"
#include "StreamVorti/utilities/logger.hpp"

#include <vector>

#include <exception>
#include <stdexcept>

namespace StreamVorti {

/*!
 *  \addtogroup Mesh
 *  @{
 */


/*!
 * \class TetraMesh
 * \brief Class implemmenting a tetrahedral mesh.
 */
class TetraMesh {
public:
    /*!
     * \brief TetraMesh constructor.
     */
    TetraMesh();


    /*!
     * \brief TetraMesh copy constructor.
     * \param [in] tetramesh The tetrahedral mesh to be copied.
     */
    TetraMesh(const TetraMesh &tetramesh);


    /*!
     * \brief TetraMesh destructor.
     */
    virtual ~TetraMesh();


    /*!
     * \brief Load a tetrahedral mesh.
     * \param [in] mesh_filename The filename (full path) for the file where the mesh should be loaded from.
     * \return [void]
     */
    void LoadFrom(const std::string &mesh_filename);


    /*!
     * \brief Save a tetrahedral mesh.
     * \param [in] mesh_filename The filename (full path) where the mesh should be saved.
     * \return [void]
     */
    void SaveTo(const std::string &mesh_filename);


    /*!
     * \brief Write access to the nodes of the mesh.
     * \return [std::vector<StreamVorti::Node>] the mesh nodes with write access.
     */
    inline std::vector<Node> & EditNodes() { return this->nodes_; }


    /*!
     * \brief Write access to the elements of the mesh.
     * \return [std::vector<StreamVorti::Tetrahedron>] The mesh elements with write access.
     */
    inline std::vector<Tetrahedron> & EditElements() { return this->tetras_; }


    /*!
     * \brief Read-only access to the nodes of the mesh.
     * \return [std::vector<StreamVorti::Node>] The mesh nodes with read-only access.
     */
    inline const std::vector<Node> & Nodes() const { return this->nodes_; }


    /*!
     * \brief Get the coordinates of the mesh's nodes.
     *
     * Performs an iteration through the nodes of the mesh when it is called and could
     * increase computational burden if called multiple times. In such a case it is prefered
     * to get the nodes of the mesh using TetraMesh::Nodes() and iterate through the nodes to get Node::Coordinates()
     *
     * \return [std::vector<StreamVorti::Vec3<double> >] The coordinates of the mesh's nodes.
     */
    inline std::vector<Vec3<double> > NodeCoordinates() const
    {
        std::vector<Vec3<double> > coordinates;
        for (auto &node : this->nodes_) {
            coordinates.emplace_back(node.Coordinates());
        }

        return coordinates;
    }


    /*!
     * \brief Get the number of nodes of the tetrahedral mesh.
     * \return [int] The number of nodes of the tetrahedral mesh.
     */
    inline const int NodesNum() const { return static_cast<int>(this->nodes_.size()); }


    /*!
     * \brief Get the sets of nodes of the mesh.
     *
     * These correspond to groups of nodes where a boundary condition could be specified.
     *
     * \return [std::vector<StreamVorti::NodeSet>] The sets of nodes of the mesh.
     */
    inline const std::vector<NodeSet> NodeSets() const { return this->node_sets_; }


    /*!
     * \brief Read-ony access to the elements of the mesh.
     * \return [std::vector<StreamVorti::Tetrahedron>] the mesh elements with read-only access.
     */
    inline const std::vector<Tetrahedron> & Elements() const { return this->tetras_; }


    /*!
     * \brief The type of the mesh (tetrahedral).
     * \return [StreamVorti::MeshType] the mesh type of the given mesh (tetrahedral).
     */
    inline const StreamVorti::MeshType & MeshType() const { return this->mesh_type_; }


    /*!
     * \brief Equal to operator.
     *
     * Compares tetrahedral meshes for equality.
     *
     * \param [in] tetramesh The tetrahedral mesh to compare.
     * \return [bool] TRUE if tetrahedral meshes are identical.
     */
    bool operator == (const TetraMesh &tetramesh) const;


    /*!
     * \brief Not equal to operator.
     *
     * Compares tetrahedral meshes for inequality.
     *
     * \param [in] tetramesh The tetrahedral mesh to compare.
     * \return [bool] TRUE if tetrahedral meshes are not identical.
     */
    bool operator != (const TetraMesh &tetramesh) const;


    /*!
     * \brief Assignment operator.
     *
     * Assigns all the properties of a given tetramesh (nodes, tetrahedra, mesh type).
     *
     * \param [in] tetramesh The tetrahedral mesh to assign.
     * \return [TetraMesh] The assigned tetrahedral mesh.
     */
    TetraMesh & operator = (const TetraMesh &tetramesh);

private:

    std::vector<Node> nodes_;                             /*!< The mesh nodes. */

    std::vector<NodeSet> node_sets_;            /*!< The mesh node sets. */

    std::vector<Tetrahedron> tetras_;                     /*!< The mesh tetrahedral elements. */

    StreamVorti::MeshType mesh_type_;                     /*!< The mesh type (tetrahedral). */
};


/*! @} End of Doxygen Groups*/
} //namespace StreamVorti

#endif //STREAMVORTI_MESH_TETRAMESH_HPP_

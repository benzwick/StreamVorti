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
   \file grid_2d.hpp
   \brief Grid2D class header file.
   \author Konstantinos A. Mountris
   \date 15/01/2018
*/

#ifndef STREAMVORTI_GRID_GRID_2D_HPP_
#define STREAMVORTI_GRID_GRID_2D_HPP_

#include "StreamVorti/elements/node.hpp"
#include "StreamVorti/mesh_io/mesh_io.hpp"

#include <vector>
#include <string>

namespace StreamVorti {

/*!
 *  \addtogroup Grid
 *  @{
 */


/*!
 * \class Grid2D
 * \brief Class implemmenting a 2D grid of nodes.
 */
class Grid2D {
public:
    /*!
     * \brief Grid2D constructor.
     */
    Grid2D();

    /*!
     * \brief Grid2D copy constructor.
     * \param [in] Grid2D The 2D grid to be copied.
     */
    Grid2D(const Grid2D &grid2D);


    /*!
     * \brief Grid2D destructor.
     */
    virtual ~Grid2D();


    void LoadFrom(const std::string &grid_filename);


    /*!
     * \brief Write access to the nodes of the mesh.
     * \return [std::vector<ExplicitSim::Node>] the mesh nodes with write access.
     */
    inline std::vector<Node> & EditNodes() { return this->nodes_; }


    /*!
     * \brief Read-only access to the nodes of the grid.
     * \return [std::vector<StreamVorti::Node>] the grid nodes with read-only access.
     */
    inline const std::vector<Node> & Nodes() const { return this->nodes_; }


    inline std::vector<Vec3<double> > NodeCoords() const
    {
        std::vector<Vec3<double> > coordinates;
        coordinates.reserve(this->nodes_.size());

        for (auto &node : this->nodes_) {
            coordinates.emplace_back(node.Coordinates());
        }

        return coordinates;
    }



    /*!
     * \brief Get the number of the grid's nodes.
     * \return [int] The number of the grid's nodes.
     */
    inline int NodesNum() const { return static_cast<int>(this->nodes_.size()); }


    /*!
     * \brief Equal to operator.
     *
     * Compares 2D grids for equality.
     *
     * \param [in] grid3D The 2D grid mesh to compare.
     * \return [bool] TRUE if 2D grids are identical.
     */
    bool operator == (const Grid2D &grid2D) const;


    /*!
     * \brief Not equal to operator.
     *
     * Compares 2D grids for inequality.
     *
     * \param [in] grid3D The 2D grids to compare.
     * \return [bool] TRUE if 2D grids are not identical.
     */
    bool operator != (const Grid2D &grid2D) const;


    /*!
     * \brief Assignment operator.
     *
     * Assigns all the properties of a given 3D grid (nodes).
     *
     * \param [in] grid3D The 2D grid to assign.
     * \return [Grid2D] The assigned 2D grid.
     */
    Grid2D & operator = (const Grid2D &grid2D);

private:
    std::vector<Node> nodes_;        /*!< The nodes of the grid. */
};



/*! @} End of Doxygen Groups*/
} // End of namespace StreamVorti

#endif //STREAMVORTI_GRID_GRID_2D_HPP_

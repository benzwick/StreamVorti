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
   \file node.hpp
   \brief Node class header file.
   \author Konstantinos A. Mountris
   \date 15/01/2018
*/

#ifndef STREAMVORTI_ELEMENTS_NODE_HPP_
#define STREAMVORTI_ELEMENTS_NODE_HPP_

#include "StreamVorti/vectors/vectors.hpp"


namespace StreamVorti {

/*!
 *  \addtogroup Elements
 *  @{
 */

/*!
 * \class Node
 * \brief Class implemmenting a node in 3 dimensions for meshfree applications.
 *
 * Node can be used to represent lower-dimension geometries too.
 *
 */
class Node {
public:
    /*!
     * \brief Node class constructor.
     */
    Node();


    /*!
     * \brief Node copy constructor.
     * \param [in] node The node to be copied on construction.
     */
    Node(const Node &node);


    /*!
     * \brief Node class destructor.
     */
    virtual ~Node();


    /*!
     * \brief Set the id of the node in a list of nodes.
     * \param id The id of the node.
     * \return [void]
     */
    void SetId(const int &id);


    /*!
     * \brief Sets the coordinates of the node.
     * \param [in] x The X coordinate.
     * \param [in] y The Y coordinate.
     * \param [in] z The Z coordinate.
     * \return [void]
     */
    void SetCoordinates(const double &x, const double &y, const double &z);


    /*!
     * \overload
     * \brief Sets the coordinates of the node.
     * \param [in] coords The X, Y, Z coordinates.
     * \return [void]
     */
    void SetCoordinates(const Vec3<double> &coords);


    /*!
     * \brief Sets the X coordinate of the node.
     * \param [in] x The X coordinate.
     * \return [void]
     */
    void SetCoordX(const double &x);


    /*!
     * \brief Sets the Y coordinate of the node.
     * \param [in] y The Y coordinate.
     * \return [void]
     */
    void SetCoordY(const double &y);


    /*!
     * \brief Sets the Z coordinate of the node.
     * \param [in] z The Z coordinate.
     * \return [void]
     */
    void SetCoordZ(const double &z);


    /*!
     * \brief Copy coordinates by a Vec3 vector.
     * \param [in] coordinates The Vec3 vector containing the coordinates to be copied.
     * \return [void]
     */
    void CopyCoordinates(const Vec3<double> &coordinates);


    /*!
     * \brief Get the Id of the node.
     * \return [int] The node's id, in other words its position in a list of nodes.
     */
    inline const int & Id() const { return this->id_; }


    /*!
     * \brief Get the node's coordinates.
     * \return [Vec3] The node's coordinates.
     */
    inline const Vec3<double> & Coordinates() const { return this->coordinates_; }


    /*!
     * \brief Equal to operator.
     *
     * Compares nodes for equality.
     *
     * \param [in] node The node to compare.
     * \return [bool] TRUE if nodes are identical.
     */
    bool operator == (const Node &node) const;


    /*!
     * \brief Not equal to operator.
     *
     * Compares nodes for inequality.
     *
     * \param [in] node The node to compare.
     * \return [bool] TRUE if nodes are not identical.
     */
    bool operator != (const Node &node) const;


    /*!
     * \brief Assignment operator.
     *
     * Assigns all the properties of a given node (id, coordinates).
     *
     * \param [in] node The node to assign.
     * \return [Node] The assigned node.
     */
    Node & operator = (const Node &node);


private:
    int id_;                            /*!< The node's id (position in a list of nodes). */

    Vec3<double> coordinates_;     /*!< The node's coordinates. */

};


/*! @} End of Doxygen Groups*/

} //end of namespace StreamVorti


#endif //STREAMVORTI_ELEMENTS_NODE_HPP_

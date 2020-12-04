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

#include "StreamVorti/elements/node.hpp"


namespace StreamVorti {


Node::Node() : id_(-1)
{
    //Id is initialized to -1 to indicate that node is not listed.

    // Initialize coordinates.
    this->coordinates_.Set(0., 0., 0.);
}


Node::Node(const Node &node)
{
    // Assign by copying node.
    *this = node;
}


Node::~Node()
{}


void Node::SetId(const int &id)
{
    // Set the nodes Id in a list of nodes.
    this->id_ = id;
}


void Node::SetCoordinates(const double &x, const double &y, const double &z)
{
    // Set the coordinates of the node.
    this->coordinates_.Set(x, y, z);
}


void Node::SetCoordinates(const Vec3<double> &coords)
{
    //Set the coordinates of the node.
    this->coordinates_ = coords;
}


void Node::SetCoordX(const double &x)
{
    // Set the node's x coordinate.
    this->coordinates_.SetX(x);
}


void Node::SetCoordY(const double &y)
{
    // Set the node's y coordinate.
    this->coordinates_.SetY(y);
}


void Node::SetCoordZ(const double &z)
{
    // Set the node's z coordinate.
    this->coordinates_.SetZ(z);
}


void Node::CopyCoordinates(const Vec3<double> &coordinates)
{
    // Set the node's coordinates by a copy.
    this->coordinates_ = coordinates;
}


bool Node::operator == (const Node &node) const
{
    // Compare nodes.
    return ((this->id_ == node.id_) &&
            (this->coordinates_ == node.coordinates_)
           );
}


bool Node::operator != (const Node &node) const
{
    // Compare nodes.
    return !(*this == node);
}


Node & Node::operator = (const Node &node)
{
    if (this != &node) {
        // Assign values from node.
        this->id_ = node.id_;
        this->coordinates_ = node.coordinates_;
    }

    return *this;

}


} //end of namespace StreamVorti



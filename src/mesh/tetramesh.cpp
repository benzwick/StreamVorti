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


#include "StreamVorti/mesh/tetramesh.hpp"

namespace StreamVorti {


TetraMesh::TetraMesh()
{}


TetraMesh::TetraMesh(const TetraMesh &tetramesh)
{
    // Copy by assigning tetramesh.
    *this = tetramesh;
}


TetraMesh::~TetraMesh()
{}


void TetraMesh::LoadFrom(const std::string &mesh_filename)
{
    // Check if mesh filename is not empty.
    if (mesh_filename.empty()) {
        throw std::invalid_argument(Logger::Error("Could not load mesh. No mesh filename was given.").c_str());
    }

    // Get the extension of the mesh filename.
    auto ext = mesh_filename.substr(mesh_filename.length()-4);

    // Clear the mesh containers.
    this->nodes_.clear();

    // Load the corresponding format.
    if (ext == ".inp") {
        AbaqusIO abaqus_io;
        abaqus_io.LoadMeshFrom(mesh_filename.c_str());
        abaqus_io.LoadNodesIn(this->nodes_);
    }
    else {
        std::string error = Logger::Error("Could not load mesh of unkown format. Expected [.inp | .feb] Check: ") + mesh_filename;
        throw std::invalid_argument(error.c_str());
    }

}


bool TetraMesh::operator == (const TetraMesh &tetramesh) const
{
    // Compare tetrahedral meshes for equality.
    return ((this->nodes_ == tetramesh.nodes_));
}


bool TetraMesh::operator != (const TetraMesh &tetramesh) const
{
    // Compare tetrahedral meshes for inequality.
    return !(*this == tetramesh);
}


TetraMesh & TetraMesh::operator = (const TetraMesh &tetramesh)
{
    if (this != &tetramesh) {
        // Assign values from tetrahedron.
        this->nodes_ = tetramesh.nodes_;
    }

    return *this;
}


}  //end of namespace StreamVorti

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


TetraMesh::TetraMesh() : mesh_type_(StreamVorti::MeshType::tetrahedral)
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
    this->tetras_.clear();
    this->node_sets_.clear();

    // Load the corresponding format.
    if (ext == ".inp") {
        AbaqusIO abaqus_io;
        abaqus_io.LoadMeshFrom(mesh_filename.c_str());
        abaqus_io.LoadNodesIn(this->nodes_);
        abaqus_io.LoadElementsIn(this->tetras_);
        if (abaqus_io.PartitionsExist()) {
            abaqus_io.LoadPartitionsIn(this->tetras_);
        }
        if (abaqus_io.NodeSetsExist()) {
            abaqus_io.LoadBoundarySetsIn(this->node_sets_);
        }

    }
    else {
        std::string error = Logger::Error("Could not load mesh of unkown format. Expected [.inp | .feb] Check: ") + mesh_filename;
        throw std::invalid_argument(error.c_str());
    }

}


void TetraMesh::SaveTo(const std::string &mesh_filename)
{
    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        std::string error = "ERROR: No filename was given to save the mesh.";
        throw std::invalid_argument(error.c_str());
    }

    // Get the extension of the mesh filename.
    auto ext = mesh_filename.substr(mesh_filename.length()-4);


    if (ext == ".inp") {
        AbaqusIO abaqus_io;
        abaqus_io.SaveMesh<TetraMesh,Tetrahedron>(*this, mesh_filename.c_str());
    }
    else {
        std::string error = "ERROR: Given mesh file: \"" + mesh_filename + "\" is of unknown format.";
        throw std::invalid_argument(error.c_str());
    }
}


bool TetraMesh::operator == (const TetraMesh &tetramesh) const
{
    // Compare tetrahedral meshes for equality.
    return ((this->nodes_ == tetramesh.nodes_) &&
            (this->node_sets_ == tetramesh.node_sets_) &&
            (this->tetras_ == tetramesh.tetras_) &&
            (this->mesh_type_ == tetramesh.mesh_type_)
           );
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
        this->node_sets_ = tetramesh.node_sets_;
        this->tetras_ = tetramesh.tetras_;
        this->mesh_type_ = tetramesh.mesh_type_;
    }

    return *this;
}


}  //end of namespace StreamVorti

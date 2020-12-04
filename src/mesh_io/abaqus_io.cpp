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



#include "StreamVorti/mesh_io/abaqus_io.hpp"

namespace StreamVorti {


AbaqusIO::AbaqusIO() : nodes_startline_(0)
{}


AbaqusIO::~AbaqusIO()
{}


void AbaqusIO::LoadMeshFrom(const std::string &mesh_filename)
{
    // Print reading status message.
    std::cout << Logger::Message("Reading mesh file: ") << mesh_filename << "\n";

    // Clear containers.
    this->input_mesh_.clear();

    // Check if mesh filename is given.
    if (mesh_filename.empty()) {
        std::string error = Logger::Error("No filename was given to read mesh.");
        throw std::invalid_argument(error.c_str());
    }

    // Check if mesh is in abaqus format (.inp).
    std::string ext = mesh_filename.substr(mesh_filename.length()-4);
    if (ext != ".inp") {
        std::string error = Logger::Error("The file: ") + mesh_filename + " is not in Abaqus format (.inp)";
        throw std::invalid_argument(error.c_str());
    }

    //Open mesh file.
    std::ifstream mesh_file(mesh_filename.c_str(), std::ios::in);

    // Check if mesh file opened successfully.
    if (!mesh_file.is_open()) {
        std::string error = Logger::Error("Could not open the file: ") + mesh_filename + " Check given path.";
        throw std::runtime_error(error.c_str());
    }

    //Available mesh information.
    bool nodes_exist = false;

    //Line of the file.
    std::string line = "";
    int line_id = 0;

    //Read mesh file line by line.
    while (std::getline(mesh_file, line)) {
        // Transform mesh file line in lowercase.
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // Read mesh file lines.
        this->input_mesh_.push_back(line);

        // Find the first line of the nodes set.
        if (line.find("*node") != std::string::npos) {
            this->nodes_startline_ = line_id;
            nodes_exist = true;
        }

        // Increase line id once it is processed.
        line_id++;

    }

    // Check if nodes are available in the mesh file.
    if (!nodes_exist) {
        std::string error = "[ExplicitSim ERROR] Mesh file: " + mesh_filename + " is incomplete. "
                "Nodes are not available...";
        std::runtime_error(error.c_str());
    }

    // Close the mesh file.
    mesh_file.close();

}


void AbaqusIO::LoadNodesIn(std::vector<Node> &nodes)
{
    // Clean the mesh nodes container.
    nodes.clear();

    // The coordinates of the nodes in the mesh.
    double x = 0.; double y = 0.;

    // The id of the nodes in the mesh. Initialized to invalid value (-1).
    int id = -1;

    // The mesh node.
    Node node;

    // Iterate though mesh starting from the nodes set starting line.
    // Skip the first line start from the first node's coordinates.
    std::string line = "";
    for (std::vector<std::string>::size_type it = this->nodes_startline_ + 1;
         it != input_mesh_.size(); ++it) {

        line = this->input_mesh_.at(it);

        // Replace comma occurence in line with space.
        std::replace(line.begin(), line.end(), ',', ' ');

        // Convert line to istringstream to be passed in the coordinates variables.
        std::stringstream ss(line);

        // Get the coordinates until reach the end of the vertices set.
        if (!(ss >> id >> x >> y)) { break; }

        // Set the id and coordinates of the mesh node. Reduce index by 1 to account for the storage offset.
        node.SetId(id-1);
        node.SetCoordinates(x, y, 0.);

        // Store the node in the nodes' container.
        nodes.emplace_back(node);
    }

    if (nodes.empty()) {
        throw std::runtime_error(Logger::Error("Could not load nodes from Abaqus. "
                                               "Check nodes in the given filename.").c_str());
    }

    // Reset the offsetted nodes indices, if any.
    int node_offset = 0;
    for (auto &node : nodes) {
        // Index of the node in the container.
        auto id = &node - &nodes[0];

        // Check consistency with stored node's id.
        if (id != node.Id()) {
            node_offset = node.Id() - id;

            // Store the offsetted node's id and it's offset from storage.
            this->offsetted_nodes_.push_back(std::make_pair(node.Id(), node_offset));

            // Reset the node's id.
            node.SetId(node.Id() - node_offset);
        }
    }


}

} // end of namespace StreamVorti

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


#include "StreamVorti/options_configuration/config_manager.hpp"


namespace StreamVorti {


ConfigManager::ConfigManager()
{
    this->description_.add_options()
            ("Model.GridFile", boost_po::value<std::string>()->required(),
             "Absolute path of model's grid file. Supported formats: [.inp]")

            ("DCPSE.CutoffRadAtNeighbor", boost_po::value<int>()->default_value(30),
             "Index of N nearest neighbor node to set cutoff radius.")
            ("DCPSE.SupportRadAtNeighbor", boost_po::value<int>()->default_value(5),
             "Index of N nearest neighbor node to set support radius.")
            ("DCPSE.SaveNeighborsToFile", boost_po::value<std::string>()->default_value(""),
             "Save the neighbor list for each node to file [.txt].")
            ("DCPSE.SaveDxToFile", boost_po::value<std::string>()->default_value(""),
             "Save first x derivative sparse matrix to file [.txt].")
            ("DCPSE.SaveDyToFile", boost_po::value<std::string>()->default_value(""),
             "Save first y derivative sparse matrix to file [.txt].")
            ("DCPSE.SaveDxxToFile", boost_po::value<std::string>()->default_value(""),
             "Save second xx derivative sparse matrix to file [.txt].")
            ("DCPSE.SaveDyyToFile", boost_po::value<std::string>()->default_value(""),
             "Save second yy derivative sparse matrix to file [.txt].")
            ("DCPSE.SaveDxyToFile", boost_po::value<std::string>()->default_value(""),
             "Save second xy derivative sparse matrix to file [.txt].");

}


ConfigManager::~ConfigManager()
{}



void ConfigManager::ReadConfigFile(const std::string &config_filename)
{
    std::ifstream input_stream(config_filename.c_str(), std::ifstream::in);
    if (!input_stream.is_open()) {
        std::string error = "[StreamVorti ERROR] cannot open StreamVorti configuration file: " + config_filename;
        throw std::invalid_argument(error.c_str());
    }

    // Store the parsed configuration file.
    boost_po::store(boost_po::parse_config_file(input_stream, this->description_), this->var_map_);
    input_stream.close();
    boost_po::notify(this->var_map_);
}


void ConfigManager::ReadConfigStream(std::istream &config_stream)
{
    boost_po::store(boost_po::parse_config_file(config_stream, this->description_), this->var_map_);
    boost_po::notify(this->var_map_);
}


std::string ConfigManager::PrintSampleFile()
{
    std::string sample = "#\n"
            "# StreamVorti\n"
            "# Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>\n"
            "#\n"
            "# This program is free software: you can redistribute it and/or modify\n"
            "# it under the terms of the GNU General Public License as published by\n"
            "# the Free Software Foundation, either version 3 of the License, or\n"
            "# (at your option) any later version.\n"
            "#\n"
            "# This program is distributed in the hope that it will be useful,\n"
            "# but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
            "# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
            "# GNU General Public License for more details.\n"
            "#\n"
            "# You should have received a copy of the GNU General Public License\n"
            "# along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
            "#\n"
            "# Contributors (alphabetically):\n"
            "#      George C. BOURANTAS\n"
            "#      Konstantinos A. MOUNTRIS\n"
            "#\n\n\n"
            "######                  StreamVorti V0.0.1 configuration file sample.                  ######\n"
            "\n"
            "# Section: Introduction\n"
            "# ---------------------\n"
            "\n"
            "# This text file demonstrates all the available arguments to launch and execute a StreamVorti\n"
            "# simulation. These arguments are grouped in sections following a logic order in this sample,\n"
            "# but they can be given in any order in general. Other than the mandatory arguments [stated later],\n"
            "# the rest of the arguments can be ommitted and the default variables used in this sample will be used.\n"
            "\n\n"
            "[Model]                                                 # Section: Model\n"
            "                                                        # --------------\n"
            "\n"
            "GridFile = /path/to/mesh/file                           # For now only grid files in Abaqus\n"
            "                                                        # format (.inp) are supported. [mandatory]\n"
            "\n\n"
            "[DCPSE]                                                 # Section: DCPSE\n"
            "                                                        # --------------\n"
            "\n"
            "CutoffRadAtNeighbor = 30                                # Index of N nearest neighbor\n"
            "                                                        # node to set cutoff radius.\n"
            "\n"
            "SupportRadAtNeighbor = 5                                # Index of N nearest neighbor\n"
            "                                                        # node to set cutoff radius.\n"
            "\n"
            "SaveNeighborsToFile =                                   # Save the neighbor list for each node to file [.txt].\n"
            "                                                        # Format: [model nodes x neigh_indices]\n"
            "\n"
            "SaveDxToFile =                                          # Save first x derivative sparse matrix to file [.txt].\n"
            "\n"
            "SaveDyToFile =                                          # Save first y derivative sparse matrix to file [.txt].\n"
            "\n"
            "SaveDxxToFile =                                         # Save second xx derivative sparse matrix to file [.txt].\n"
            "\n"
            "SaveDyyToFile =                                         # Save second yy derivative sparse matrix to file [.txt].\n"
            "\n"
            "SaveDxyToFile =                                         # Save second xy derivative sparse matrix to file [.txt].\n"
            "\n\n";

    return sample;

}


} // End of namespace StreamVorti

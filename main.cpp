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
 *      Konstantinos A. MOUNTRIS
 */

#include <StreamVorti/stream_vorti.hpp>

#include <boost/filesystem.hpp>

#include <cstddef>
#include <string>
#include <fstream>

using namespace StreamVorti;

int main(int argc, char *argv[]) {

    try {

        // Initialize configuration manager.
        ConfigManager *config = new ConfigManager();

        // Initialize configuration filename empty to be read during the StreamVorti execution.
        std::string config_filename = "";

        // Check if configuration file was provided during execution.
        if (argc == 1) {
            std::cout << Logger::Warning("No input file was specified. Type configuration filename with absolute path.\n"
                         "Otherwise tap '-g' to generate sample configuration file or '-q' to exit StreamVorti.\nInput filename: ");

            // Read configuration file given by the user.
            std::cin >> config_filename;

            if (config_filename == "-g") {
                std::cout << Logger::Message("Give text file name [.ini] with absolute path to store the sample configuration "
                             "file\nor tap '-t' to print in terminal.\nSample filename: ");

                std::string sample_filename = "";

                std::cin >> sample_filename;

                if (sample_filename == "-t") {
                    std::cout << config->PrintSampleFile() << std::endl;
                    std::cout << Logger::Message("Save the sample configuration in a text file [.ini], edit it according to your simulation,"
                                 " and relaunch StreamVorti passing your configuration file as argument.\n");
                    return EXIT_SUCCESS;
                }
                else {
                    // Initialize the path of the sample configuration file.
                    std::string sample_path = "";

                    // Position of the last slash in the sample configuration file.
                    std::size_t last_slash = sample_filename.find_last_of("/\\");

                    // Get the path directory of the sample configuration file.
                    if (last_slash != std::string::npos) {
                        sample_path = sample_filename.substr(0, last_slash);
                    }

                    // Create the path's directory if it doesn't exist.
                    boost::filesystem::path dir(sample_path);
                    if (!sample_path.empty() && !boost::filesystem::exists(dir)) {
                        boost::filesystem::create_directories(dir);
                    }

                    // Position of the last slash in the exporting file's name.
                    std::size_t last_dot = sample_filename.find_last_of(".");

                    // Initialize sample configuration file extension.
                    std::string sample_ext = "";

                    // Get sample configuration file extension.
                    if (last_dot != std::string::npos) {
                        sample_ext = sample_filename.substr(last_dot);
                    }

                    // Add extension if necessary
                    if (sample_ext != ".ini") { sample_filename += ".ini"; }

                    // Output the sample file.
                    std::ofstream sample_output(sample_filename, std::ios_base::out | std::ios_base::trunc);
                    sample_output << config->PrintSampleFile();
                    std::cout << Logger::Message("Sample configuration file saved at: ") << sample_filename << std::endl;
                    std::cout << Logger::Message("Edit the sample configuration according to your simulation and "
                                 "relaunch StreamVorti passing your configuration file as argument.\n");
                    return EXIT_SUCCESS;
                }

            }

            if (config_filename == "-q") { std::cout << Logger::Message("User requested termination. See you soon!\n"); exit(0); }
            std::cout << config_filename << std::endl;
            exit(0);
        }
        else { config_filename = argv[1]; }

        // Read configuration file.
        config->ReadConfigFile(config_filename);

        std::cout << "\t<<< Welcome to StreamVorti >>>\n";
        std::cout << Logger::Message("Loading configuration file: ") << config_filename << std::endl;

        // Profiling spent time in StreamVorti
        Timer timer;

        // Create 2d strong model.
        StrongModel2d model;
        model.LoadGrid(config->RetrieveArgument<std::string>("Model.GridFile"));
        std::cout << Logger::Message("Model has nodes: ") << model.Grid().Nodes().size() << std::endl;

        // Set support domain.
        SupportDomain support;
        support.SetSupportNodes(model.Grid().Nodes());

        timer.Reset();
        support.ComputeCutOffRadiuses(config->RetrieveArgument<int>("DCPSE.CutoffRadAtNeighbor"));
        std::cout << Logger::Message("Execution time for cut-off radiuses computation for all nodes: ")
                  << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        support.ComputeSupportRadiuses(config->RetrieveArgument<int>("DCPSE.SupportRadAtNeighbor"));
        std::cout << Logger::Message("Execution time for support radiuses computation for all nodes: ")
                  << timer.PrintElapsedTime() << "\n";

        timer.Reset();
        auto neighs = support.NeighborIndices();
        std::cout << Logger::Message("Execution time for neighbor indices computation for all nodes: ")
                  << timer.PrintElapsedTime() << "\n";

        if (config->RetrieveArgument<std::string>("DCPSE.SaveNeighborsToFile") != "") {
            support.SaveNeighsToFile(neighs, config->RetrieveArgument<std::string>("DCPSE.SaveNeighborsToFile"));
        }

        Dcpse2d derivs;
        timer.Reset();
        derivs.ComputeDerivs(model.Grid().Nodes(), neighs, support.SupportRadiuses());
        std::cout << Logger::Message("Execution time for DCPSE derivatives calculation: ")
                  << timer.PrintElapsedTime() << "\n";

        if (config->RetrieveArgument<std::string>("DCPSE.SaveDxToFile") != "") {
            derivs.SaveDerivToFile("dx", config->RetrieveArgument<std::string>("DCPSE.SaveDxToFile"));
        }

        if (config->RetrieveArgument<std::string>("DCPSE.SaveDyToFile") != "") {
            derivs.SaveDerivToFile("dy", config->RetrieveArgument<std::string>("DCPSE.SaveDyToFile"));
        }

        if (config->RetrieveArgument<std::string>("DCPSE.SaveDxxToFile") != "") {
            derivs.SaveDerivToFile("dxx", config->RetrieveArgument<std::string>("DCPSE.SaveDxxToFile"));
        }

        if (config->RetrieveArgument<std::string>("DCPSE.SaveDyyToFile") != "") {
            derivs.SaveDerivToFile("dyy", config->RetrieveArgument<std::string>("DCPSE.SaveDyyToFile"));
        }

        if (config->RetrieveArgument<std::string>("DCPSE.SaveDxyToFile") != "") {
            derivs.SaveDerivToFile("dxy", config->RetrieveArgument<std::string>("DCPSE.SaveDxyToFile"));
        }

        std::cout << Logger::Message("Simulation terminated successfully.") << std::endl;

        // Release memory from configuration manager.
        delete config;

    }
    catch (const std::invalid_argument &e) {
        std::cerr << e.what() << std::endl;
    }
    catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
    }
    catch (const std::out_of_range &e) {
       std::cerr << e.what() << std::endl;
    }
    catch (const std::bad_alloc &e) {
       std::cerr << e.what() << std::endl;
    }
    catch (const boost::program_options::error &e) {
        std::cerr << "[Boost program_options error] " << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "[ExplicitSim Unknown exception]" << std::endl;
    }

    return EXIT_SUCCESS;
}

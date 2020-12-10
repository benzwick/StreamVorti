
#include <StreamVorti/stream_vorti.hpp>

#include <filesystem>
#include <cstddef>
#include <string>
#include <fstream>

#include "mfem.hpp"
#include <iostream>


using namespace StreamVorti;

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // Mesh mesh(10, 10, Element::QUADRILATERAL);
    Mesh mesh(10, 10, Element::QUADRILATERAL, false, 1.0, 1.0, false);

    ofstream mesh_ofs("mfem_square10x10.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);

    const int dim = mesh.Dimension();
    const int order = 1;
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec, dim);

    GridFunction nodes(&fes);
    mesh.GetNodes(nodes);

    // cout << nodes;

    // const int nNodes = nodes.Size() / dim;
    // double coord[dim]; // coordinates of a node
    // for (int i = 0; i < nNodes; ++i)
    // {
    //     for (int j = 0; j < dim; ++j)
    //     {
    //         coord[j] = nodes(j * nNodes + i);
    //         cout << coord[j] << " ";
    //     }
    //     cout << endl;
    // }

    cout << "Initialize configuration manager." << endl;
    ConfigManager *config = new ConfigManager();

    cout << "Initialize configuration filename empty to be read during the StreamVorti execution." << endl;
    std::string config_filename = "";

    cout << "Check if configuration file was provided during execution." << endl;
    if (argc == 1) {
        std::cout << "WARNING: No input file was specified. "
                  << "Type configuration filename with absolute path." << std::endl
                  << "Otherwise tap '-g' to generate sample configuration file "
                  << "or '-q' to exit StreamVorti.\nInput filename: ";

        // Read configuration file given by the user.
        std::cin >> config_filename;

        if (config_filename == "-g") {
            std::cout << "Give text file name [.ini] with absolute path to store the sample configuration "
                      << "file\nor tap '-t' to print in terminal.\nSample filename: " << std::endl;

            std::string sample_filename = "";

            std::cin >> sample_filename;

            if (sample_filename == "-t") {
                std::cout << config->PrintSampleFile() << std::endl;
                std::cout << "Save the sample configuration in a text file [.ini], edit it according to your simulation,"
                          << " and relaunch StreamVorti passing your configuration file as argument." << std::endl;
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
                std::filesystem::path dir(sample_path);
                if (!sample_path.empty() && !std::filesystem::exists(dir)) {
                    std::filesystem::create_directories(dir);
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
                std::cout << "Sample configuration file saved at: " << sample_filename << std::endl;
                std::cout << "Edit the sample configuration according to your simulation and "
                          << "relaunch StreamVorti passing your configuration file as argument." << std::endl;
                return EXIT_SUCCESS;
            }

        }

        if (config_filename == "-q") {
            std::cout << "User requested termination. See you soon!\n" << std::endl;
            exit(0);
        }
        std::cout << config_filename << std::endl;
        exit(0);
    }
    else { config_filename = argv[1]; }

    cout << "Read configuration file." << endl;
    config->ReadConfigFile(config_filename);

    std::cout << "\t<<< Welcome to StreamVorti >>>\n";
    std::cout << "Loading configuration file: " << config_filename << std::endl;

    // Profiling spent time in StreamVorti
    mfem::StopWatch timer;

    cout << "Set support domain." << endl;
    SupportDomain support(nodes);

    timer.Start();
    cout << "support: compute cutoff radiuses" << endl;
    support.ComputeCutOffRadiuses(config->RetrieveArgument<int>("DCPSE.CutoffRadAtNeighbor"));
    std::cout << "Execution time for cut-off radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    cout << "support: compute support radiuses" << endl;
    support.ComputeSupportRadiuses(config->RetrieveArgument<int>("DCPSE.SupportRadAtNeighbor"));
    std::cout << "Execution time for support radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    cout << "support: compute neighbor indices" << endl;
    auto neighs = support.NeighborIndices();
    std::cout << "Execution time for neighbor indices computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    if (config->RetrieveArgument<std::string>("DCPSE.SaveNeighborsToFile") != "") {
        cout << "support: save neighbor indices to file" << endl;
        support.SaveNeighsToFile(neighs, config->RetrieveArgument<std::string>("DCPSE.SaveNeighborsToFile"));
    }

    cout << "DC PSE derivatives." << endl;
    Dcpse2d derivs(nodes);
    timer.Clear();
    derivs.Update();
    // derivs.ComputeDerivs(nodes, neighs, support.SupportRadiuses());
    std::cout << "Execution time for DCPSE derivatives calculation: "
              << timer.RealTime() << " s" << std::endl;

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

    std::cout << "Simulation terminated successfully." << std::endl;

    cout << "Release memory from configuration manager." << endl;
    delete config;

    return EXIT_SUCCESS;
}

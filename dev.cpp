#include <StreamVorti/stream_vorti.hpp>

#include <cstddef>
#include <string>
#include <fstream>
#include <algorithm>

using namespace StreamVorti;

int main()
{
    try {

       std::cout << "\t<<< Welcome to StreamVorti >>>\n";

       // Load grid.
       Grid grid;
       grid.LoadFrom("/home/mood/DATABASE/mesh/geometric/2D/square10x10.inp");

       // Set support domain.
       SupportDomain support;
       support.SetSupportNodes(grid.Nodes());
       support.ComputeCutOffRadiuses(30);
       support.ComputeSupportRadiuses(5);
       auto neighs = support.NeighborIndices();


       Dcpse2d derivs;
       derivs.ComputeDerivs(grid.Nodes(), neighs, support.SupportRadiuses());


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
    catch (...) {
        std::cerr << "[StreamVorti Unknown exception]" << std::endl;
    }

    return EXIT_SUCCESS;
}

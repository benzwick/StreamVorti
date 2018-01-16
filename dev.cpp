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

       // Create 2d strong model.
       StrongModel2d model;
       model.LoadGrid("/home/mood/DATABASE/mesh/geometric/2D/square10x10.inp");

       // Set support domain.
       SupportDomain support;
       support.SetSupportNodes(model.Grid().Nodes());
       support.ComputeCutOffRadiuses(30);
       support.ComputeSupportRadiuses(5);
       auto neighs = support.NeighborIndices();

       for (auto &o_it : neighs) {
           std::sort(o_it.begin(), o_it.end());
           for (auto &i_it : o_it) { std::cout << i_it+1 << " "; }
           std::cout << std::endl << std::endl;
       }


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

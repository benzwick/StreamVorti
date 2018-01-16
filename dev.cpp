#include <StreamVorti/stream_vorti.hpp>

#include <cstddef>
#include <string>
#include <fstream>

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

       for (auto &radius : support.SupportRadiuses()) {
           std::cout << std::setprecision(15) << radius << std::endl;
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

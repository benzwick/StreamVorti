#include <StreamVorti/stream_vorti.hpp>

#include <cstddef>
#include <string>
#include <fstream>

using namespace StreamVorti;

int main()
{
    try {

       std::cout << "<<< Welcome to StreamVorti >>>" << std::endl;

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
        std::cerr << "[ExplicitSim Unknown exception]" << std::endl;
    }

    return EXIT_SUCCESS;
}

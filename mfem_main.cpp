
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
    // Output filename prefix and extension
    std::string fname = "mfem_square10x10";
    std::string fext = ".dat";

    // DC PSE parameters
    int CutoffRadAtNeighbor = 30;
    int SupportRadAtNeighbor = 5;

    // Output requests
    bool save_neighbors = true;
    bool save_dx = true;
    bool save_dy = true;
    bool save_dz = true;
    bool save_dxx = true;
    bool save_dxy = true;
    bool save_dyy = true;
    bool save_dyz = true;
    bool save_dzz = true;

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

    // Profiling spent time in StreamVorti
    mfem::StopWatch timer;

    cout << "Set support domain." << endl;
    SupportDomain support(nodes);

    timer.Start();
    cout << "support: compute cutoff radiuses" << endl;
    support.ComputeCutOffRadiuses(CutoffRadAtNeighbor);
    std::cout << "Execution time for cut-off radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    cout << "support: compute support radiuses" << endl;
    support.ComputeSupportRadiuses(SupportRadAtNeighbor);
    std::cout << "Execution time for support radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    cout << "support: compute neighbor indices" << endl;
    auto neighs = support.NeighborIndices();
    std::cout << "Execution time for neighbor indices computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    if (save_neighbors)
    {
        cout << "support: save neighbor indices to file" << endl;
        support.SaveNeighsToFile(neighs, fname + ".neighbors" + fext);
    }

    cout << "DC PSE derivatives." << endl;
    Dcpse2d derivs(nodes);
    timer.Clear();
    derivs.Update();
    std::cout << "Execution time for DCPSE derivatives calculation: "
              << timer.RealTime() << " s" << std::endl;

    if (save_dx)  {derivs.SaveDerivToFile("dx", fname + ".dx" + fext);}
    if (save_dy)  {derivs.SaveDerivToFile("dy", fname + ".dy" + fext);}
    if (save_dxx) {derivs.SaveDerivToFile("dxx", fname + ".dxx" + fext);}
    if (save_dxy) {derivs.SaveDerivToFile("dxy", fname + ".dxy" + fext);}
    if (save_dyy) {derivs.SaveDerivToFile("dyy", fname + ".dyy" + fext);}

    std::cout << "Simulation terminated successfully." << std::endl;

    return EXIT_SUCCESS;
}

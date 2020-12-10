
#include <StreamVorti/stream_vorti.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "mfem.hpp"

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
    mfem::Mesh mesh(10, 10, mfem::Element::QUADRILATERAL, false, 1.0, 1.0, false);

    std::ofstream mesh_ofs("mfem_square10x10.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);

    const int dim = mesh.Dimension();
    const int order = 1;
    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fes(&mesh, &fec, 1);

    mfem::GridFunction gf(&fes);

    // Profiling spent time in StreamVorti
    mfem::StopWatch timer;

    std::cout << "Set support domain." << std::endl;
    StreamVorti::SupportDomain support(gf);

    timer.Start();
    std::cout << "support: compute cutoff radiuses" << std::endl;
    support.ComputeCutOffRadiuses(CutoffRadAtNeighbor);
    std::cout << "Execution time for cut-off radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    std::cout << "support: compute support radiuses" << std::endl;
    support.ComputeSupportRadiuses(SupportRadAtNeighbor);
    std::cout << "Execution time for support radiuses computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    std::cout << "support: compute neighbor indices" << std::endl;
    auto neighs = support.NeighborIndices();
    std::cout << "Execution time for neighbor indices computation for all nodes: "
              << timer.RealTime() << " s" << std::endl;

    if (save_neighbors)
    {
        std::cout << "support: save neighbor indices to file" << std::endl;
        support.SaveNeighsToFile(neighs, fname + ".neighbors" + fext);
    }

    std::cout << "DC PSE derivatives." << std::endl;
    StreamVorti::Dcpse2d derivs(gf);
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

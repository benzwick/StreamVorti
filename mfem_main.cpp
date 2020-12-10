
#include <StreamVorti/stream_vorti.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "mfem.hpp"

int main(int argc, char *argv[])
{
    // Options
    const char *mesh_file = "";
    int order = 1;
    bool static_cond = false;
    bool visualization = 1;

    // Output filename prefix and extension
    std::string fname = "mfem_square10x10";
    std::string fext = ".dat";

    // DC PSE parameters
    int CutoffRadAtNeighbor = 30;
    int SupportRadAtNeighbor = 5;

    // Output requests
    bool save_mesh = true;
    bool save_neighbors = true;
    bool save_d = true;
    bool save_dd = true;
    bool save_dx = true;
    bool save_dy = true;
    bool save_dz = true;
    bool save_dxx = true;
    bool save_dxy = true;
    bool save_dxz = true;
    bool save_dyy = true;
    bool save_dyz = true;
    bool save_dzz = true;

    // Parse command-line options
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&save_mesh,
                   "-sm", "--save-mesh", "-no-sm", "--no-save-mesh",
                   "Save mesh to file.");
    args.AddOption(&save_neighbors,
                   "-sn", "--save-neighbors", "-no-sn", "--no-save-neighbors",
                   "Save neighbors to file.");
    args.AddOption(&save_d,
                   "-sd", "--save-1st-derivative", "-no-sd", "--no-save-1st-derivative",
                   "Save 1st derivatives to file.");
    args.AddOption(&save_dd,
                   "-sdd", "--save-2nd-derivative", "-no-sdd", "--no-save-2nd-derivative",
                   "Save 2nd derivatives to file.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    // mesh_file and save_mesh are mutually exclusive options
    if (mesh_file[0] != '\0' && save_mesh)
    {
        MFEM_ABORT( "This would overwrite mesh already exists!");
    }

    // Timer
    mfem::StopWatch timer;
    timer.Start();

    // Generate (if no file provided) or read mesh from file
    mfem::Mesh *mesh;
    if (mesh_file[0] == '\0')
    {
        std::cout << "Generating a new mesh... " << std::flush;
        mesh = new mfem::Mesh(10, 10, mfem::Element::QUADRILATERAL, false, 1.0, 1.0, false);
        if (save_mesh)
        {
            std::ofstream mesh_ofs("mfem_square10x10.mesh");
            mesh_ofs.precision(8);
            mesh->Print(mesh_ofs);
        }
    }
    else
    {
        std::cout << "Read the mesh from the given mesh file... " << std::flush;
        mesh = new mfem::Mesh(mesh_file, 1, 1);
    }
    std::cout << "done." << std::endl;

    const int dim = mesh->Dimension();
    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fes(mesh, &fec, 1);

    mfem::GridFunction gf(&fes);

    std::cout << "DC PSE derivatives." << std::endl;
    timer.Clear();
    StreamVorti::Dcpse2d derivs(gf);
    std::cout << "Execution time for DCPSE derivatives initialization: "
              << timer.RealTime() << " s" << std::endl;

    if (save_neighbors)
    {
        std::cout << "support: save neighbor indices to file" << std::endl;
        derivs.SaveNeighsToFile(derivs.NeighborIndices(), fname + ".neighbors" + fext);
    }

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

    // Free the used memory
    delete mesh;

    return EXIT_SUCCESS;
}

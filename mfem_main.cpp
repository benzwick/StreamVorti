/**
 * @file mfem_main.cpp
 * @brief MFEM integration demo for DCPSE derivative operators
 *
 * Demonstrates computation of DCPSE derivative matrices on MFEM meshes.
 * Generates or loads a mesh, computes support domains and DCPSE operators,
 * and optionally saves derivative matrices and neighbor information.
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
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
 *      Benjamin F. ZWICK
 *
 * @section usage Usage
 * @code
 * MfemRun -sd -sn
 * @endcode
 */

#include <StreamVorti/mfem_main.hpp>

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
    int NumNeighbors = 35;

    // Mesh generation
    int dim = 2;
    int nx = 10;
    int ny = 10;
    int nz = 10;
    double sx = 1.0;
    double sy = 1.0;
    double sz = 1.0;

    // Output requests
    bool save_mesh = false;
    bool save_neighbors = false;
    bool save_d = false;        // all 1st derivatives (gradient)
    bool save_dd = false;       // all 2nd derivatives (Hessian)
    bool save_dx = false;
    bool save_dy = false;
    bool save_dz = false;
    bool save_dxx = false;
    bool save_dxy = false;
    bool save_dxz = false;
    bool save_dyy = false;
    bool save_dyz = false;
    bool save_dzz = false;

    // Parse command-line options
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&dim, "-dim", "--dimension",
                   "Dimension of mesh (applies to generated mesh only).");
    args.AddOption(&nx, "-nx", "--num-x-divisions",
                   "Number of x divisions.");
    args.AddOption(&ny, "-ny", "--num-y-divisions",
                   "Number of y divisions.");
    args.AddOption(&nz, "-nz", "--num-z-divisions",
                   "Number of z divisions.");
    args.AddOption(&sx, "-sx", "--size-x",
                   "Mesh size in x direction.");
    args.AddOption(&sy, "-sy", "--size-y",
                   "Mesh size in y direction.");
    args.AddOption(&sz, "-sz", "--size-z",
                   "Mesh size in z direction.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.AddOption(&NumNeighbors, "-nn", "--num-neighbors",
                   "Number of neighbors for DCPSE.");
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

    if (save_d)
    {
        save_dx = save_dy = save_dz = save_d;
    }
    if (save_dd)
    {
        save_dxx = save_dxy = save_dxz = save_dyy = save_dyz = save_dzz = save_dd;
    }

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
        std::cout << "main: Generating a new mesh... " << std::flush;
        if (dim == 2)
        {
            mesh = new mfem::Mesh(nx, ny, mfem::Element::QUADRILATERAL,
                                  false, sx, sy, false);
        }
        else if (dim == 3)
        {
            mesh = new mfem::Mesh(nx, ny, nz, mfem::Element::HEXAHEDRON,
                                  false, sx, sy, sz, false);
        }
        else
        {
            MFEM_ABORT( "Unsupported mesh dimension: " << dim );
        }
        if (save_mesh)
        {
            std::ofstream mesh_ofs("mfem_square10x10.mesh");
            mesh_ofs.precision(8);
            mesh->Print(mesh_ofs);
        }
    }
    else
    {
        std::cout << "main: Read the mesh from the given mesh file... " << std::flush;
        mesh = new mfem::Mesh(mesh_file, 1, 1);
    }
    std::cout << "done." << std::endl;

    dim = mesh->Dimension();    // The actual mesh dimension
    std::cout << "main: Mesh dimension: " << dim << std::endl;
    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fes(mesh, &fec, 1);

    mfem::GridFunction gf(&fes);

    std::cout << "main: DC PSE derivatives." << std::endl;
    timer.Clear();
    StreamVorti::Dcpse *derivs;
    if (dim == 2)
    {
        derivs = new StreamVorti::Dcpse2d(
            gf, NumNeighbors);
    }
    else if (dim == 3)
    {
        derivs = new StreamVorti::Dcpse3d(
            gf, NumNeighbors);
    }
    else
    {
        MFEM_ABORT( "Unsupported dimension: " << dim << "." );
    }
    std::cout << "main: Execution time for DCPSE derivatives initialization: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    derivs->Update();
    std::cout << "main: Execution time for DCPSE derivatives calculation: "
              << timer.RealTime() << " s" << std::endl;

    if (save_neighbors)
    {
        std::cout << "main: Save neighbor indices to file... " << std::endl;
        derivs->SaveNeighsToFile(derivs->NeighborIndices(), fname + ".neighbors" + fext);
        std::cout << "done." << std::endl;
    }

    std::cout << "main: Save derivative operator matrices to file... " << std::flush;
    if (dim > 1) {if (save_dx)  {derivs->SaveDerivToFile("dx",  fname + ".dx"  + fext);}}
    if (dim > 1) {if (save_dy)  {derivs->SaveDerivToFile("dy",  fname + ".dy"  + fext);}}
    if (dim > 2) {if (save_dz)  {derivs->SaveDerivToFile("dz",  fname + ".dz"  + fext);}}
    if (dim > 1) {if (save_dxx) {derivs->SaveDerivToFile("dxx", fname + ".dxx" + fext);}}
    if (dim > 1) {if (save_dxy) {derivs->SaveDerivToFile("dxy", fname + ".dxy" + fext);}}
    if (dim > 2) {if (save_dxz) {derivs->SaveDerivToFile("dxz", fname + ".dxz" + fext);}}
    if (dim > 1) {if (save_dyy) {derivs->SaveDerivToFile("dyy", fname + ".dyy" + fext);}}
    if (dim > 2) {if (save_dyz) {derivs->SaveDerivToFile("dyz", fname + ".dyz" + fext);}}
    if (dim > 2) {if (save_dzz) {derivs->SaveDerivToFile("dzz", fname + ".dzz" + fext);}}
    std::cout << "done." << std::endl;

    // Free the used memory
    delete mesh;

    std::cout << "main: success!" << std::endl;

    return EXIT_SUCCESS;
}

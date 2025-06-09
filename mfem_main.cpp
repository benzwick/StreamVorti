/*
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
 */

// Demo usage:
//     ./MfemRun -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -nn 25 -sm -sn -sd -sdd

#include <StreamVorti/stream_vorti.hpp>
#include "mfem.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


// modular main structure
// 1. init parameters: mesh, FES, DCPSE derives
// 2. CreateOrLoadMesh()
// 3. Set up FES
// 3. Initialise DCPSE
// 4. Save derivative matrices
// 5. save neighbours



// Stream-Vorti simulation (vorticity, Gersc time step, boundary condition)
// 1. get DCPSE derivative matrices
// Define boundaries (LDC)
// Initialize vorticity and stream function
// 3. update vorticty (eq. 11)
// 4. apply boundary condition
// 5. compute Gershgorin time step (2.6)
// 5. solve streamfunction (eq. 12)
// 6. update velocity (eq. 13)
// run simulation (main time steps for loop)
// save solution


// Simulation parameters structure
struct SimulationParams {
    double final_time = 60.0;
    double dt = 1e-3;
    double reynolds_number = 1000.0;

    std::string output_prefix = "mfem_square10x10";
    std::string output_extension = ".dat";

};

// Function declarations
// Refactored functions
mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh);
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors);
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                            int dim, bool save_d, bool save_dd);

// New functions to match "Explicit_streamfunction_vorticity_meshless.m"
void IdentifyBoundaryNodesLDC(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                            std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                            std::vector<int>& left_nodes, std::vector<int>& interior_nodes);


int main(int argc, char *argv[])
{
    // Options
    const char *mesh_file = "";
    int order = 1;
    bool visualization = 1;

    // Output filename prefix and extension
    std::string fname = "mfem_square10x10";
    std::string fext = ".dat";

    // DC PSE parameters
    int NumNeighbors = 25;

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


    // Simulation parameters
    SimulationParams params;

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

    // Create or load mesh
    mfem::Mesh* mesh = CreateOrLoadMesh(mesh_file, dim, nx, ny, nz, sx, sy, sz, save_mesh);

    // Save mesh if requested
    if (save_mesh) {
        std::ofstream mesh_ofs(fname+fext);
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
    }

    // Set up finite element space
    dim = mesh->Dimension();
    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fes(mesh, &fec, 1);
    mfem::GridFunction gf(&fes);

    // Initialise DCPSE derivatives
    StreamVorti::Dcpse* derivs = InitialiseDCPSE(gf, dim, NumNeighbors);
    std::cout << "main: DC PSE derivatives." << std::endl;

    // save derivs matrices
    SaveDerivativeMatrices(derivs, params, dim, save_d, save_dd);

    // Save neighbors if requested
    if (save_neighbors) {
        std::cout << "main: Save neighbor indices to file... " << std::endl;
        derivs->SaveNeighsToFile(derivs->NeighborIndices(), params.output_prefix + ".neighbors" + params.output_extension);
        std::cout << "done." << std::endl;
    }

    /*********************** Simulation *************************/

    // Ensure we have a 2D DCPSE object first
    StreamVorti::Dcpse2d* dcpse2d = dynamic_cast<StreamVorti::Dcpse2d*>(derivs);
    if (!dcpse2d) {
        MFEM_ABORT("RunSimulation: Only 2D simulations are currently supported.");
    }

    // Get derivative matrices dx,dy,dxx,dyy
    const mfem::SparseMatrix& dx_matrix = dcpse2d->ShapeFunctionDx();
    const mfem::SparseMatrix& dy_matrix = dcpse2d->ShapeFunctionDy();
    const mfem::SparseMatrix& dxx_matrix = dcpse2d->ShapeFunctionDxx();
    const mfem::SparseMatrix& dyy_matrix = dcpse2d->ShapeFunctionDxy();

    std::cout << "RunSimulation: Retrieved DCPSE derivative matrices successfully." << std::endl;

    // Identify boundary and interior nodes
    std::vector<int> bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes;
    IdentifyBoundaryNodesLDC(mesh, bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes);

    // Combine all boundary nodes for streamfunction boundary conditions
    std::vector<int> all_boundary_nodes;
    all_boundary_nodes.insert(all_boundary_nodes.end(), bottom_nodes.begin(), bottom_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), right_nodes.begin(), right_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), top_nodes.begin(), top_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), left_nodes.begin(), left_nodes.end());

    /****************** Streamfunction *******************/
    // Initialise solution fields
    mfem::Vector vorticity, streamfunction;
    const int num_nodes = mesh->GetNV();
    InitialiseFields(num_nodes, vorticity, streamfunction);

    /* TODO: GLVis
    // Send the above data by socket to a GLVis server.  Use the "n" and "b"
    // keys in GLVis to visualize the displacements.
    if (visualization) {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << flush;
    }
    */


    // Free the used memory
    delete mesh;

    std::cout << "main: success!" << std::endl;

    return EXIT_SUCCESS;
}

/**
 * @brief Create or load mesh object
 *
 * @param mesh_file
 * @param dim
 * @param nx
 * @param ny
 * @param nz
 * @param sx
 * @param sy
 * @param sz
 * @param save_mesh
 * @return mfem::Mesh*
 */
mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh)
{
    mfem::Mesh* mesh;

    if (mesh_file[0] == '\0') {
        std::cout << "CreateOrLoadMesh: Generating a new mesh... " << std::flush;

        if (dim == 2) {
            mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, false, sx, sy, false));
        } else if (dim == 3) {
            mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian3D(nx, ny, nz, mfem::Element::HEXAHEDRON, false, sx, sy, sz));
        } else {
                MFEM_ABORT("Unsupported mesh dimension: " << dim);
        }

    } else {
        std::cout << "CreateOrLoadMesh: Read the mesh from the given mesh file... " << std::flush;
        mesh = new mfem::Mesh(mesh_file, 1, 1);
    }
    std::cout << "done." << std::endl;

    return mesh;
}

/**
 * @brief Initialise DC PSE
 *
 * @param gf
 * @param dim
 * @param NumNeighbors
 * @return StreamVorti::Dcpse*
 */
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors)
{
    std::cout << "InitialiseDCPSE: Initialising DC PSE derivatives." << std::endl;
    mfem::StopWatch timer;
    timer.Start();

    StreamVorti::Dcpse* derivs;
    if (dim == 2) {
        derivs = new StreamVorti::Dcpse2d(gf, NumNeighbors);
    } else if (dim == 3) {
        derivs = new StreamVorti::Dcpse3d(gf, NumNeighbors);
    } else {
        MFEM_ABORT("Unsupported dimension: " << dim << ".");
    }

    std::cout << "InitialiseDCPSE: Execution time for DCPSE derivatives initialisation: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    derivs->Update();
    std::cout << "InitialiseDCPSE: Execution time for DCPSE derivatives calculation: "
              << timer.RealTime() << " s" << std::endl;

    return derivs;
}

/**
 * @brief Save derivative matrices
 *
 * @param derivs
 * @param params
 * @param dim
 * @param save_d
 * @param save_dd
 */
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                           int dim, bool save_d, bool save_dd)
{
    if (!save_d && !save_dd) return;

    std::cout << "SaveDerivativeMatrices: Save derivative operator matrices to file... " << std::flush;

    if (save_d) {
        if (dim > 1) derivs->SaveDerivToFile("dx", params.output_prefix + ".dx" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dy", params.output_prefix + ".dy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dz", params.output_prefix + ".dz" + params.output_extension);
    }

    if (save_dd) {
        if (dim > 1) derivs->SaveDerivToFile("dxx", params.output_prefix + ".dxx" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dxy", params.output_prefix + ".dxy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dxz", params.output_prefix + ".dxz" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dyy", params.output_prefix + ".dyy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dyz", params.output_prefix + ".dyz" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dzz", params.output_prefix + ".dzz" + params.output_extension);
    }

    std::cout << "done." << std::endl;
}


/* Functions to match "Explicit_streamfunction_vorticity_meshless.m" */
//TODO: try define boundaries for other benchmark problems: Backward-Facing Step, Unbounded Flow Past a Cylinder
/**
 * @brief Identify boundary ndoes for Lid-Driven Cavity problem
 *
 * @param mesh
 * @param bottom_nodes
 * @param right_nodes
 * @param top_nodes
 * @param left_nodes
 * @param interior_nodes
 */
void IdentifyBoundaryNodesLDC(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                          std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                          std::vector<int>& left_nodes, std::vector<int>& interior_nodes)
{
    const double tolerance = 1e-10;
    const int num_vertices = mesh->GetNV();

    // Clear all vectors
    bottom_nodes.clear();
    right_nodes.clear();
    top_nodes.clear();
    left_nodes.clear();
    interior_nodes.clear();

    for (int i = 0; i < num_vertices; ++i) {
        const double* vertex = mesh->GetVertex(i);
        double x = vertex[0];
        double y = vertex[1];

        if (std::abs(y) < tolerance) {
            bottom_nodes.push_back(i);  // y = 0
        } else if (std::abs(x - 1.0) < tolerance && y > tolerance && y < 1.0 - tolerance) {
            right_nodes.push_back(i);   // x = 1, 0 < y < 1
        } else if (std::abs(y - 1.0) < tolerance) {
            top_nodes.push_back(i);     // y = 1
        } else if (std::abs(x) < tolerance && y > tolerance && y < 1.0 - tolerance) {
            left_nodes.push_back(i);    // x = 0, 0 < y < 1
        } else {
            interior_nodes.push_back(i);
        }
    }

    std::cout << "IdentifyBoundaryNodesLDC: Found " << bottom_nodes.size() << " bottom, "
              << right_nodes.size() << " right, " << top_nodes.size() << " top, "
              << left_nodes.size() << " left, " << interior_nodes.size() << " interior nodes."
              << std::endl;
}


void InitialiseFields(int num_nodes, mfem::Vector& vorticity, mfem::Vector& streamfunction)
{
    vorticity.SetSize(num_nodes);
    streamfunction.SetSize(num_nodes);
    vorticity = 0.0;
    streamfunction = 0.0;

    std::cout << "InitialiseFields: Initialised vorticity and streamfunction fields with "
              << num_nodes << " nodes." << std::endl;
}
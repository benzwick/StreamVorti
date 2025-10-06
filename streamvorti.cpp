/**
 * @file streamvorti.cpp
 * @brief StreamVorti lid-driven cavity flow simulation
 *
 * Implements stream function-vorticity formulation for incompressible
 * Navier-Stokes equations using DCPSE meshless method on MFEM meshes.
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
 * ./StreamVorti -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -nn 25 -sm -sn -sd -sdd -ss -pv
 * @endcode
 */

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
// 1. get DCPSE derivative matrices - DONE
// 2. Define boundaries (use MFEM method) - DONE
// Initialize vorticity and stream function -DONE
// Linear solver (Mfem example)
// 3. update vorticty (eq. 11)
// apply boundary condition
// 4. compute Gershgorin time step (2.6)
// 5. solve streamfunction (eq. 12)
// 6. update velocity (eq. 13)
// run simulation (main time steps for loop)
// save solution


/**
 * @brief Simulation parameters structure
 *
 * Container for all simulation parameters including time stepping,
 * Reynolds number, and output settings.
 */
struct SimulationParams {
    double final_time = 60.0;                      ///< Final simulation time
    double dt = 1e-3;                              ///< Initial time step size
    double reynolds_number = 1000.0;               ///< Reynolds number for flow simulation

    std::string output_prefix = "mfem_square10x10"; ///< Prefix for output filenames
    std::string output_extension = ".dat";          ///< Extension for output files

    int output_frequency = 100;                     ///< Output solution every N timesteps
    int gershgorin_frequency = 1000;                ///< Recompute Gershgorin timestep every N steps
    bool enable_adaptive_timestep = true;           ///< Enable adaptive timestepping
    double gershgorin_safety_factor = 5.0;          ///< Safety factor for Gershgorin timestep

};

/**
 * @brief Create or load a mesh
 * @param mesh_file Path to mesh file (empty string to generate mesh)
 * @param dim Mesh dimension (2 or 3)
 * @param nx Number of divisions in x direction
 * @param ny Number of divisions in y direction
 * @param nz Number of divisions in z direction
 * @param sx Domain size in x direction
 * @param sy Domain size in y direction
 * @param sz Domain size in z direction
 * @param save_mesh Whether to save the generated mesh
 * @return Pointer to created or loaded mesh
 */
mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh);

/**
 * @brief Initialize DCPSE derivative operators
 * @param gf Grid function containing mesh and DOF information
 * @param dim Spatial dimension (2 or 3)
 * @param NumNeighbors Number of neighbors for DCPSE stencil
 * @return Pointer to initialized DCPSE object (Dcpse2d or Dcpse3d)
 */
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors);

/**
 * @brief Save DCPSE derivative matrices to files
 * @param derivs Pointer to DCPSE object containing derivative operators
 * @param params Simulation parameters for output naming
 * @param dim Spatial dimension
 * @param save_d Save first derivative operators (dx, dy, dz)
 * @param save_dd Save second derivative operators (dxx, dxy, etc.)
 * @param dat_dir Directory for output files
 */
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                            int dim, bool save_d, bool save_dd, std::string dat_dir);

/**
 * @brief Identify boundary nodes for lid-driven cavity problem
 * @param mesh Pointer to the mesh
 * @param bottom_nodes Output vector of bottom boundary node indices (y=0)
 * @param right_nodes Output vector of right boundary node indices (x=1)
 * @param top_nodes Output vector of top boundary node indices (y=1)
 * @param left_nodes Output vector of left boundary node indices (x=0)
 * @param interior_nodes Output vector of interior node indices
 */
void IdentifyBoundaryNodesLDC(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                            std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                            std::vector<int>& left_nodes, std::vector<int>& interior_nodes);

/**
 * @brief Create Laplacian matrix from second derivative operators
 * @param dxx_matrix Second derivative operator in x direction
 * @param dyy_matrix Second derivative operator in y direction
 * @return Sparse matrix representing Laplacian operator (∇²)
 */
mfem::SparseMatrix CreateLaplacianMatrix(const mfem::SparseMatrix dxx_matrix, const mfem::SparseMatrix dyy_matrix);

/**
 * @brief Apply Dirichlet boundary conditions to a sparse matrix
 * @param matrix Matrix to modify (in-place)
 * @param boundary_nodes Indices of boundary nodes
 */
static void ApplyDirichletBC(mfem::SparseMatrix& matrix, const std::vector<int>& boundary_nodes);

/**
 * @brief Create Laplacian matrix with Dirichlet boundary conditions applied
 * @param boundary_nodes Indices of boundary nodes
 * @param dxx_matrix Second derivative operator in x direction
 * @param dyy_matrix Second derivative operator in y direction
 * @return Laplacian matrix with boundary conditions
 */
mfem::SparseMatrix CreateLaplacianWithBC(const std::vector<int>& boundary_nodes,
                                        const mfem::SparseMatrix dxx_matrix,
                                        const mfem::SparseMatrix dyy_matrix);

/**
 * @brief Update vorticity field using explicit Euler time stepping
 * @param dx_matrix First derivative operator in x direction
 * @param dy_matrix First derivative operator in y direction
 * @param dxx_matrix Second derivative operator in x direction
 * @param dyy_matrix Second derivative operator in y direction
 * @param streamfunction Current streamfunction field
 * @param vorticity Vorticity field (updated in-place)
 * @param dt Time step size
 * @param reynolds_number Reynolds number
 */
void UpdateVorticity(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                    const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                    const mfem::Vector& streamfunction, mfem::Vector& vorticity,
                    double dt, double reynolds_number);

/**
 * @brief Apply boundary conditions for vorticity in lid-driven cavity
 * @param bottom_nodes Bottom boundary node indices
 * @param right_nodes Right boundary node indices
 * @param top_nodes Top boundary node indices (moving lid)
 * @param left_nodes Left boundary node indices
 * @param dx_matrix First derivative operator in x direction
 * @param dy_matrix First derivative operator in y direction
 * @param streamfunction Current streamfunction field
 * @param vorticity Vorticity field (updated in-place)
 */
void ApplyBoundaryConditions(const std::vector<int>& bottom_nodes, const std::vector<int>& right_nodes,
                           const std::vector<int>& top_nodes, const std::vector<int>& left_nodes,
                           const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                           const mfem::Vector& streamfunction, mfem::Vector& vorticity);

/**
 * @brief Compute critical time step using Gershgorin circle theorem
 * @param dx_matrix First derivative operator in x direction
 * @param dy_matrix First derivative operator in y direction
 * @param dxx_matrix Second derivative operator in x direction
 * @param dyy_matrix Second derivative operator in y direction
 * @param streamfunction Current streamfunction field
 * @param reynolds_number Reynolds number
 * @return Critical time step for stability
 */
double ComputeGershgorinTimeStep(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                                const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                                const mfem::Vector& streamfunction, double reynolds_number);

/**
 * @brief Save vorticity and streamfunction fields to files
 * @param vorticity Vorticity field vector
 * @param streamfunction Streamfunction field vector
 * @param filename Base filename for output
 * @param timestep Current timestep number
 * @param dat_dir Directory for output files
 */
void SaveSolutionToFile(const mfem::Vector& vorticity, const mfem::Vector& streamfunction, const std::string& filename, int timestep, std::string dat_dir);

// visulisation
#include <iomanip>
#include <sstream>


/**
 * @brief Main function for StreamVorti lid-driven cavity simulation
 *
 * Implements the stream function-vorticity formulation for solving the
 * incompressible Navier-Stokes equations using DCPSE (DC Particle Strength Exchange)
 * meshless method on MFEM meshes.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return EXIT_SUCCESS on successful completion
 */
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
    int NumNeighbors = 10;

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
    bool save_solutions = false;

    int vtu_frequency = 1000;
    bool save_vtu = false;

    bool paraview_output = false;
    std::string paraview_filename = fname;


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
    args.AddOption(&save_solutions,
                   "-ss", "--save-solutions", "-no-ss", "--no-save-solutions",
                   "Save simuation solutions to file.");
    args.AddOption(&vtu_frequency, "-vf", "--vtu-frequency",
                   "Frequency for VTU output (every N timesteps).");
    args.AddOption(&save_vtu, "-vtu", "--save-vtu", "-no-vtu", "--no-save-vtu",
                   "Enable VTU output for ParaView.");
    args.AddOption(&paraview_output, "-pv", "--paraview", "-no-pv",
                  "--no-paraview",
                  "Enable or disable ParaView output.");
    args.Parse();

    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    // Create or load mesh
    mfem::Mesh* mesh = CreateOrLoadMesh(mesh_file, dim, nx, ny, nz, sx, sy, sz, save_mesh);

    // dat files for Matlab
    std::string dat_dir;
    // dat_dir = params.output_prefix + "_dat";
    dat_dir = "output_dat";
    std::filesystem::create_directories(dat_dir);
    std::cout << "Created DAT output directory: " << dat_dir << std::endl;

    // Save mesh if requested
    if (save_mesh) {
        // Save mesh as MFEM format
        std::ofstream mesh_ofs(fname+".mesh");
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
        // Save mesh as VTU format for ParaView
        //mesh->PrintVTU(fname+".elem");
        //mesh->PrintBdrVTU(fname+".bdr");
    }

    // Set up finite element space
    dim = mesh->Dimension();
    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fes(mesh, &fec, 1);
    std::cout << "main: Number of finite element unknowns: " << fes.GetTrueVSize() << std::endl;
    mfem::GridFunction gf(&fes); // as 'x' in ex1.cpp

    // Initialise DCPSE derivatives
    StreamVorti::Dcpse* derivs = InitialiseDCPSE(gf, dim, NumNeighbors);
    std::cout << "main: DC PSE derivatives." << std::endl;

    // save derivs matrices
    SaveDerivativeMatrices(derivs, params, dim, save_d, save_dd, dat_dir);

    // Save neighbors if requested
    if (save_neighbors) {
        std::cout << "main: Save neighbor indices to file... " << std::endl;
        derivs->SaveNeighsToFile(derivs->NeighborIndices(), dat_dir + "/" + params.output_prefix + ".neighbors" + params.output_extension);
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
    const mfem::SparseMatrix& dyy_matrix = dcpse2d->ShapeFunctionDyy();
    std::cout << "RunSimulation: Retrieved DCPSE derivative matrices successfully." << std::endl;


    // Identify boundary and interior nodes
    // TODO: use Mfem attributes
    std::vector<int> bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes;
    IdentifyBoundaryNodesLDC(mesh, bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes);


    /****************** Streamfunction *******************/
    // Initialise solution fields
    mfem::Vector vorticity, streamfunction, STREAMFUNCTION;
    const int num_nodes = mesh->GetNV();
    int num_timesteps = static_cast<int>(params.final_time / params.dt);
    double current_dt = params.dt;

    vorticity.SetSize(num_nodes);
    streamfunction.SetSize(num_nodes);
    vorticity = 0.0;
    streamfunction = 0.0;

    std::cout << "InitialiseFields: Initialised vorticity and streamfunction fields with "
              << num_nodes << " nodes." << std::endl;

    // Combine all boundary nodes for streamfunction boundary conditions
    std::vector<int> all_boundary_nodes;
    all_boundary_nodes.insert(all_boundary_nodes.end(), bottom_nodes.begin(), bottom_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), right_nodes.begin(), right_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), top_nodes.begin(), top_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), left_nodes.begin(), left_nodes.end());

    // Create Laplacian matrix with boundary conditions for streamfunction equation
    mfem::SparseMatrix laplacian_matrix = CreateLaplacianWithBC(all_boundary_nodes, dxx_matrix, dyy_matrix);
    std::cout << "RunSimulation: Created Laplacian matrix with boundary conditions." << std::endl;
    //std::cout << "Size of linear system 'Laplacian matrix': " << A->Height() << std::endl;


    std::cout << "RunSimulation: Running " << num_timesteps << " time steps..." << std::endl;
    mfem::StopWatch sim_timer;
    sim_timer.Start();

     // Paraview output file

    // Calculate velocity from streamfunction
    mfem::Vector u_velocity(streamfunction.Size());
    mfem::Vector v_velocity(streamfunction.Size());

    dy_matrix.Mult(streamfunction, u_velocity);   // u = ∂ψ/∂y
    dx_matrix.Mult(streamfunction, v_velocity);   // v = -∂ψ/∂x
    v_velocity *= -1.0;

    // Create finite element spaces
    //mfem::H1_FECollection fec(1, mesh->Dimension());
    mfem::FiniteElementSpace scalar_fes(const_cast<mfem::Mesh*>(mesh), &fec, 1);
    mfem::FiniteElementSpace vector_fes(const_cast<mfem::Mesh*>(mesh), &fec, mesh->Dimension());

    // Create grid functions
    mfem::GridFunction vorticity_gf(&scalar_fes);
    mfem::GridFunction streamfunction_gf(&scalar_fes);
    mfem::GridFunction velocity_gf(&vector_fes);

    // Create ParaView data collection
    // mfem::ParaViewDataCollection paraview_dc(filename_base, const_cast<mfem::Mesh*>(mesh));
    mfem::ParaViewDataCollection paraview_dc(paraview_filename, mesh);

    if (paraview_output) {
        // paraview_dc.SetPrecision(8);
        paraview_dc.SetPrefixPath("ParaView");
        paraview_dc.SetLevelsOfDetail(order);
        // paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
        paraview_dc.SetDataFormat(mfem::VTKFormat::BINARY);
        // paraview_dc.SetHighOrderOutput(true);
        // Register fields
        paraview_dc.RegisterField("Vorticity", &vorticity_gf);
        paraview_dc.RegisterField("StreamFunction", &streamfunction_gf);
        paraview_dc.RegisterField("Velocity", &velocity_gf);
        // Set time and cycle, then save
        paraview_dc.SetCycle(0);
        paraview_dc.SetTime(0.0);
        paraview_dc.Save();

        std::cout << "Created ParaView output directory: " << "ParaView" << std::endl;
    }

#ifndef MFEM_USE_SUITESPARSE
    std::cout << "Linear Solver: CGSolver with GSSmoother" << std::endl;
    // Iterative solver
    mfem::CGSolver solver;
    mfem::GSSmoother preconditioner;
    // Setup preconditioner

    // Setup CG solver
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(1000);
    solver.SetPrintLevel(0);  // Set to 1 for debugging
    solver.SetPreconditioner(preconditioner);
    solver.SetOperator(laplacian_matrix);
#else
    std::cout << "Linear Solver: UMFPackSolver" << std::endl;
    std::cout << "Using UMFPackSolver linear solver" << std::endl;
    // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    // Direct solver
    mfem::UMFPackSolver solver;
    solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    solver.SetOperator(laplacian_matrix);
#endif

    mfem::Vector rhs;

    // Main simulation loop (implementing Algorithm from Bourantas et al. 2019)
    for (int time_step = 1; time_step <= num_timesteps; ++time_step) {
        double current_time = time_step * current_dt;

        if (time_step % params.output_frequency == 0) {
            // std::cout << "Time step: " << time_step << " / " << num_timesteps
            //           << ", t = " << time_step * current_dt
            //           << ", dt = " << current_dt << std::endl;
        }
        /****************** Vorticity *******************/
        // Step 1: Update vorticity using explicit Euler scheme (Equation 11)

        UpdateVorticity(dx_matrix, dy_matrix, dxx_matrix, dyy_matrix, streamfunction, vorticity, current_dt, params.reynolds_number);

        // Step 2: Apply vorticity boundary conditions (Equation 33)

        ApplyBoundaryConditions(bottom_nodes, right_nodes, top_nodes, left_nodes,dx_matrix, dy_matrix, streamfunction, vorticity);

        // Step 3: Solve Poisson equation for streamfunction (Equation 12)

        // SolveStreamFunction(laplacian_matrix,vorticity, all_boundary_nodes, streamfunction);
	// Set up the Poisson equation: ∇²ψ = -ω (Equation 12 from paper)
	rhs = vorticity;
	rhs *= -1.0;

	// Apply homogeneous Dirichlet boundary conditions: ψ = 0 on all boundaries
	for (int idx : all_boundary_nodes) {
	  rhs[idx] = 0.0;
	}

	solver.Mult(rhs, streamfunction);

        // Step 4: Adaptive time stepping using Gershgorin circle theorem
        if (params.enable_adaptive_timestep && time_step % params.gershgorin_frequency == 0) {
            double dt_critical = ComputeGershgorinTimeStep(dx_matrix, dy_matrix, dxx_matrix, dyy_matrix,
                                                          streamfunction, params.reynolds_number);

            if (current_dt > dt_critical / params.gershgorin_safety_factor) {
                current_dt = dt_critical / params.gershgorin_safety_factor;
                num_timesteps = static_cast<int>(params.final_time / current_dt);
                std::cout << "Adjusted time step to dt = " << current_dt
                          << ", new total steps = " << num_timesteps << std::endl;
            }
        }

        // Step 5: Output solution periodically
        if (save_solutions && (time_step % (params.output_frequency * 10) == 0)) {
            SaveSolutionToFile(vorticity, streamfunction,
                             params.output_prefix + "_solution", time_step, dat_dir);
        }

        // Step 6 Save Paraview output
        if (paraview_output && (time_step % 100 == 0)) {

            // Copy data to grid functions
            for (int i = 0; i < vorticity.Size(); ++i) {
                vorticity_gf[i] = vorticity[i];
                streamfunction_gf[i] = streamfunction[i];
            }

            // Set velocity components
            for (int i = 0; i < u_velocity.Size(); ++i) {
                velocity_gf[i] = u_velocity[i];                    // x-component
                if (mesh->Dimension() > 1) {
                    velocity_gf[i + u_velocity.Size()] = v_velocity[i]; // y-component
                }
            }

            // save paraview data
            paraview_dc.SetTime(current_time);
            paraview_dc.SetCycle(time_step);
            paraview_dc.Save();
        }

        // Step 7: Check for convergence (steady state)
        if (time_step > 1000 && time_step % 1000 == 0) {
            // Compute L2 norm of vorticity change
            static mfem::Vector prev_vorticity;
            if (prev_vorticity.Size() == vorticity.Size()) {
                mfem::Vector diff = vorticity;
                diff -= prev_vorticity;
                double change_norm = diff.Norml2() / vorticity.Norml2();

                if (change_norm < 1e-8) {
                    std::cout << "Converged to steady state at timestep " << time_step
                              << " (relative change = " << change_norm << ")" << std::endl;
                    break;
                }
            }
            prev_vorticity = vorticity;
        }
    }

    std::cout << "RunSimulation: Simulation completed successfully in "
              << sim_timer.RealTime() << " seconds." << std::endl;

    // Save final solution
    if (save_solutions) {
            SaveSolutionToFile(vorticity, streamfunction, params.output_prefix + "_final", num_timesteps, dat_dir);
    }


    // Output simulation statistics
    std::cout << "Final simulation statistics:" << std::endl;
    std::cout << "  - Final time step size: " << current_dt << std::endl;
    std::cout << "  - Total time steps completed: " << std::min(num_timesteps, static_cast<int>(params.final_time / current_dt)) << std::endl;
    std::cout << "  - Maximum vorticity: " << vorticity.Max() << std::endl;
    std::cout << "  - Minimum vorticity: " << vorticity.Min() << std::endl;
    std::cout << "  - Maximum streamfunction: " << streamfunction.Max() << std::endl;
    std::cout << "  - Minimum streamfunction: " << streamfunction.Min() << std::endl;

    /*
    TODO: visualization on GLVis or Paraview
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
    delete derivs;
    delete mesh;

    std::cout << "main: success!" << std::endl;
    return EXIT_SUCCESS;
}


/**
 * @brief Create a new mesh or load from file
 *
 * Generates a Cartesian mesh (2D quadrilateral or 3D hexahedral) if no mesh file
 * is provided, otherwise loads the mesh from the specified file.
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
 * @brief Initialize DCPSE derivative operators from mesh
 *
 * Creates a Dcpse2d or Dcpse3d object based on the spatial dimension and
 * computes all derivative operators (first and second derivatives).
 * Reports timing for initialization and computation phases.
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
 * @brief Save derivative operator matrices to files
 *
 * Exports the DCPSE derivative operators to data files for analysis or
 * verification. Saves first derivatives (dx, dy, dz) and/or second
 * derivatives (dxx, dxy, dxz, dyy, dyz, dzz) based on flags.
 */
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                           int dim, bool save_d, bool save_dd, std::string dat_dir)
{
    if (!save_d && !save_dd) return;

    std::cout << "SaveDerivativeMatrices: Save derivative operator matrices to file... " << std::flush;

    if (save_d) {
        if (dim > 1) derivs->SaveDerivToFile("dx", dat_dir + "/" + params.output_prefix + ".dx" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dy", dat_dir + "/" + params.output_prefix + ".dy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dz", dat_dir + "/" + params.output_prefix + ".dz" + params.output_extension);
    }

    if (save_dd) {
        if (dim > 1) derivs->SaveDerivToFile("dxx", dat_dir + "/" + params.output_prefix + ".dxx" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dxy", dat_dir + "/" + params.output_prefix + ".dxy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dxz", dat_dir + "/" + params.output_prefix + ".dxz" + params.output_extension);
        if (dim > 1) derivs->SaveDerivToFile("dyy", dat_dir + "/" + params.output_prefix + ".dyy" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dyz", dat_dir + "/" + params.output_prefix + ".dyz" + params.output_extension);
        if (dim > 2) derivs->SaveDerivToFile("dzz", dat_dir + "/" + params.output_prefix + ".dzz" + params.output_extension);
    }

    std::cout << "done." << std::endl;
}


/**
 * @brief Identify boundary nodes for lid-driven cavity problem
 *
 * Classifies mesh vertices into boundary segments (bottom, right, top, left)
 * and interior nodes based on their coordinates. Assumes a unit square domain [0,1]×[0,1].
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



/**
 * @brief Create Laplacian matrix from second derivative operators
 *
 * Constructs the discrete Laplacian operator ∇² = ∂²/∂x² + ∂²/∂y² by summing
 * the dxx and dyy matrices element-wise.
 */
mfem::SparseMatrix CreateLaplacianMatrix(const mfem::SparseMatrix dxx_matrix, const mfem::SparseMatrix dyy_matrix) {
    const mfem::SparseMatrix& dxx = dxx_matrix;
    const mfem::SparseMatrix& dyy = dyy_matrix;

    // Proper MFEM matrix addition
    mfem::SparseMatrix laplacian(dxx.Height(), dxx.Width());

    // Manual addition since MFEM Add() might not work as expected
    for (int i = 0; i < dxx.Height(); ++i) {
        for (int j = dxx.GetI()[i]; j < dxx.GetI()[i + 1]; ++j) {
            int col = dxx.GetJ()[j];
            double val = dxx.GetData()[j];
            laplacian.Add(i, col, val);
        }
    }

    for (int i = 0; i < dyy.Height(); ++i) {
        for (int j = dyy.GetI()[i]; j < dyy.GetI()[i + 1]; ++j) {
            int col = dyy.GetJ()[j];
            double val = dyy.GetData()[j];
            laplacian.Add(i, col, val);
        }
    }

    laplacian.Finalize();
    return laplacian;
}


/**
 * @brief Apply Dirichlet boundary conditions to a sparse matrix
 *
 * Modifies the matrix in-place by zeroing out rows corresponding to boundary nodes
 * and setting diagonal entries to 1. This enforces homogeneous Dirichlet BCs.
 */
static void ApplyDirichletBC(mfem::SparseMatrix& matrix, const std::vector<int>& boundary_nodes) {
    for (int boundary_idx : boundary_nodes) {
        // Zero out the row
        for (int j = matrix.GetI()[boundary_idx]; j < matrix.GetI()[boundary_idx + 1]; ++j) {
            matrix.GetData()[j] = 0.0;
        }

        // Set diagonal entry to 1
        for (int j = matrix.GetI()[boundary_idx]; j < matrix.GetI()[boundary_idx + 1]; ++j) {
            if (matrix.GetJ()[j] == boundary_idx) {
                matrix.GetData()[j] = 1.0;
                break;
            }
        }
    }
}

/**
 * @brief Create Laplacian matrix with boundary conditions applied
 *
 * Combines CreateLaplacianMatrix and ApplyDirichletBC into a single operation,
 * returning a Laplacian operator ready for solving the Poisson equation.
 */
mfem::SparseMatrix CreateLaplacianWithBC(const std::vector<int>& boundary_nodes,
                                        const mfem::SparseMatrix dxx_matrix,
                                        const mfem::SparseMatrix dyy_matrix) {
    mfem::SparseMatrix laplacian = CreateLaplacianMatrix(dxx_matrix, dyy_matrix);
    ApplyDirichletBC(laplacian, boundary_nodes);
    return laplacian;
}

/**
 * @brief Update vorticity field using explicit Euler time stepping
 *
 * Implements Equation 11 from Bourantas et al. (2019):
 * ω^(n+1) = ω^n + dt*[(1/Re)*∇²ω - (∂ψ/∂y)(∂ω/∂x) + (∂ψ/∂x)(∂ω/∂y)]
 * where the second and third terms represent convection.
 */
void UpdateVorticity(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                    const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                    const mfem::Vector& streamfunction, mfem::Vector& vorticity,
                    double dt, double reynolds_number)
{
    // Compute spatial derivatives
    mfem::Vector dpsi_dy(streamfunction.Size());
    mfem::Vector dpsi_dx(streamfunction.Size());
    mfem::Vector domega_dx(vorticity.Size());
    mfem::Vector domega_dy(vorticity.Size());
    mfem::Vector d2omega_dx2(vorticity.Size());
    mfem::Vector d2omega_dy2(vorticity.Size());

    // Compute derivatives
    dy_matrix.Mult(streamfunction, dpsi_dy);      // ∂ψ/∂y
    dx_matrix.Mult(streamfunction, dpsi_dx);      // ∂ψ/∂x
    dx_matrix.Mult(vorticity, domega_dx);         // ∂ω/∂x
    dy_matrix.Mult(vorticity, domega_dy);         // ∂ω/∂y
    dxx_matrix.Mult(vorticity, d2omega_dx2);      // ∂²ω/∂x²
    dyy_matrix.Mult(vorticity, d2omega_dy2);      // ∂²ω/∂y²

    // Update vorticity using explicit Euler scheme (Equation 11 from paper)
    // ω^(n+1) = ω^n + dt*[(1/Re)*∇²ω - (∂ψ/∂y)(∂ω/∂x) + (∂ψ/∂x)(∂ω/∂y)]
    for (int i = 0; i < vorticity.Size(); ++i) {
        double convection = dpsi_dy[i] * domega_dx[i] - dpsi_dx[i] * domega_dy[i];
        double diffusion = (1.0 / reynolds_number) * (d2omega_dx2[i] + d2omega_dy2[i]);
        vorticity[i] += dt * (diffusion - convection);
    }
}

/**
 * @brief Apply vorticity boundary conditions for lid-driven cavity
 *
 * Implements Equation 33 from Bourantas et al. (2019): ω = ∂v/∂x - ∂u/∂y
 * First applies velocity BCs (u=v=0 on walls, u=1 on top lid), then
 * computes boundary vorticity from the velocity field.
 */
void ApplyBoundaryConditions(const std::vector<int>& bottom_nodes, const std::vector<int>& right_nodes,
                           const std::vector<int>& top_nodes, const std::vector<int>& left_nodes,
                           const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                           const mfem::Vector& streamfunction, mfem::Vector& vorticity)
{
    // Compute velocity components from streamfunction
    mfem::Vector u(streamfunction.Size());  // u = ∂ψ/∂y
    mfem::Vector v(streamfunction.Size());  // v = -∂ψ/∂x

    dy_matrix.Mult(streamfunction, u);
    dx_matrix.Mult(streamfunction, v);
    v *= -1.0;

    // Apply velocity boundary conditions for lid-driven cavity
    for (int idx : bottom_nodes) u[idx] = 0.0;  // u = 0 on bottom wall
    for (int idx : right_nodes)  u[idx] = 0.0;  // u = 0 on right wall
    for (int idx : left_nodes)   u[idx] = 0.0;  // u = 0 on left wall
    for (int idx : top_nodes)    u[idx] = 1.0;  // u = 1 on top wall (lid-driven)

    for (int idx : bottom_nodes) v[idx] = 0.0;  // v = 0 on all walls
    for (int idx : right_nodes)  v[idx] = 0.0;
    for (int idx : top_nodes)    v[idx] = 0.0;
    for (int idx : left_nodes)   v[idx] = 0.0;

    // Compute boundary vorticity: ω = ∂v/∂x - ∂u/∂y (Equation 33 from paper)
    mfem::Vector du_dy(u.Size());
    mfem::Vector dv_dx(v.Size());

    dy_matrix.Mult(u, du_dy);
    dx_matrix.Mult(v, dv_dx);

    // Apply vorticity boundary conditions on all boundary nodes
    auto apply_bc_to_nodes = [&](const std::vector<int>& nodes) {
        for (int idx : nodes) {
            vorticity[idx] = dv_dx[idx] - du_dy[idx];
        }
    };

    apply_bc_to_nodes(bottom_nodes);
    apply_bc_to_nodes(right_nodes);
    apply_bc_to_nodes(top_nodes);
    apply_bc_to_nodes(left_nodes);
}


/**
 * @brief Compute critical time step using Gershgorin circle theorem
 *
 * Implements Equation 42 from Bourantas et al. (2019) to estimate the maximum
 * eigenvalue of the linearized system matrix. Returns dt ≤ 2/|λ_max| for stability.
 */
double ComputeGershgorinTimeStep(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                                const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                                const mfem::Vector& streamfunction, double reynolds_number) {
    // Implement Gershgorin circle theorem estimation (Equation 42 from paper)
    double max_eigenvalue_bound = 0.0;
    const int num_nodes = streamfunction.Size();

    // Compute velocity components
    mfem::Vector u(num_nodes), v(num_nodes);
    dy_matrix.Mult(streamfunction, u);   // u = ∂ψ/∂y
    dx_matrix.Mult(streamfunction, v);   // v = -∂ψ/∂x
    v *= -1.0;

    // For each node, compute the Gershgorin bound
    for (int i = 0; i < num_nodes; ++i) {
        double row_sum = 0.0;

        // Diffusion terms: (1/Re) * (∂²/∂x² + ∂²/∂y²)
        for (int j = dxx_matrix.GetI()[i]; j < dxx_matrix.GetI()[i + 1]; ++j) {
            row_sum += std::abs((1.0 / reynolds_number) * dxx_matrix.GetData()[j]);
        }
        for (int j = dyy_matrix.GetI()[i]; j < dyy_matrix.GetI()[i + 1]; ++j) {
            row_sum += std::abs((1.0 / reynolds_number) * dyy_matrix.GetData()[j]);
        }

        // Convection terms: u*(∂/∂x) + v*(∂/∂y)
        for (int j = dx_matrix.GetI()[i]; j < dx_matrix.GetI()[i + 1]; ++j) {
            row_sum += std::abs(u[i] * dx_matrix.GetData()[j]);
        }
        for (int j = dy_matrix.GetI()[i]; j < dy_matrix.GetI()[i + 1]; ++j) {
            row_sum += std::abs(v[i] * dy_matrix.GetData()[j]);
        }

        max_eigenvalue_bound = std::max(max_eigenvalue_bound, row_sum);
    }

    // Critical time step: dt ≤ 2/|λ_max|
    return (max_eigenvalue_bound > 0.0) ? 2.0 / max_eigenvalue_bound : 1e-4;
}

/**
 * @brief Save vorticity and streamfunction fields to data files
 *
 * Exports solution fields as text files with one value per line for analysis
 * in MATLAB or other post-processing tools.
 */
void SaveSolutionToFile(const mfem::Vector& vorticity, const mfem::Vector& streamfunction,
                       const std::string& filename, int timestep, std::string dat_dir) {
    //std::string vort_filename = filename + "_vorticity_" + std::to_string(timestep) + ".dat";
    //std::string stream_filename = filename + "_streamfunction_" + std::to_string(timestep) + ".dat";
    std::string vort_filename = dat_dir + "/" + filename + "_vorticity_" + ".dat";
    std::string stream_filename = dat_dir + "/" + filename + "_streamfunction_" +  ".dat";

    // Save vorticity field
    std::ofstream vort_file(vort_filename);
    if (vort_file.is_open()) {
        vort_file.precision(12);
        for (int i = 0; i < vorticity.Size(); ++i) {
            vort_file << vorticity[i] << std::endl;
        }
        vort_file.close();
    }

    // Save streamfunction field
    std::ofstream stream_file(stream_filename);
    if (stream_file.is_open()) {
        stream_file.precision(12);
        for (int i = 0; i < streamfunction.Size(); ++i) {
            stream_file << streamfunction[i] << std::endl;
        }
        stream_file.close();
    }

    std::cout << "SaveSolutionToFile: Saved solution at timestep " << timestep << std::endl;
}
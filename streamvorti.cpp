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
//     ./StreamVorti -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -nn 25 -Re 1000 -solver umfpack -sm -sn -sd -sdd -ss -pv

#include <StreamVorti/stream_vorti.hpp>
#include "mfem.hpp"
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif


// Stream-Vorti simulation (vorticity, Gersc time step, boundary condition)
// 1. get DCPSE derivative matrices - DONE
// 2. Define boundaries (use MFEM method) - DONE
// Initialize vorticity and stream function -DONE
// set up Linear solver (Mfem example)
// 3. update vorticty (eq. 11)
// apply boundary condition
// 4. compute Gershgorin time step (2.6)
// 5. solve streamfunction (eq. 12)
// 6. update velocity (eq. 13)
// run simulation (main time steps for loop)
// save solution

// Simulation parameters structure
struct SimulationParams {
    double final_time = 60000.0;
    double dt = 1e-2;
    double reynolds_number = 1000.0;

    std::string output_prefix = "mfem_square10x10";
    std::string output_extension = ".dat";

    // Solution output
    int output_frequency = 100;

    // Adaptive timestep options
    int gershgorin_frequency = 1000;
    bool enable_adaptive_timestep = true;
    double gershgorin_safety_factor = 5.0;

    // Linear solver options
    std::string linear_solver = "umfpack";
    int solver_max_iter = 500;
    double solver_rel_tol = 1e-12;
    int solver_print_level = 0;

    // Convergence parameters
    double steady_state_tol = 1e-6;
    int steady_state_freq = 1000;
    int steady_state_checks = 3;
};

// Function declarations
// Refactored functions
mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh);
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors);
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                            int dim, bool save_d, bool save_dd, std::string dat_dir);

// Lid-driven cavity boundaries
void IdentifyBoundaryNodesLDC(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                            std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                            std::vector<int>& left_nodes, std::vector<int>& interior_nodes);


// simulation
void UpdateVorticity(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                    const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                    const mfem::Vector& streamfunction, mfem::Vector& vorticity,
                    double dt, double reynolds_number);
void ApplyBoundaryConditions(const std::vector<int>& bottom_nodes, const std::vector<int>& right_nodes,
                           const std::vector<int>& top_nodes, const std::vector<int>& left_nodes,
                           const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                           const mfem::Vector& streamfunction, mfem::Vector& vorticity);
double ComputeGershgorinTimeStep(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                                const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                                const mfem::Vector& streamfunction, double reynolds_number);
void SaveSolutionToFile(const mfem::Vector& vorticity, const mfem::Vector& streamfunction, const std::string& filename, int timestep, std::string dat_dir);
mfem::Solver* CreateLinearSolver(const std::string& solver_type, const mfem::SparseMatrix& A, const SimulationParams& params);
bool CheckSteadyState(const mfem::Vector& vorticity, const mfem::Vector& streamfunction,
                      int timestep, double dt,
                      double tolerance, int check_freq, int required_checks);


int main(int argc, char *argv[])
{
    // Options
    const char *mesh_file = "";
    int order = 1;
    // bool visualization = 1;

    // Output filename prefix and extension
    std::string fname = "mfem_square10x10";
    std::string fext = ".dat";

    // DC PSE parameters
    double reynolds_number = 1000.0;
    int NumNeighbors = 25;

    // Mesh generation
    int dim = 2;
    int nx = 20;
    int ny = 20;
    int nz = 0;
    double sx = 1.0;
    double sy = 1.0;
    double sz = 1.0;

    // Output requests
    bool save_mesh = false;
    bool save_neighbors = false;
    bool save_d = false;        // all 1st derivatives (gradient)
    bool save_dd = false;       // all 2nd derivatives (Hessian)
    bool save_solutions = false;

    // Solver options
    std::string solver_type = "umfpack";
    int solver_max_iter = 500;
    double solver_rel_tol = 1e-10;
    int solver_print_level = 1;

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
    // args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
    //                "--no-visualization",
    //                "Enable or disable GLVis visualization.");
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
    args.AddOption(&paraview_output, "-pv", "--paraview", "-no-pv",
                  "--no-paraview",
                  "Enable or disable ParaView output.");
    args.AddOption(&reynolds_number, "-Re", "--reynolds-number",
                    "Reynolds number of simulating fluid.");
    args.AddOption(&solver_type, "-solver", "--linear-solver",
                   "Linear solver: umfpack, cg, gmres, minres, bicgstab, hypre (if available)");
    args.AddOption(&solver_max_iter, "-maxit", "--max-iterations",
                   "Maximum iterations for iterative solvers");
    args.AddOption(&solver_rel_tol, "-tol", "--solver-tolerance",
                   "Relative tolerance for iterative solvers");
    args.AddOption(&solver_print_level, "-spl", "--solver-print-level",
                   "Solver print level (0=quiet, 1=summary, 2=detailed)");
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


    /*********************** Simulation  Setup  *************************/

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


    /*********************  Streamfunction Initialization   ***************************/
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
    // Create Laplacian and flip sign
    mfem::SparseMatrix* laplacian_matrix = new mfem::SparseMatrix(dxx_matrix);
    laplacian_matrix->Add(1.0, dyy_matrix);
    *laplacian_matrix *= -1.0;

    // Pre-sort boundary nodes for cache efficiency
    std::vector<int> boundary_sorted = all_boundary_nodes;
    std::sort(boundary_sorted.begin(), boundary_sorted.end());

    // Apply BC
    for (int idx : all_boundary_nodes) {
        laplacian_matrix->EliminateRow(idx);  // Zeros the row except diagonal
        laplacian_matrix->Set(idx, idx, 1.0); // Set diagonal
    }
    laplacian_matrix->Finalize();
    std::cout << "RunSimulation: Created Laplacian matrix with boundary conditions." << std::endl;

    // Set up linear solver
    params.linear_solver = solver_type;
    params.solver_max_iter = solver_max_iter;
    params.solver_rel_tol = solver_rel_tol;

    mfem::Solver* linear_solver = CreateLinearSolver(params.linear_solver, *laplacian_matrix, params);
    if (!linear_solver) {
        MFEM_ABORT("Failed to create linear solver");
    }

    mfem::Vector rhs;

    /*********************** Paraview Output Setup *************************/

    // Calculate velocity from streamfunction
    mfem::Vector u_velocity(streamfunction.Size());
    mfem::Vector v_velocity(streamfunction.Size());

    dy_matrix.Mult(streamfunction, u_velocity);   // u = ∂ψ/∂y
    dx_matrix.Mult(streamfunction, v_velocity);   // v = -∂ψ/∂x
    v_velocity *= -1.0;

    // Create finite element spaces
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

    /*********************** Simulation Timestepping Loop *************************/

    std::cout << "RunSimulation: Running " << num_timesteps << " time steps..." << std::endl;
    mfem::StopWatch sim_timer;
    sim_timer.Start();

    // #ifdef _OPENMP
    // int max_threads = omp_get_max_threads();
    // std::cout << "Using OpenMP with " << max_threads << " threads" << std::endl;
    // #else
    // std::cout << "OpenMP NOT enabled (serial execution)" << std::endl;
    // #endif

    // volatile bool flag=false;

    // #pragma omp parallel for schedule(dynamic, 64)
    // Main simulation loop (implementing Algorithm from Bourantas et al. 2019)
    for (int time_step = 1; time_step <= num_timesteps; ++time_step) {
        // if(flag) continue;
        double current_time = time_step * current_dt;

        // Step 1: Update vorticity using explicit Euler scheme (Equation 11)
        UpdateVorticity(dx_matrix, dy_matrix, dxx_matrix, dyy_matrix, streamfunction, vorticity, current_dt, params.reynolds_number);

        // Step 2: Solve Poisson equation for streamfunction (Equation 12)
        // // Set up the Poisson equation: ∇²ψ = -ω (Equation 12 from paper)
        // Solve Poisson equation: -∇²ψ = ω
        rhs = vorticity;
        //rhs *= -1.0;

        // Apply homogeneous Dirichlet boundary conditions: ψ = 0 on all boundaries
        for (int idx : all_boundary_nodes) {
            rhs[idx] = 0.0;
        }
        // Time the solve
        // mfem::StopWatch solve_timer;
        // solve_timer.Start();

        linear_solver->Mult(rhs, streamfunction);
        // Check if solver converged
        if (auto* iterative_solver = dynamic_cast<mfem::IterativeSolver*>(linear_solver)) {
            if (!iterative_solver->GetConverged()) {
                std::cerr << "Warning: Iterative Solver did not converge at timestep " << time_step
                        << " (iterations: " << iterative_solver->GetNumIterations() << ")" << std::endl;
            }
        }

        // solve_timer.Stop();
        // std::cout << "Solve time: " << solve_timer.RealTime() << " seconds" << std::endl;

        // Step 3: Apply vorticity boundary conditions (Equation 33)
        ApplyBoundaryConditions(bottom_nodes, right_nodes, top_nodes, left_nodes,dx_matrix, dy_matrix, streamfunction, vorticity);

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

        // Step 6: Save Paraview output (optional)
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
        if (CheckSteadyState(vorticity, streamfunction, time_step, current_dt,
                            params.steady_state_tol,
                            params.steady_state_freq,
                            params.steady_state_checks)) {
            std::cout << "Early termination - steady state reached" << std::endl;
            break;
            // flag=true;
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

    // Free the used memory
    delete derivs;
    delete mesh;
    delete linear_solver;
    delete laplacian_matrix;

    std::cout << "main: success!" << std::endl;
    return EXIT_SUCCESS;
}


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


// Functions matching MATLAB  "Explicit_streamfunction_vorticity_meshless.m"

// Identify boundary ndoes for Lid-Driven Cavity problem
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

    // std::cout << "SaveSolutionToFile: Saved solution at timestep " << timestep << std::endl;
}


// Create appropriate linear solver based on user selection
// Options: umfpack, klu, cg, minres, gmres, bicgstab, hypre
mfem::Solver* CreateLinearSolver(const std::string& solver_type,
                                  const mfem::SparseMatrix& A,
                                  const SimulationParams& params)
{
    mfem::Solver* solver = nullptr;

    std::cout << "Creating linear solver: " << solver_type << std::endl;

    if (solver_type == "umfpack") {
#ifdef MFEM_USE_SUITESPARSE
        std::cout << "Using UMFPackSolver (direct solver)" << std::endl;
        mfem::UMFPackSolver* umf_solver = new mfem::UMFPackSolver();
        umf_solver->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
        umf_solver->SetOperator(A);
        solver = umf_solver;
#else
        std::cerr << "ERROR: UMFPACK not available. Rebuild MFEM with SuiteSparse." << std::endl;
        return nullptr;
#endif
    }
    else if (solver_type == "klu") {
#ifdef MFEM_USE_SUITESPARSE
        std::cout << "Using KLUSolver (direct solver)" << std::endl;
        mfem::KLUSolver* klu_solver = new mfem::KLUSolver();
        klu_solver->SetOperator(A);
        solver = klu_solver;
#else
        std::cerr << "ERROR: KLU not available. Rebuild MFEM with SuiteSparse." << std::endl;
        return nullptr;
#endif
    }
    else if (solver_type == "cg") {
        std::cout << "Using CGSolver (iterative, symmetric positive definite)" << std::endl;
        mfem::CGSolver* cg_solver = new mfem::CGSolver();
        cg_solver->SetRelTol(params.solver_rel_tol);
        cg_solver->SetMaxIter(params.solver_max_iter);
        cg_solver->SetPrintLevel(params.solver_print_level);

        // Optional: Add preconditioner
        mfem::GSSmoother* precond = new mfem::GSSmoother();
        cg_solver->SetPreconditioner(*precond);
        cg_solver->SetOperator(A);
        solver = cg_solver;
    }
    else if (solver_type == "gmres") {
        std::cout << "Using GMRESSolver (iterative, general matrices)" << std::endl;
        mfem::GMRESSolver* gmres_solver = new mfem::GMRESSolver();
        gmres_solver->SetRelTol(params.solver_rel_tol);
        gmres_solver->SetMaxIter(params.solver_max_iter);
        gmres_solver->SetPrintLevel(params.solver_print_level);
        gmres_solver->SetKDim(200); // Restart parameter

        // Optional: Add preconditioner
        // Option 3: Incomplete LU (better than op1 & 2)
#ifdef MFEM_USE_SUITESPARSE
        mfem::UMFPackSolver* ilu_precond = new mfem::UMFPackSolver();
        ilu_precond->SetOperator(A);
        gmres_solver->SetPreconditioner(*ilu_precond);
#else
        // Option 1: Gauss-Seidel (simple, effective)
        mfem::GSSmoother* gs_precond = new mfem::GSSmoother();
        gmres_solver->SetPreconditioner(*gs_precond);

        // Option 2: Jacobi (parallel-friendly)
        // mfem::DSmoother* jacobi_precond = new mfem::DSmoother();
#endif
        gmres_solver->SetOperator(A);
        solver = gmres_solver;
    }
    else if (solver_type == "minres") {
        std::cout << "Using MINRESSolver (iterative, symmetric indefinite)" << std::endl;
        mfem::MINRESSolver* minres_solver = new mfem::MINRESSolver();
        minres_solver->SetRelTol(params.solver_rel_tol);
        minres_solver->SetMaxIter(params.solver_max_iter);
        minres_solver->SetPrintLevel(params.solver_print_level);
        minres_solver->SetOperator(A);
        solver = minres_solver;
    }
    else if (solver_type == "bicgstab") {
        std::cout << "Using BiCGSTABSolver (iterative, nonsymmetric)" << std::endl;
        mfem::BiCGSTABSolver* bicg_solver = new mfem::BiCGSTABSolver();
        bicg_solver->SetRelTol(params.solver_rel_tol);
        bicg_solver->SetAbsTol(1e-14);
        bicg_solver->SetMaxIter(params.solver_max_iter);
        bicg_solver->SetPrintLevel(params.solver_print_level);

        // Optional: Add preconditioner
        mfem::GSSmoother* precond = new mfem::GSSmoother();
        bicg_solver->SetPreconditioner(*precond);
        bicg_solver->SetOperator(A);
        solver = bicg_solver;
    }
    else if (solver_type == "hypre") { // Todo: HypreBoomerAMG does not work (need HypreParMatrix)
#ifdef MFEM_USE_HYPRE
        std::cout << "Using HypreBoomerAMG (parallel algebraic multigrid)" << std::endl;
        mfem::HypreBoomerAMG* amg = new mfem::HypreBoomerAMG();
        amg->SetPrintLevel(params.solver_print_level);
        amg->SetOperator(A);
        solver = amg;
#else
        std::cerr << "ERROR: HYPRE not available. Rebuild MFEM with HYPRE support." << std::endl;
        return nullptr;
#endif
    }
    else if (solver_type == "amg") {
        std::cout << "Using serial geometric multigrid" << std::endl;

        // Use MFEM's serial multigrid preconditioner with CG
        mfem::CGSolver* cg = new mfem::CGSolver();

        // Simple Jacobi preconditioner (serial alternative to AMG)
        mfem::DSmoother* jacobi = new mfem::DSmoother(A);

        cg->SetPreconditioner(*jacobi);
        cg->SetRelTol(params.solver_rel_tol);
        cg->SetMaxIter(params.solver_max_iter);
        cg->SetPrintLevel(params.solver_print_level);
        cg->SetOperator(A);

        solver = cg;
    }
    else {
        std::cerr << "ERROR: Unknown solver type: " << solver_type << std::endl;
        std::cerr << "Available solvers: umfpack, klu, cg, gmres, minres, bicgstab, hypre, amg" << std::endl;
        return nullptr;
    }

    return solver;
}

// Simple convergence check function
bool CheckSteadyState(const mfem::Vector& vorticity, const mfem::Vector& streamfunction,
                      int timestep, double dt,
                      double tolerance, int check_freq, int required_checks) {

    // Static variables to maintain state between calls
    static mfem::Vector prev_vorticity;
    static mfem::Vector prev_streamfunction;
    static int consecutive_passes = 0;
    static bool initialized = false;

    // Only check at specified intervals and after initial period
    if (timestep < 100 || timestep % check_freq != 0) {
        return false;
    }

    // Initialize on first check
    if (!initialized) {
        prev_vorticity = vorticity;
        prev_streamfunction = streamfunction;
        initialized = true;
        std::cout << "Steady state monitoring started at timestep " << timestep << std::endl;
        return false;
    }

    // Compute relative changes
    mfem::Vector diff_vort(vorticity.Size());
    mfem::Vector diff_stream(streamfunction.Size());

    diff_vort = vorticity;
    diff_vort -= prev_vorticity;

    diff_stream = streamfunction;
    diff_stream -= prev_streamfunction;

    double vort_change = diff_vort.Norml2() / std::max(vorticity.Norml2(), 1e-12);
    double stream_change = diff_stream.Norml2() / std::max(streamfunction.Norml2(), 1e-12);
    double max_change = std::max(vort_change, stream_change);

    // // Report progress
    // std::cout << "Step " << timestep << " (t=" << std::fixed << std::setprecision(6)
    //           << timestep * dt << "): ‖Δ‖/‖·‖ = " << std::scientific
    //           << max_change;

    // Check convergence
    if (max_change < tolerance) {
        consecutive_passes++;
        //std::cout << " ✓ [" << consecutive_passes << "/" << required_checks << "]" << std::endl;

        if (consecutive_passes >= required_checks) {
            std::cout << "\n" << std::string(50, '=') << std::endl;
            std::cout << "CONVERGED TO STEADY STATE" << std::endl;
            std::cout << std::string(50, '=') << std::endl;
            std::cout << "Final relative change: " << max_change << std::endl;
            std::cout << "Convergence tolerance: " << tolerance << std::endl;
            std::cout << "Total timesteps: " << timestep << std::endl;
            std::cout << "Physical time: " << timestep * dt << std::endl;
            std::cout << std::string(50, '=') << "\n" << std::endl;
            return true;
        }
    } else {
        //std::cout << " ✗" << std::endl;
        consecutive_passes = 0;  // Reset counter
    }

    // Update previous values
    prev_vorticity = vorticity;
    prev_streamfunction = streamfunction;

    return false;
}
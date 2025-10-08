/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 * Copyright (C) 2025 Weizheng Li
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
 *      Weizheng Li
 */

// Demo usage: (save everything)
//     ./StreamVorti -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -nn 25 -Re 1000 -solver umfpack -sm -sn -sd -sdd -ss -pv -cr

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
#include <chrono>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

// Simulation parameters structure
struct SimulationParams {
    double final_time = 60000.0;
    double dt = 1e-2;
    double reynolds_number = 1000.0;

    std::string output_prefix = "mfem_square10x10";
    std::string output_extension = ".dat";

    // Solution output
    int output_frequency = 100;
    int paraview_output_frequency = 1000;  // Default value

    // Adaptive timestep options
    int gershgorin_frequency = 1000;
    bool enable_adaptive_timestep = true;
    double gershgorin_safety_factor = 5.0;

    // Linear solver options
    std::string solver_type = "umfpack";
    int solver_max_iter = 500;
    double solver_rel_tol = 1e-12;
    int solver_print_level = 0;

    // Convergence parameters
    double steady_state_tol = 1e-6;
    int steady_state_freq = 500;
    int steady_state_checks = 2;

    // Residual monitoring
    bool check_residuals = false;
    int residual_check_freq = 100;

    // Grid-adaptive timestep computation
    void ComputeStableTimestep(int num_nodes, double reynolds) {
        // Estimate grid size (assuming square grid)
        int N = static_cast<int>(std::sqrt(num_nodes));
        double h = 1.0 / (N - 1);  // Grid spacing for unit square

        // CFL condition for advection
        double dt_cfl = h;  // Since U_max = 1.0 for lid-driven cavity

        // Diffusion stability limit
        double dt_diffusion = 0.25 * h * h * reynolds;

        // Use more restrictive limit
        double dt_stability = std::min(dt_cfl, dt_diffusion);

        // Apply safety factor for explicit scheme
        double safety_factor = 0.5;
        double dt_recommended = safety_factor * dt_stability;

        // Update if current dt is too large
        if (dt > dt_recommended) {
            std::cout << "WARNING: Default dt=" << dt << " exceeds stability limit!" << std::endl;
            std::cout << "  Grid: " << N << "x" << N << " (h=" << h << ")" << std::endl;
            std::cout << "  CFL limit: " << dt_cfl << std::endl;
            std::cout << "  Diffusion limit: " << dt_diffusion << std::endl;
            std::cout << "  Adjusting dt to: " << dt_recommended << std::endl;
            dt = dt_recommended;
        }

        // Additional safety for very fine grids
        if (N >= 200) {
            // Based on paper's empirical values
            dt = std::min(dt, 5e-3);
            std::cout << "Fine grid detected. Limiting dt to: " << dt << std::endl;
        }
    }
};

// Structure to hold residual history
struct ResidualHistory {
    std::vector<int> timesteps;
    std::vector<double> times;
    std::vector<double> vorticity_rms;
    std::vector<double> vorticity_max;
    std::vector<double> streamfunction_rms;
    std::vector<double> streamfunction_max;
    std::vector<double> vorticity_change;
    std::vector<double> streamfunction_change;

    void Add(int step, double time,
             double vort_rms, double vort_max,
             double psi_rms, double psi_max,
             double vort_rel_change = 0.0,
             double psi_rel_change = 0.0) {
        timesteps.push_back(step);
        times.push_back(time);
        vorticity_rms.push_back(vort_rms);
        vorticity_max.push_back(vort_max);
        streamfunction_rms.push_back(psi_rms);
        streamfunction_max.push_back(psi_max);
        vorticity_change.push_back(vort_rel_change);
        streamfunction_change.push_back(psi_rel_change);
    }

    void SaveToCSV(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        file << "timestep,time,vorticity_rms,vorticity_max,"
             << "streamfunction_rms,streamfunction_max,"
             << "vorticity_rel_change,streamfunction_rel_change\n";

        file << std::scientific << std::setprecision(12);

        for (size_t i = 0; i < timesteps.size(); ++i) {
            file << timesteps[i] << ","
                 << times[i] << ","
                 << vorticity_rms[i] << ","
                 << vorticity_max[i] << ","
                 << streamfunction_rms[i] << ","
                 << streamfunction_max[i] << ","
                 << vorticity_change[i] << ","
                 << streamfunction_change[i] << "\n";
        }

        file.close();
        std::cout << "Residual history saved to: " << filename << std::endl;
    }

    void PrintSummary() {
        if (timesteps.empty()) {
            std::cout << "No residual data collected." << std::endl;
            return;
        }

        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "RESIDUAL HISTORY SUMMARY" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        std::cout << "Total data points: " << timesteps.size() << std::endl;
        std::cout << "Initial vorticity RMS: " << std::scientific
                  << vorticity_rms.front() << std::endl;
        std::cout << "Final vorticity RMS: " << vorticity_rms.back() << std::endl;
        std::cout << "Initial streamfunction RMS: " << streamfunction_rms.front() << std::endl;
        std::cout << "Final streamfunction RMS: " << streamfunction_rms.back() << std::endl;

        if (timesteps.size() > 10) {
            size_t n = timesteps.size();
            double last_vort = vorticity_rms.back();
            double prev_vort = vorticity_rms[n-10];
            double vort_reduction = (prev_vort - last_vort) / prev_vort;

            if (vort_reduction < 0.01) {
                std::cout << "\n⚠ WARNING: Residuals appear stagnant!" << std::endl;
                std::cout << "   Reduction over last 10 checks: "
                          << (vort_reduction * 100.0) << "%" << std::endl;
            }
        }
        std::cout << std::string(70, '=') << "\n" << std::endl;
    }
};

// Structure to hold timing results
struct PerformanceMetrics {
    double derivative_time = 0.0;
    double factorization_time = 0.0;
    double solution_time_per_iter = 0.0;
    int grid_size = 0;
    int total_iterations = 0;

    void PrintPerformanceTable(const PerformanceMetrics& metrics) {
        std::cout << "\n";
        std::cout << "========================================================================\n";
        std::cout << "                     PERFORMANCE METRICS TABLE                          \n";
        std::cout << "========================================================================\n";
        std::cout << std::setw(15) << "Grid Resolution"
                << std::setw(15) << "Derivatives"
                << std::setw(20) << "Factorization"
                << std::setw(20) << "Time/Iteration" << "\n";
        std::cout << "------------------------------------------------------------------------\n";

        // Print the single row for this run
        std::cout << std::setw(10) << std::to_string(metrics.grid_size) + " x " + std::to_string(metrics.grid_size)
                << std::setw(17) << std::fixed << std::setprecision(4) << metrics.derivative_time << " s"
                << std::setw(17) << std::fixed << std::setprecision(4) << metrics.factorization_time << " s"
                << std::setw(20) << std::fixed << std::setprecision(4) << metrics.solution_time_per_iter << " s" << "\n";
        std::cout << "========================================================================\n\n";
    }
};

// Function declarations

mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh);
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors);
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                            int dim, bool save_d, bool save_dd, std::string dat_dir);

// Lid-driven cavity boundaries
void IdentifyBoundaryNodes(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                            std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                            std::vector<int>& left_nodes, std::vector<int>& interior_nodes);

double ComputeGershgorinTimeStep(const mfem::SparseMatrix& dx_matrix, const mfem::SparseMatrix& dy_matrix,
                                const mfem::SparseMatrix& dxx_matrix, const mfem::SparseMatrix& dyy_matrix,
                                const mfem::Vector& streamfunction, double reynolds_number);

mfem::Solver* CreateLinearSolver(const std::string& solver_type, const mfem::SparseMatrix& A, const SimulationParams& params);

void SaveSolutionToFile(const mfem::Vector& vorticity, const mfem::Vector& streamfunction,
                       const mfem::Vector& u_velocity, const mfem::Vector& v_velocity,
                       const std::string& filename, int timestep,
                       std::string dat_dir, bool is_final);
// Helper function
std::string GetCurrentDateTime();



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

    bool paraview_output = false;
    std::string paraview_filename = fname;

    // Simulation parameters
    SimulationParams params;
    // Performance Metrics parameters
    PerformanceMetrics perf_metrics;

    // Parse command-line options
    mfem::OptionsParser args(argc, argv);
    // Similation options
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
    args.AddOption(&NumNeighbors, "-nn", "--num-neighbors",
                    "Number of neighbors for DCPSE (default: 25).");
    args.AddOption(&reynolds_number, "-Re", "--reynolds-number",
                    "Reynolds number of simulating fluid (default: 1000).");
    // Output options
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
    args.AddOption(&paraview_output, "-pv", "--paraview", "-no-pv", "--no-paraview",
                  "Enable or disable ParaView output.");
    args.AddOption(&params.paraview_output_frequency, "-pvf", "--paraview-freq",
               "Frequency of ParaView output (in timesteps).");
    // Timesteping options
    args.AddOption(&params.gershgorin_frequency, "-gf", "--gershgorin-freq",
                   "Frequency of Gershgorin timestep checks (default: 1000).");
    args.AddOption(&params.enable_adaptive_timestep,
                   "-adaptive", "--adaptive-timestep", "-no-adaptive", "--no-adaptive-timestep",
                   "Enable or disable adaptive timestepping based on Gershgorin criterion.");
    args.AddOption(&params.gershgorin_safety_factor, "-gsf", "--gershgorin-safety",
                   "Safety factor for Gershgorin timestep (default: 5.0).");
    // Solver options
    args.AddOption(&params.solver_type, "-solver", "--linear-solver",
                   "Linear solver: umfpack, cg, gmres, minres, bicgstab, hypre (if available)");
    args.AddOption(&params.solver_print_level, "-spl", "--solver-print-level",
                   "Solver print level (0=quiet, 1=summary, 2=verbose).");
    args.AddOption(&params.solver_max_iter, "-smi", "--solver-max-iter",
                   "Maximum iterations for iterative solvers.");
    args.AddOption(&params.solver_rel_tol, "-srt", "--solver-rel-tol",
                   "Relative tolerance for iterative solvers.");
    // Residual monitoring options
    args.AddOption(&params.check_residuals,
                "-cr", "--check-residuals", "-no-cr", "--no-check-residuals",
                "Enable or disable residual monitoring and logging.");
    args.AddOption(&params.residual_check_freq, "-rcf", "--residual-check-freq",
                "Frequency of residual checks (default: 100).");
    args.Parse();

    if (!args.Good()) {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    // Create or load mesh
    mfem::Mesh* mesh = CreateOrLoadMesh(mesh_file, dim, nx, ny, nz, sx, sy, sz, save_mesh);
    perf_metrics.grid_size = nx;

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

// ====================================================================
// SIMULATION SETUP
// ====================================================================
    const int num_nodes = mesh->GetNV();

    // Ensure we have a 2D DCPSE object first
    StreamVorti::Dcpse2d* dcpse2d = dynamic_cast<StreamVorti::Dcpse2d*>(derivs);
    if (!dcpse2d) {
        MFEM_ABORT("Setup: Only 2D simulations are currently supported.");
    }

    // ================== Timing PHASE 1: DERIVATIVE COMPUTATION ==================
    mfem::StopWatch deriv_timer;
    deriv_timer.Start();

    // Get DCPSE derivative matrices
    const mfem::SparseMatrix& dx_matrix = dcpse2d->ShapeFunctionDx();
    const mfem::SparseMatrix& dy_matrix = dcpse2d->ShapeFunctionDy();
    const mfem::SparseMatrix& dxx_matrix = dcpse2d->ShapeFunctionDxx();
    const mfem::SparseMatrix& dyy_matrix = dcpse2d->ShapeFunctionDyy();

    std::cout << "Setup: Retrieved DCPSE derivative matrices successfully." << std::endl;

    deriv_timer.Stop();
    perf_metrics.derivative_time = deriv_timer.RealTime();

    std::cout << "Phase 1 - Derivatives computation: " << perf_metrics.derivative_time << " s" << std::endl;

    // Identify boundary and interior nodes
    std::vector<int> bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes;
    IdentifyBoundaryNodes(mesh, bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes);

    // Combine all boundary nodes for streamfunction boundary conditions
    std::vector<int> all_boundary_nodes;
    all_boundary_nodes.insert(all_boundary_nodes.end(), bottom_nodes.begin(), bottom_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), right_nodes.begin(), right_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), top_nodes.begin(), top_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), left_nodes.begin(), left_nodes.end());

    // Create Laplacian matrix for streamfunction equation
    // Flip sign to Ensures matrix positive definiteness for iterative solvers
    mfem::SparseMatrix* laplacian_matrix = new mfem::SparseMatrix(dxx_matrix);
    laplacian_matrix->Add(1.0, dyy_matrix);
    *laplacian_matrix *= -1.0;
    // Apply BC
    for (int idx : all_boundary_nodes) {
        laplacian_matrix->EliminateRow(idx);  // Zeros the row except diagonal
        laplacian_matrix->Set(idx, idx, 1.0); // Set diagonal
    }
    laplacian_matrix->Finalize();

    std::cout << "Setup: Created Laplacian matrix with boundary conditions." << std::endl;

    // ================== Timing PHASE 2: FACTORIZATION ==================
    mfem::StopWatch factor_timer;
    factor_timer.Start();

    // Create and setup linear solver (factorization happens here for direct solvers)
    mfem::Solver* linear_solver = CreateLinearSolver(params.solver_type, *laplacian_matrix, params);

    if (!linear_solver) {
        MFEM_ABORT("Failed to create linear solver");
    }

    factor_timer.Stop();
    perf_metrics.factorization_time = factor_timer.RealTime();
    std::cout << "Phase 2 - Factorization: " << perf_metrics.factorization_time << " s" << std::endl;

    // Initialize solution vectors
    mfem::Vector vorticity(num_nodes);
    mfem::Vector streamfunction(num_nodes);
    mfem::Vector vorticity_old(num_nodes); // for residual monitor
    mfem::Vector rhs(num_nodes);

    vorticity = 0.0;
    streamfunction = 0.0;

    // Derivative vectors (allocated once, reused throughout)
    mfem::Vector dpsi_dy(num_nodes);      // ∂ψ/∂y = u-velocity =
    mfem::Vector dpsi_dx(num_nodes);      // ∂ψ/∂x (needs negation for v-velocity)
    mfem::Vector domega_dx(num_nodes);
    mfem::Vector domega_dy(num_nodes);
    mfem::Vector d2omega_dx2(num_nodes);
    mfem::Vector d2omega_dy2(num_nodes);

// ====================================================================
// PARAVIEW OUTPUT SETUP
// ====================================================================
    // Create finite element spaces
    mfem::FiniteElementSpace scalar_fes(const_cast<mfem::Mesh*>(mesh), &fec, 1);
    mfem::FiniteElementSpace vector_fes(const_cast<mfem::Mesh*>(mesh), &fec, mesh->Dimension());

    // Create grid functions
    mfem::GridFunction vorticity_gf(&scalar_fes);
    mfem::GridFunction streamfunction_gf(&scalar_fes);
    mfem::GridFunction velocity_gf(&vector_fes);

    // Create ParaView data collection
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

// ====================================================================
// TIME-STEPPING PARAMETERS
// ====================================================================
    // Todo:  check ComputeStableTimestep
    // params.ComputeStableTimestep(num_nodes, reynolds_number);
    double current_dt = params.dt;
    int num_timesteps = static_cast<int>(params.final_time / params.dt);

    // Initialize steady state monitoring
    mfem::Vector steady_prev_vorticity(num_nodes);
    mfem::Vector steady_prev_streamfunction(num_nodes);
    int steady_consecutive_passes = 0;
    bool steady_initialized = false;

    ResidualHistory residual_history;  // Always create (empty if not used)
    if (params.check_residuals) {
        std::cout << "\nResidual monitoring ENABLED" << std::endl;
        std::cout << "  Check frequency: every " << params.residual_check_freq
                << " iterations" << std::endl;
        std::cout << "  Output file: " << dat_dir << "/residual_history.csv" << std::endl;
    } else {
        std::cout << "\nResidual monitoring DISABLED" << std::endl;
    }

    std::cout << "Starting simulation: " << num_nodes << " nodes, "
              << num_timesteps << " timesteps\n\n";


// ====================================================================
// TIME-STEPPING LOOP
// ====================================================================
    // ================== Timing PHASE 3: SIMULATION  ==================
    mfem::StopWatch main_loop_timer;
    main_loop_timer.Start();

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
        perf_metrics.total_iterations = time_step;

        vorticity_old = vorticity;

        // ================================================================
        // COMPUTE ALL DERIVATIVES ONCE (reused by all subsequent steps)
        // ================================================================
        dy_matrix.Mult(streamfunction, dpsi_dy);      // ∂ψ/∂y = u
        dx_matrix.Mult(streamfunction, dpsi_dx);      // ∂ψ/∂x = -v (will negate later)
        dx_matrix.Mult(vorticity, domega_dx);         // ∂ω/∂x
        dy_matrix.Mult(vorticity, domega_dy);         // ∂ω/∂y
        dxx_matrix.Mult(vorticity, d2omega_dx2);      // ∂²ω/∂x²
        dyy_matrix.Mult(vorticity, d2omega_dy2);      // ∂²ω/∂y²

        // ================================================================
        // STEP 1: UPDATE VORTICITY (using derivatives computed above)
        // ================================================================
        // Use explicit Euler scheme (Eq. 11)
        for (int i = 0; i < num_nodes; ++i) {
            double convection = dpsi_dy[i] * domega_dx[i] - dpsi_dx[i] * domega_dy[i];
            double diffusion = (1.0 / params.reynolds_number) * (d2omega_dx2[i] + d2omega_dy2[i]);
            vorticity[i] += current_dt * (diffusion - convection);
        }

        // ================================================================
        // STEP 2: SOLVE STREAMFUNCTION POISSON EQUATION (Equation 12)
        // ================================================================
        // Solve -∇²ψ = ω instead ∇²ψ = -ω, so that RHS stays positive for linear solvers (cg)
        rhs = vorticity;

        // Apply homogeneous Dirichlet boundary conditions: ψ = 0 on all boundaries
        for (int idx : all_boundary_nodes) {
            rhs[idx] = 0.0;
        }

        linear_solver->Mult(rhs, streamfunction);

        // ================================================================
        // STEP 3: APPLY BOUNDARY CONDITIONS
        // ================================================================
        // We already have u = dpsi_dy and need v = -dpsi_dx
        // Apply velocity BCs first
        // Lid driven cavity
        for (int idx : bottom_nodes) { dpsi_dy[idx] = 0.0; dpsi_dx[idx] = 0.0; }  // u=0, v=0
        for (int idx : right_nodes)  { dpsi_dy[idx] = 0.0; dpsi_dx[idx] = 0.0; }  // u=0, v=0
        for (int idx : left_nodes)   { dpsi_dy[idx] = 0.0; dpsi_dx[idx] = 0.0; }  // u=0, v=0
        for (int idx : top_nodes)    { dpsi_dy[idx] = 1.0; dpsi_dx[idx] = 0.0; }  // u=1, v=0 (lid)

        // Compute boundary vorticity: ω = ∂v/∂x - ∂u/∂y
        // We need derivatives of velocity (which are already enforced at boundaries)
        mfem::Vector du_dy(num_nodes);
        mfem::Vector dv_dx(num_nodes);
        dy_matrix.Mult(dpsi_dy, du_dy);  // ∂u/∂y (using BC-enforced u)
        dx_matrix.Mult(dpsi_dx, dv_dx);  // ∂(-v)/∂x = -∂v/∂x
        dv_dx *= -1.0;                   // Now it's ∂v/∂x

        // Apply to all boundary nodes
        for (int idx : all_boundary_nodes) {
            vorticity[idx] = dv_dx[idx] - du_dy[idx];
        }

        // ================================================================
        // STEP 4: COMPUTE AND LOG RESIDUALS (if enabled)
        // ================================================================
        if (params.check_residuals && (time_step % params.residual_check_freq == 0)) {

            // Vorticity residual: ∂ω/∂t + K(ψ,ω) - L(ω) = 0
            double vort_sum_sq = 0.0, vort_max = 0.0;
            for (int idx : interior_nodes) {
                double domega_dt = (vorticity[idx] - vorticity_old[idx]) / current_dt;
                double conv = dpsi_dy[idx] * domega_dx[idx] - dpsi_dx[idx] * domega_dy[idx];
                double diff = (1.0 / params.reynolds_number) *
                            (d2omega_dx2[idx] + d2omega_dy2[idx]);
                double res = domega_dt + conv - diff;
                vort_sum_sq += res * res;
                vort_max = std::max(vort_max, std::abs(res));
            }
            double vort_rms = std::sqrt(vort_sum_sq / interior_nodes.size());

            // Streamfunction residual: ∇²ψ + ω = 0
            mfem::Vector d2psi_dx2(num_nodes), d2psi_dy2(num_nodes);
            dxx_matrix.Mult(streamfunction, d2psi_dx2);
            dyy_matrix.Mult(streamfunction, d2psi_dy2);

            double psi_sum_sq = 0.0, psi_max = 0.0;
            for (int idx : interior_nodes) {
                double res = d2psi_dx2[idx] + d2psi_dy2[idx] + vorticity[idx];
                psi_sum_sq += res * res;
                psi_max = std::max(psi_max, std::abs(res));
            }
            double psi_rms = std::sqrt(psi_sum_sq / interior_nodes.size());

            // Compute relative changes for steady state context
            double vort_rel_change = 0.0;
            double psi_rel_change = 0.0;
            if (steady_initialized) {
                mfem::Vector diff_vort(num_nodes);
                mfem::Vector diff_stream(num_nodes);
                diff_vort = vorticity;
                diff_vort -= steady_prev_vorticity;
                diff_stream = streamfunction;
                diff_stream -= steady_prev_streamfunction;
                vort_rel_change = diff_vort.Norml2() / std::max(vorticity.Norml2(), 1e-12);
                psi_rel_change = diff_stream.Norml2() / std::max(streamfunction.Norml2(), 1e-12);
            }

            // Store in history
            residual_history.Add(time_step, current_time,
                                vort_rms, vort_max, psi_rms, psi_max,
                                vort_rel_change, psi_rel_change);

            // Print periodic updates
            // if (time_step % (params.residual_check_freq * 10) == 0) {
            //     std::cout << "Step " << time_step
            //             << " (t=" << std::fixed << std::setprecision(2) << current_time << "): "
            //             << "ω_RMS=" << std::scientific << std::setprecision(2) << vort_rms
            //             << ", ψ_RMS=" << psi_rms;
            //     if (steady_initialized) {
            //         std::cout << ", Δω=" << vort_rel_change
            //                 << ", Δψ=" << psi_rel_change;
            //     }
            //     std::cout << std::endl;
            // }
        }

        // ================================================================
        // STEP 5: ADAPTIVE TIMESTEP (if enabled)
        // ================================================================
        if (params.enable_adaptive_timestep &&
            time_step % params.gershgorin_frequency == 0) {
            // Gershgorin critical timestep (inline to avoid recomputing derivatives)
            double max_row_sum = 0.0;

            for (int i = 0; i < num_nodes; ++i) {
                // Convection operator contribution
                double K_sum = std::abs(dpsi_dx[i]) + std::abs(dpsi_dy[i]);

                // Diffusion operator contribution
                double L_sum = (1.0 / params.reynolds_number) *
                              (std::abs(d2omega_dx2[i]) + std::abs(d2omega_dy2[i]));

                max_row_sum = std::max(max_row_sum, K_sum + L_sum);
            }

            double dt_critical = 2.0 / max_row_sum;

            if (current_dt > dt_critical / params.gershgorin_safety_factor) {
                current_dt = dt_critical / params.gershgorin_safety_factor;
                num_timesteps = static_cast<int>(params.final_time / current_dt);
                std::cout << "Adjusted dt = " << current_dt << std::endl;
            }
        }

        // ================================================================
        // STEP 6: OUTPUT SOLUTION PERIODICALLY (if enabled)
        // ================================================================
        if (save_solutions && (time_step % params.output_frequency == 0)) {
            // dpsi_dy = u-velocity
            // Create v-velocity from dpsi_dx with sign correction
            mfem::Vector v_velocity = dpsi_dx;
            v_velocity *= -1.0;

            SaveSolutionToFile(vorticity, streamfunction, dpsi_dy, v_velocity,
                            params.output_prefix, time_step, dat_dir, false);
        }

        // ================================================================
        // STEP 7: PARAVIEW OUTPUT (if enabled)
        // ================================================================
        if (paraview_output && (time_step % params.paraview_output_frequency == 0)) {
            // dpsi_dy is already u-velocity
            // dpsi_dx needs to be negated for v-velocity

            for (int i = 0; i < num_nodes; ++i) {
                vorticity_gf[i] = vorticity[i];
                streamfunction_gf[i] = streamfunction[i];
                velocity_gf[i] = dpsi_dy[i];                    // u
                velocity_gf[i + num_nodes] = -dpsi_dx[i];      // v = -∂ψ/∂x
            }

            paraview_dc.SetTime(current_time);
            paraview_dc.SetCycle(time_step);
            paraview_dc.Save();
        }

        // ================================================================
        // STEP 8: CHECK STEADY STATE
        // ================================================================
        if (time_step >= 100 && time_step % params.steady_state_freq == 0) {

            // Initialize on first check
            if (!steady_initialized) {
                steady_prev_vorticity = vorticity;
                steady_prev_streamfunction = streamfunction;
                steady_initialized = true;
                std::cout << "Steady state monitoring started at timestep "
                        << time_step << std::endl;
            } else {
                // Compute relative changes since last check
                mfem::Vector diff_vort(num_nodes);
                mfem::Vector diff_stream(num_nodes);

                diff_vort = vorticity;
                diff_vort -= steady_prev_vorticity;

                diff_stream = streamfunction;
                diff_stream -= steady_prev_streamfunction;

                double vort_change = diff_vort.Norml2() /
                                    std::max(vorticity.Norml2(), 1e-12);
                double stream_change = diff_stream.Norml2() /
                                    std::max(streamfunction.Norml2(), 1e-12);
                double max_change = std::max(vort_change, stream_change);

                // Report progress
                std::cout << "Step " << time_step << " (t=" << std::fixed
                        << std::setprecision(6) << current_time
                        << "): ‖Δ‖/‖·‖ = " << std::scientific << max_change;

                // Show current residuals if available
                if (params.check_residuals && !residual_history.timesteps.empty()) {
                    size_t last_idx = residual_history.timesteps.size() - 1;
                    std::cout << " [ω_RMS=" << residual_history.vorticity_rms[last_idx]
                            << ", ψ_RMS=" << residual_history.streamfunction_rms[last_idx] << "]";
                }

                // Check convergence
                if (max_change < params.steady_state_tol) {
                    steady_consecutive_passes++;
                    std::cout << " ✓ [" << steady_consecutive_passes << "/"
                            << params.steady_state_checks << "]" << std::endl;

                    if (steady_consecutive_passes >= params.steady_state_checks) {
                        std::cout << "\n" << std::string(70, '=') << "\n";
                        std::cout << "STEADY STATE REACHED\n";
                        std::cout << "Time: " << GetCurrentDateTime() << "\n";
                        std::cout << std::string(70, '=') << std::endl;
                        std::cout << "Final relative change: " << max_change << std::endl;
                        std::cout << "Convergence tolerance: " << params.steady_state_tol << std::endl;
                        std::cout << "Total timesteps: " << time_step << std::endl;
                        std::cout << "Physical time: " << current_time << std::endl;

                        if (params.check_residuals && !residual_history.timesteps.empty()) {
                            size_t last_idx = residual_history.timesteps.size() - 1;
                            std::cout << "Final vorticity residual: "
                                    << residual_history.vorticity_rms[last_idx] << std::endl;
                            std::cout << "Final streamfunction residual: "
                                    << residual_history.streamfunction_rms[last_idx] << std::endl;
                        }

                        std::cout << std::string(70, '=') << "\n" << std::endl;

                        // Save final solution before breaking
                        if (save_solutions) {
                            mfem::Vector v_final = dpsi_dx;
                            v_final *= -1.0;
                            SaveSolutionToFile(vorticity, streamfunction,
                                            dpsi_dy, v_final,
                                            params.output_prefix, time_step,
                                            dat_dir, true);
                        }

                        break;  // Exit time-stepping loop
                    }
                } else {
                    std::cout << " ✗" << std::endl;
                    steady_consecutive_passes = 0; // Reset counter
                }

                // Update previous values for next check
                steady_prev_vorticity = vorticity;
                steady_prev_streamfunction = streamfunction;
            }
        }
    } // End of time-stepping loop

    main_loop_timer.Stop();
    double total_solve_time = main_loop_timer.RealTime();

    std::cout << "\nSimulation completed successfully in "
            << main_loop_timer.RealTime() << " seconds." << std::endl;

// ====================================================================
// POST-PROCESSING
// ====================================================================
    // Save residual history if enabled and data exists
    if (params.check_residuals && !residual_history.timesteps.empty()) {
        std::string residual_file = dat_dir + "/residual_history.csv";
        residual_history.SaveToCSV(residual_file);
        residual_history.PrintSummary();
    }

    // Print performance table
    perf_metrics.solution_time_per_iter = total_solve_time / perf_metrics.total_iterations;
    perf_metrics.PrintPerformanceTable(perf_metrics);

    // Additional statistics
    std::cout << "Additional Performance Statistics:\n";
    std::cout << "  Total simulation time: " << total_solve_time << " s\n";
    std::cout << "  Total iterations completed: " << perf_metrics.total_iterations << "\n";
    std::cout << "  Average time per iteration: " << total_solve_time / perf_metrics.total_iterations << " s\n\n";


    // Final solution statistics
    std::cout << "Final simulation statistics:" << std::endl;
    std::cout << "  - Final time step size: " << current_dt << std::endl;
    std::cout << "  - Total time steps completed: " << std::min(num_timesteps, static_cast<int>(params.final_time / current_dt)) << std::endl;
    std::cout << "  - Maximum vorticity: " << vorticity.Max() << std::endl;
    std::cout << "  - Minimum vorticity: " << vorticity.Min() << std::endl;
    std::cout << "  - Maximum streamfunction: " << streamfunction.Max() << std::endl;
    std::cout << "  - Minimum streamfunction: " << streamfunction.Min() << std::endl;

    // Cleanup
    delete derivs;
    delete mesh;
    delete linear_solver;
    delete laplacian_matrix;

    std::cout << "main: success!" << std::endl;
    return EXIT_SUCCESS;
}


mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh) {
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


StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int NumNeighbors) {
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
                           int dim, bool save_d, bool save_dd, std::string dat_dir) {
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


// Functions matching MATLAB code "Explicit_streamfunction_vorticity_meshless.m"
void IdentifyBoundaryNodes(mfem::Mesh* mesh,
                                          std::vector<int>& bottom_nodes,
                                          std::vector<int>& right_nodes,
                                          std::vector<int>& top_nodes,
                                          std::vector<int>& left_nodes,
                                          std::vector<int>& interior_nodes) {
    const double tolerance = 1e-10;
    const int num_vertices = mesh->GetNV();

    // Clear all vectors
    bottom_nodes.clear();
    right_nodes.clear();
    top_nodes.clear();
    left_nodes.clear();
    interior_nodes.clear();

    // First pass: collect all boundary nodes
    std::vector<int> all_boundary_nodes;

    for (int i = 0; i < num_vertices; ++i) {
        const double* vertex = mesh->GetVertex(i);
        double x = vertex[0];
        double y = vertex[1];

        // Bottom boundary: y == 0 (includes corners (0,0) and (1,0))
        if (std::abs(y) < tolerance) {
            bottom_nodes.push_back(i);
            all_boundary_nodes.push_back(i);
        }
        // Right boundary: x == 1 AND 0 < y < 1 (excludes corners)
        else if (std::abs(x - 1.0) < tolerance &&
                 y > tolerance &&
                 y < (1.0 - tolerance)) {
            right_nodes.push_back(i);
            all_boundary_nodes.push_back(i);
        }
        // Top boundary: y == 1 (includes corners (0,1) and (1,1))
        else if (std::abs(y - 1.0) < tolerance) {
            top_nodes.push_back(i);
            all_boundary_nodes.push_back(i);
        }
        // Left boundary: x == 0 AND 0 < y < 1 (excludes corners)
        else if (std::abs(x) < tolerance &&
                 y > tolerance &&
                 y < (1.0 - tolerance)) {
            left_nodes.push_back(i);
            all_boundary_nodes.push_back(i);
        }
        // Interior nodes: everything else
        else {
            interior_nodes.push_back(i);
        }
    }

    // Create combined boundary nodes vector (nb in MATLAB)
    // This matches: nb=[bottom; right; top; left];

    // Diagnostic output
    std::cout << "IdentifyBoundaryNodes:" << std::endl;
    std::cout << "  Bottom: " << bottom_nodes.size() << " nodes" << std::endl;
    std::cout << "  Right:  " << right_nodes.size() << " nodes" << std::endl;
    std::cout << "  Top:    " << top_nodes.size() << " nodes" << std::endl;
    std::cout << "  Left:   " << left_nodes.size() << " nodes" << std::endl;
    std::cout << "  Total boundary: " << all_boundary_nodes.size() << " nodes" << std::endl;
    std::cout << "  Interior: " << interior_nodes.size() << " nodes" << std::endl;

    // Expected counts for an N×N grid:
    // - Bottom and Top: N nodes each (including corners)
    // - Right and Left: (N-2) nodes each (excluding corners)
    // - Total boundary: 4N - 4
    // - Interior: (N-2)²

    // For a 20×20 grid: 20, 18, 20, 18 (total boundary: 76)
    // For a 40×40 grid: 40, 38, 40, 38 (total boundary: 156)
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
    // v *= -1.0;

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


void SaveSolutionToFile(const mfem::Vector& vorticity,
                       const mfem::Vector& streamfunction,
                       const mfem::Vector& u_velocity,
                       const mfem::Vector& v_velocity,
                       const std::string& filename,
                       int timestep,
                       std::string dat_dir,
                       bool is_final = false)
{
    // Determine suffix based on whether this is final or intermediate save
    std::string suffix = is_final ? "_final" : "_solution";

    // Save vorticity
    std::string vort_filename = dat_dir + "/" + filename + "_vorticity" + suffix + ".dat";
    std::ofstream vort_file(vort_filename);
    if (vort_file.is_open()) {
        vort_file.precision(12);
        for (int i = 0; i < vorticity.Size(); ++i) {
            vort_file << vorticity[i] << std::endl;
        }
        vort_file.close();
    }

    // Save streamfunction
    std::string stream_filename = dat_dir + "/" + filename + "_streamfunction" + suffix + ".dat";
    std::ofstream stream_file(stream_filename);
    if (stream_file.is_open()) {
        stream_file.precision(12);
        for (int i = 0; i < streamfunction.Size(); ++i) {
            stream_file << streamfunction[i] << std::endl;
        }
        stream_file.close();
    }

    // Save u-velocity
    std::string u_filename = dat_dir + "/" + filename + "_u_velocity" + suffix + ".dat";
    std::ofstream u_file(u_filename);
    if (u_file.is_open()) {
        u_file.precision(12);
        for (int i = 0; i < u_velocity.Size(); ++i) {
            u_file << u_velocity[i] << std::endl;
        }
        u_file.close();
    }

    // Save v-velocity
    std::string v_filename = dat_dir + "/" + filename + "_v_velocity" + suffix + ".dat";
    std::ofstream v_file(v_filename);
    if (v_file.is_open()) {
        v_file.precision(12);
        for (int i = 0; i < v_velocity.Size(); ++i) {
            v_file << v_velocity[i] << std::endl;
        }
        v_file.close();
    }

    if (is_final) {
        std::cout << "SaveSolutionToFile: Saved final solution at timestep " << timestep << std::endl;
        std::cout << "  Files saved to: " << dat_dir << std::endl;
    }
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
        // Option 1: Gauss-Seidel (simple, effective)
        mfem::GSSmoother* gs_precond = new mfem::GSSmoother();
        gmres_solver->SetPreconditioner(*gs_precond);

        // Option 2: Jacobi (parallel-friendly)
        // mfem::DSmoother* jacobi_precond = new mfem::DSmoother();
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
    else if (solver_type == "hypre") { // Todo: HypreBoomerAMG does not work (require HypreParMatrix)
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
    else {
        std::cerr << "ERROR: Unknown solver type: " << solver_type << std::endl;
        std::cerr << "Available solvers: umfpack, klu, cg, gmres, minres, bicgstab, hypre, amg" << std::endl;
        return nullptr;
    }

    return solver;
}

std::string GetCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

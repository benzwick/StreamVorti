/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 * Copyright (C) 2026 Weizheng (Will) Li
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
 *      Weizheng (Will) LI
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

// PARALLEL VERSION - requires MFEM built with MPI and HYPRE
//
// Demo usage: (save everything, 4 MPI processes)
//     mpirun -np 4 ./StreamVorti_par -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -nn 25 -Re 1000 -sm -sn -sd -sdd -ss -pv -cr
//
// SDL usage (requires ECL):
//     mpirun -np 4 ./StreamVorti_par -f demo/cavity.lisp -pv

// Related header
#include <StreamVorti/streamvorti_par.hpp>

// SDL/Lisp support (conditional)
#ifdef STREAMVORTI_WITH_ECL
#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_loader.hpp>
#endif

// C++ standard library headers (alphabetically sorted)
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Other libraries' headers
#include "mfem.hpp"

// MPI header (required for parallel version)
#include <mpi.h>

// Conditional includes
#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// FUNCTION IMPLEMENTATIONS
// ============================================================================
// All struct definitions and function declarations have been moved to
// include/StreamVorti/stream_vorti.hpp for proper header/implementation split.

// ============================================================================
// PARALLEL-SPECIFIC BOUNDARY NODE HELPERS
// ============================================================================

#ifdef MFEM_USE_MPI

/**
 * @brief Parallel version of IdentifyBoundaryNodes using proper DOF ownership.
 *
 * Unlike the serial version, this iterates over ALL local vertices (GetNV())
 * and uses GetLocalTDofNumber() to skip non-owned vertices. Node arrays are
 * filled with TRUE DOF indices (0..TrueVSize-1), matching solution vector
 * indexing.
 *
 * In MFEM's ParMesh the owned vertices are NOT guaranteed to be the first
 * TrueVSize entries; shared vertices can appear at any position in the local
 * array.  The previous approach of iterating 0..TrueVSize-1 missed boundary
 * nodes whose local index happened to exceed TrueVSize.
 */
static void IdentifyBoundaryNodesPar(
    mfem::ParMesh* pmesh,
    mfem::ParFiniteElementSpace& pfes,
    std::vector<int>& bottom_nodes,
    std::vector<int>& right_nodes,
    std::vector<int>& top_nodes,
    std::vector<int>& left_nodes,
    std::vector<int>& interior_nodes)
{
    const double tol = 1e-10;
    bottom_nodes.clear(); right_nodes.clear(); top_nodes.clear();
    left_nodes.clear();   interior_nodes.clear();

    for (int i = 0; i < pmesh->GetNV(); ++i) {
        int tdof = pfes.GetLocalTDofNumber(i);
        if (tdof < 0) continue;  // Not owned by this rank

        const double* v = pmesh->GetVertex(i);
        double x = v[0], y = v[1];

        if (std::abs(y) < tol) {
            bottom_nodes.push_back(tdof);
        } else if (std::abs(x - 1.0) < tol && y > tol && y < (1.0 - tol)) {
            right_nodes.push_back(tdof);
        } else if (std::abs(y - 1.0) < tol) {
            top_nodes.push_back(tdof);
        } else if (std::abs(x) < tol && y > tol && y < (1.0 - tol)) {
            left_nodes.push_back(tdof);
        } else {
            interior_nodes.push_back(tdof);
        }
    }

    std::cout << "IdentifyBoundaryNodesPar (rank " << pfes.GetMyRank() << "):"
              << " Bottom=" << bottom_nodes.size()
              << " Right=" << right_nodes.size()
              << " Top=" << top_nodes.size()
              << " Left=" << left_nodes.size()
              << " Interior=" << interior_nodes.size() << std::endl;
}

/**
 * @brief Parallel version of IdentifyBoundaryNodesByAttribute using proper DOF ownership.
 *
 * Same ownership fix as IdentifyBoundaryNodesPar. Node indices stored are
 * TRUE DOF indices.
 */
static void IdentifyBoundaryNodesByAttributePar(
    mfem::ParMesh* pmesh,
    mfem::ParFiniteElementSpace& pfes,
    std::map<int, std::vector<int>>& boundary_nodes,
    std::vector<int>& interior_nodes)
{
    boundary_nodes.clear();
    interior_nodes.clear();

    int nv = pmesh->GetNV();
    std::vector<bool> is_boundary(nv, false);

    // Collect vertices from each boundary element, grouped by attribute.
    // MFEM assigns attributes to boundary *elements* (edges in 2D), not
    // vertices. At corners, two edges with different attributes share a
    // vertex, so that vertex appears in both attribute lists. This is
    // intentional — the BC application loop processes BCs in order and
    // last-writer-wins at shared corners (standard FEM practice).
    // The std::find check only prevents duplicates within a single attribute.
    for (int be = 0; be < pmesh->GetNBE(); ++be) {
        int attr = pmesh->GetBdrAttribute(be);
        mfem::Array<int> vertices;
        pmesh->GetBdrElementVertices(be, vertices);

        for (int v = 0; v < vertices.Size(); ++v) {
            int vi = vertices[v];
            int tdof = pfes.GetLocalTDofNumber(vi);
            if (tdof < 0) continue;  // Not owned by this rank

            is_boundary[vi] = true;
            auto& nodes = boundary_nodes[attr];
            if (std::find(nodes.begin(), nodes.end(), tdof) == nodes.end()) {
                nodes.push_back(tdof);
            }
        }
    }

    // Collect interior owned vertices
    for (int i = 0; i < nv; ++i) {
        int tdof = pfes.GetLocalTDofNumber(i);
        if (tdof < 0) continue;
        if (!is_boundary[i]) {
            interior_nodes.push_back(tdof);
        }
    }
}

#endif  // MFEM_USE_MPI

// ============================================================================
// MAIN FUNCTION
// ============================================================================

/**
 * @brief PARALLEL version: 2D lid-driven cavity flow solver with MPI/Hypre
 *
 * Solves incompressible Navier-Stokes equations using stream function-vorticity
 * formulation with Meshless Point Collocation (MPC) and DC PSE operators.
 * Distributed across multiple MPI processes using MFEM's parallel data structures.
 *
 * Algorithm (from Bourantas et al. 2019):
 * 1. Compute DC PSE derivative operators (pre-processing)
 * 2. Time-stepping loop:
 *    a) Update vorticity: explicit Euler for transport equation
 *    b) Solve streamfunction: ∇²ψ = -ω (Poisson equation with Hypre)
 *    c) Apply boundary conditions for lid-driven cavity
 *    d) Check steady-state convergence
 *    e) Output solutions (DAT files, ParaView with parallel I/O)
 *
 * @param argc Number of command-line arguments
 * @param argv Command-line argument values
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int main(int argc, char *argv[])
{
#ifndef MFEM_USE_MPI
    #error "StreamVorti_par requires MFEM built with MPI support (MFEM_USE_MPI)"
#endif

    // Initialize MPI (required for parallel execution)
    MPI_Init(&argc, &argv);
    int num_procs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0)
    {
        std::cout << "MPI initialized with " << num_procs << " process(es)" << std::endl;
    }

    // Options
    const char *mesh_file = "";
    const char *sdl_file = "";  // SDL file for Lisp-based configuration
    const char *lisp_path = ""; // Path to SDL Lisp source files
    int order = 1;
    const char *method = "dcpse"; // Spatial discretization method: "dcpse" or "fdm"

    // Output filename prefix and extension
    std::string fname = "mfem_square10x10";
    std::string fext = ".dat";

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
    params.solver_type = "gmres";  // DCPSE Laplacian is non-symmetric; GMRES required (CG stalls for some partitions)
    // Performance Metrics parameters
    PerformanceMetrics perf_metrics;

    // Parse command-line options
    mfem::OptionsParser args(argc, argv);
#ifdef STREAMVORTI_WITH_ECL
    args.AddOption(&sdl_file, "-f", "--sdl-file",
                   "SDL (Simulation Definition Language) file to load.");
    args.AddOption(&lisp_path, "-lp", "--lisp-path",
                   "Path to SDL Lisp source files.");
#endif
    // Simulation options
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
    args.AddOption(&method, "-method", "--spatial-method",
                   "Spatial discretization method: dcpse (default: dcpse). FDM not yet supported in parallel.");
    args.AddOption(&params.num_neighbors, "-nn", "--num-neighbors",
                    "Number of neighbors for DCPSE (default: 25).");
    args.AddOption(&params.reynolds_number, "-Re", "--reynolds-number",
                    "Reynolds number of simulating fluid (default: 1000).");
    args.AddOption(&params.dt, "-dt", "--timestep",
                    "Time step size (default: 1e-4).");
    args.AddOption(&params.final_time, "-ft", "--final-time",
                    "Final simulation time (default: 60000).");
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
        if (myid == 0) {
            args.PrintUsage(std::cout);
        }
        return 1;
    }

    // Check for unsupported parallel methods
    std::string spatial_method(method);
    if (spatial_method == "fdm") {
        if (myid == 0) {
            std::cerr << "Error: FDM (finite differences) is not yet supported in parallel.\n"
                      << "Use the serial StreamVorti executable with -method fdm, or use -method dcpse.\n";
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (myid == 0) {
        args.PrintOptions(std::cout);
    }

    // Mesh pointer - serial mesh for initial loading, then converted to ParMesh
    mfem::Mesh* serial_mesh = nullptr;
#ifdef STREAMVORTI_WITH_ECL
    bool sdl_mode = false;
    // Store SDL boundary conditions for runtime evaluation
    std::vector<StreamVorti::Lisp::BoundaryCondition> sdl_boundaries;
    std::vector<StreamVorti::Lisp::LineProbe> sdl_probes;
    bool use_sdl_bcs = false;
#endif

#ifdef STREAMVORTI_WITH_ECL
    // ===== SDL Mode =====
    // If an SDL file is specified, use the Lisp-based configuration
    if (sdl_file[0] != '\0')
    {
        sdl_mode = true;
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "SDL MODE - Loading simulation from: " << sdl_file << std::endl;
        std::cout << std::string(70, '=') << "\n" << std::endl;

        try {
            // Initialize ECL runtime
            std::string lpath = lisp_path[0] != '\0' ? lisp_path : "lisp";
            StreamVorti::Lisp::Runtime::init(lpath);

            // Load simulation from SDL file
            auto config = StreamVorti::Lisp::Loader::load(sdl_file);

            std::cout << "Loaded simulation: " << config.name << std::endl;
            std::cout << "  Dimension: " << config.dimension << std::endl;
            std::cout << "  Mesh vertices: " << config.mesh->GetNV() << std::endl;
            std::cout << "  Mesh elements: " << config.mesh->GetNE() << std::endl;
            std::cout << "  Boundaries defined: " << config.boundaries.size() << std::endl;

            // Store SDL boundary conditions for runtime evaluation
            if (!config.boundaries.empty()) {
                sdl_boundaries = std::move(config.boundaries);
                use_sdl_bcs = true;
                std::cout << "  Using SDL boundary conditions (" << sdl_boundaries.size() << " regions):\n";
                for (const auto& bc : sdl_boundaries) {
                    std::cout << "    - " << bc.name << " (attr=" << bc.attribute
                              << ", type=" << bc.type;
                    if (bc.u_function) std::cout << ", u-func";
                    if (bc.v_function) std::cout << ", v-func";
                    if (bc.function && !bc.u_function) std::cout << ", scalar-func";
                    std::cout << ")\n";
                }
            } else {
                std::cout << "  No SDL boundary conditions defined, using hardcoded lid-driven cavity BCs.\n";
            }

            // Store SDL line probes for post-processing
            if (!config.probes.empty()) {
                sdl_probes = std::move(config.probes);
                std::cout << "  Line probes (" << sdl_probes.size() << "):\n";
                for (const auto& p : sdl_probes) {
                    std::cout << "    - " << p.name << " (" << p.axis << " = " << p.position << ")\n";
                }
            }

            // Map SDL config to SimulationParams
            params.num_neighbors = config.dcpse.num_neighbors;
            params.reynolds_number = static_cast<int>(config.physics.reynolds);
            params.dt = config.solver.dt;
            params.final_time = config.solver.end_time;
            params.output_prefix = config.name;

            // Map solver tolerance and iterations
            params.solver_rel_tol = config.solver.tolerance;
            params.solver_max_iter = config.solver.max_iterations;

            // Map output settings
            if (config.output.format == "vtk") {
                paraview_output = true;
            }
            // Map output interval to paraview frequency (convert time to steps)
            if (config.output.interval > 0 && params.dt > 0) {
                params.paraview_output_frequency = static_cast<int>(config.output.interval / params.dt);
            }

            std::cout << "\nMapped SDL parameters:" << std::endl;
            std::cout << "  Reynolds number: " << params.reynolds_number << std::endl;
            std::cout << "  Timestep: " << params.dt << std::endl;
            std::cout << "  Final time: " << params.final_time << std::endl;
            std::cout << "  DCPSE neighbors: " << params.num_neighbors << std::endl;
            std::cout << "  Solver tolerance: " << params.solver_rel_tol << std::endl;
            std::cout << "  ParaView output frequency: " << params.paraview_output_frequency << std::endl;

            // Get mesh dimensions from SDL config (for grid_size reporting)
            // Note: GetNV() = (nx+1)*(ny+1) for a structured mesh, so sqrt gives nx+1
            int nv = config.mesh->GetNV();
            nx = ny = static_cast<int>(std::sqrt(static_cast<double>(nv))) - 1;
            dim = config.dimension;

            // Transfer ownership of mesh (serial for now, will convert to ParMesh)
            serial_mesh = config.mesh.release();

            std::cout << "\nContinuing with SDL-configured simulation...\n" << std::endl;

        } catch (const StreamVorti::Lisp::EclException& e) {
            if (myid == 0) {
                std::cerr << "ECL Error: " << e.what() << std::endl;
            }
            StreamVorti::Lisp::Runtime::shutdown();
            return EXIT_FAILURE;
        } catch (const std::exception& e) {
            if (myid == 0) {
                std::cerr << "Error loading SDL file: " << e.what() << std::endl;
            }
            StreamVorti::Lisp::Runtime::shutdown();
            return EXIT_FAILURE;
        }
    }
#endif // STREAMVORTI_WITH_ECL

    // Create or load mesh (traditional mode - only if not in SDL mode)
    if (serial_mesh == nullptr) {
        serial_mesh = CreateOrLoadMesh(mesh_file, dim, nx, ny, nz, sx, sy, sz, save_mesh);
    }
    perf_metrics.grid_size = nx;

#ifdef STREAMVORTI_WITH_ECL
    // ================================================================
    // Compute wall ψ constants using serial mesh boundary connectivity.
    // This must happen BEFORE delete serial_mesh, while the global mesh
    // is still available (same algorithm as serial streamvorti.cpp).
    // ================================================================
    std::map<int, double> wall_psi_by_attr;
    if (use_sdl_bcs) {
        // Build vertex → set of boundary attributes map
        std::map<int, std::set<int>> vertex_attrs;
        for (int be = 0; be < serial_mesh->GetNBE(); ++be) {
            int attr = serial_mesh->GetBdrAttribute(be);
            mfem::Array<int> verts;
            serial_mesh->GetBdrElementVertices(be, verts);
            for (int j = 0; j < verts.Size(); ++j) {
                vertex_attrs[verts[j]].insert(attr);
            }
        }

        // Build attribute → BC map
        std::map<int, const StreamVorti::Lisp::BoundaryCondition*> attr_to_bc;
        for (const auto& bc : sdl_boundaries) {
            attr_to_bc[bc.attribute] = &bc;
        }

        // For each wall BC, find adjacent velocity BC via corner vertices
        constexpr int nsub_wall = 20;
        for (const auto& wall_bc : sdl_boundaries) {
            if (wall_bc.type != "no-slip" && wall_bc.type != "slip") continue;

            double wall_psi = 0.0;
            bool found = false;
            for (const auto& [vi, attrs] : vertex_attrs) {
                if (attrs.find(wall_bc.attribute) == attrs.end()) continue;
                if (attrs.size() < 2) continue;  // not a corner

                for (int other_attr : attrs) {
                    if (other_attr == wall_bc.attribute) continue;
                    auto bc_it = attr_to_bc.find(other_attr);
                    if (bc_it == attr_to_bc.end()) continue;
                    const auto* vel_bc = bc_it->second;
                    if (vel_bc->type != "velocity") continue;

                    const double* corner = serial_mesh->GetVertex(vi);
                    double psi_val = 0.0;
                    if (vel_bc->predicate_axis == 'x' && vel_bc->u_function) {
                        double x0 = vel_bc->predicate_value;
                        double y = corner[1];
                        double h = y / nsub_wall;
                        for (int k = 0; k < nsub_wall; ++k) {
                            double a = k * h, b = (k + 1) * h, m = 0.5 * (a + b);
                            psi_val += (h / 6.0) * (
                                vel_bc->u_function->evaluateAt(x0, a, 0.0) +
                                4.0 * vel_bc->u_function->evaluateAt(x0, m, 0.0) +
                                vel_bc->u_function->evaluateAt(x0, b, 0.0));
                        }
                    } else if (vel_bc->predicate_axis == 'y' && vel_bc->v_function) {
                        double y0 = vel_bc->predicate_value;
                        double x = corner[0];
                        double h = x / nsub_wall;
                        for (int k = 0; k < nsub_wall; ++k) {
                            double a = k * h, b = (k + 1) * h, m = 0.5 * (a + b);
                            psi_val -= (h / 6.0) * (
                                vel_bc->v_function->evaluateAt(a, y0, 0.0) +
                                4.0 * vel_bc->v_function->evaluateAt(m, y0, 0.0) +
                                vel_bc->v_function->evaluateAt(b, y0, 0.0));
                        }
                    }
                    wall_psi = psi_val;
                    found = true;
                    break;
                }
                if (found) break;
            }

            if (found) {
                wall_psi_by_attr[wall_bc.attribute] = wall_psi;
            }
        }
    }
#endif

    // Convert serial mesh to parallel mesh (distributed across MPI ranks)
    mfem::ParMesh* mesh = new mfem::ParMesh(MPI_COMM_WORLD, *serial_mesh);
    delete serial_mesh;  // No longer needed after creating ParMesh

    // dat files for Matlab
    std::string dat_dir = "output_dat";
    std::filesystem::create_directories(dat_dir);
    std::cout << "Created DAT output directory: " << dat_dir << std::endl;

    // Save mesh if requested
    if (save_mesh) {
        // Save mesh as MFEM format
        std::ofstream mesh_ofs(fname+".mesh");
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
    }

    // Set up parallel finite element space
    dim = mesh->Dimension();
    mfem::H1_FECollection fec(order, dim);
    mfem::ParFiniteElementSpace fes(mesh, &fec, 1);

    HYPRE_BigInt global_dofs = fes.GlobalTrueVSize();
    if (myid == 0) {
        std::cout << "Number of finite element unknowns: " << global_dofs << std::endl;
    }

    mfem::ParGridFunction gf(&fes); // Parallel grid function

    // ================== Timing PHASE 1: DERIVATIVE COMPUTATION ==================
    mfem::StopWatch deriv_timer;
    deriv_timer.Start();
    // Initialise parallel DCPSE derivatives
    StreamVorti::ParDcpse2d* derivs = InitialiseParDCPSE(gf, dim, params.num_neighbors);

    deriv_timer.Stop();
    perf_metrics.derivative_time = deriv_timer.RealTime();

    std::cout << "Phase 1 - Derivatives computation: " << perf_metrics.derivative_time << " s" << std::endl;

    // save derivs matrices
    // TODO: SaveDerivativeMatrices not yet implemented for ParDcpse2d (HypreParMatrix output)
    // SaveDerivativeMatrices(derivs, params, dim, save_d, save_dd, dat_dir);

    // Save neighbors if requested
    if (save_neighbors) {
        std::cout << "main: Save neighbor indices to file... " << std::endl;
        derivs->SaveNeighsToFile(derivs->NeighborIndices(), dat_dir + "/" + params.output_prefix + ".neighbors" + params.output_extension);
        std::cout << "done." << std::endl;
    }


// ================================================================
// PART 1: INITIALIZATION
// ================================================================

    const int num_nodes = fes.GetTrueVSize();

    // Get parallel DCPSE derivative matrices (HypreParMatrix)
    const mfem::HypreParMatrix& dx_matrix = derivs->ShapeFunctionDx();
    const mfem::HypreParMatrix& dy_matrix = derivs->ShapeFunctionDy();
    const mfem::HypreParMatrix& dxx_matrix = derivs->ShapeFunctionDxx();
    const mfem::HypreParMatrix& dyy_matrix = derivs->ShapeFunctionDyy();

    std::cout << "Setup: Retrieved parallel DCPSE derivative matrices successfully." << std::endl;

    // Identify boundary and interior nodes
    // Use the parallel-specific version which correctly handles DOF ownership.
    std::vector<int> bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes;
    IdentifyBoundaryNodesPar(mesh, fes, bottom_nodes, right_nodes, top_nodes, left_nodes, interior_nodes);

    // Also identify boundary nodes by MFEM boundary attributes (for SDL BC integration)
    std::map<int, std::vector<int>> boundary_nodes_by_attr;
    std::vector<int> interior_nodes_by_attr;
    IdentifyBoundaryNodesByAttributePar(mesh, fes, boundary_nodes_by_attr, interior_nodes_by_attr);

    // Combine all boundary nodes for streamfunction boundary conditions
    std::vector<int> all_boundary_nodes;
    all_boundary_nodes.insert(all_boundary_nodes.end(), bottom_nodes.begin(), bottom_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), right_nodes.begin(), right_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), top_nodes.begin(), top_nodes.end());
    all_boundary_nodes.insert(all_boundary_nodes.end(), left_nodes.begin(), left_nodes.end());

#ifdef STREAMVORTI_WITH_ECL
    // Map true DOF indices back to vertex indices (needed for GetVertex in SDL BC application)
    std::map<int, int> tdof_to_vertex;
    for (int i = 0; i < mesh->GetNV(); ++i) {
        int tdof = fes.GetLocalTDofNumber(i);
        if (tdof >= 0) {
            tdof_to_vertex[tdof] = i;
        }
    }
#endif

    // ================================================================
    // STREAM FUNCTION BOUNDARY CONDITIONS
    // ================================================================
    // The Poisson equation for the stream function is: -∇²ψ = ω
    // (negated for positive definiteness with CG/direct solvers).
    //
    // Boundary conditions on ψ depend on BC type:
    //   no-slip/slip:      Dirichlet ψ = constant (determined by flow rate)
    //   velocity:          Dirichlet ψ = ∫u dy  (integrated from velocity profile)
    //   pressure/outflow:  Neumann ∂ψ/∂n = 0 (natural outflow)
    //
    // For Dirichlet nodes, EliminateRowsCols replaces the row with identity.
    // For Neumann nodes, the row is replaced with the normal derivative
    // operator via GetDiag/GetOffd: (∂ψ/∂n)_i = 0

    // Prescribed stream function values at boundary nodes.
    // Initialized to 0 (backward-compatible with cavity where ψ=0 everywhere).
    int num_true_dofs = fes.GetTrueVSize();
    mfem::Vector psi_bc(num_true_dofs);
    psi_bc = 0.0;

    // Lists of Dirichlet and Neumann boundary nodes for the ψ equation.
    std::vector<int> dirichlet_psi_nodes;
    std::vector<std::pair<int, char>> neumann_psi_info;  // (tdof_idx, axis)

#ifdef STREAMVORTI_WITH_ECL
    if (use_sdl_bcs) {
        // ============================================================
        // Classify boundary nodes and compute ψ values
        // ============================================================
        constexpr int nsub = 20;
        for (const auto& bc : sdl_boundaries) {
            auto it = boundary_nodes_by_attr.find(bc.attribute);
            if (it == boundary_nodes_by_attr.end()) continue;
            const auto& nodes = it->second;

            if (bc.type == "velocity") {
                // Compute ψ at each node by numerically integrating the
                // velocity function from 0 to the node's coordinate.
                for (int tdof : nodes) {
                    auto vit = tdof_to_vertex.find(tdof);
                    if (vit == tdof_to_vertex.end()) continue;
                    const double* vtx = mesh->GetVertex(vit->second);
                    double psi_val = 0.0;

                    if (bc.predicate_axis == 'x' && bc.u_function) {
                        double x0 = bc.predicate_value;
                        double y = vtx[1];
                        double h = y / nsub;
                        for (int k = 0; k < nsub; ++k) {
                            double a = k * h;
                            double b = (k + 1) * h;
                            double m = 0.5 * (a + b);
                            psi_val += (h / 6.0) * (
                                bc.u_function->evaluateAt(x0, a, 0.0) +
                                4.0 * bc.u_function->evaluateAt(x0, m, 0.0) +
                                bc.u_function->evaluateAt(x0, b, 0.0));
                        }
                    } else if (bc.predicate_axis == 'y' && bc.v_function) {
                        double y0 = bc.predicate_value;
                        double x = vtx[0];
                        double h = x / nsub;
                        for (int k = 0; k < nsub; ++k) {
                            double a = k * h;
                            double b = (k + 1) * h;
                            double m = 0.5 * (a + b);
                            psi_val -= (h / 6.0) * (
                                bc.v_function->evaluateAt(a, y0, 0.0) +
                                4.0 * bc.v_function->evaluateAt(m, y0, 0.0) +
                                bc.v_function->evaluateAt(b, y0, 0.0));
                        }
                    }
                    psi_bc[tdof] = psi_val;
                }

                for (int tdof : nodes) {
                    dirichlet_psi_nodes.push_back(tdof);
                }
            } else if (bc.type == "pressure" || bc.type == "outflow") {
                // Neumann ∂ψ/∂n = 0: natural outflow condition.
                char axis = bc.predicate_axis;
                if (axis == '\0') {
                    MFEM_ABORT("Pressure/outflow BC '" << bc.name
                               << "' requires a simple predicate (= x val) or (= y val) "
                               << "to determine the normal direction for Neumann ∂ψ/∂n = 0");
                }
                for (int tdof : nodes) {
                    neumann_psi_info.push_back({tdof, axis});
                }
            } else {
                // no-slip, slip: Dirichlet ψ = constant
                for (int tdof : nodes) {
                    dirichlet_psi_nodes.push_back(tdof);
                }
            }
        }

        // Apply wall ψ constants computed from serial mesh boundary connectivity
        for (const auto& bc : sdl_boundaries) {
            if (bc.type != "no-slip" && bc.type != "slip") continue;
            auto wit = wall_psi_by_attr.find(bc.attribute);
            if (wit == wall_psi_by_attr.end()) continue;
            auto nit = boundary_nodes_by_attr.find(bc.attribute);
            if (nit == boundary_nodes_by_attr.end()) continue;
            for (int tdof : nit->second) {
                psi_bc[tdof] = wit->second;
            }
        }
    } else
#endif
    {
        // CLI mode: all boundary nodes are Dirichlet with ψ = 0 (cavity)
        dirichlet_psi_nodes = all_boundary_nodes;
    }

    // ================================================================
    // LAPLACIAN MATRIX ASSEMBLY
    // ================================================================
    mfem::HypreParMatrix* laplacian_matrix = mfem::ParAdd(&dxx_matrix, &dyy_matrix);
    *laplacian_matrix *= -1.0;

    // STEP 1: Replace Neumann rows with derivative operators FIRST
    // (before Dirichlet column elimination).
    // For outlet at x=const: row ← dx_matrix row (enforces ∂ψ/∂x = 0)
    // For outlet at y=const: row ← dy_matrix row (enforces ∂ψ/∂y = 0)
    if (!neumann_psi_info.empty()) {
        mfem::Array<int> neumann_rows;
        for (const auto& [idx, axis] : neumann_psi_info) {
            neumann_rows.Append(idx);
        }
        laplacian_matrix->EliminateRows(neumann_rows);

        // Replace zeroed rows with derivative operators via GetDiag/GetOffd.
        // All DCPSE matrices share identical sparsity, so SetRow always
        // matches GetRow column count.
        for (const auto& [idx, axis] : neumann_psi_info) {
            const mfem::HypreParMatrix& deriv = (axis == 'x') ? dx_matrix : dy_matrix;

            // Diagonal block
            mfem::SparseMatrix A_diag, D_diag;
            laplacian_matrix->GetDiag(A_diag);
            deriv.GetDiag(D_diag);
            mfem::Array<int> cols;
            mfem::Vector vals;
            D_diag.GetRow(idx, cols, vals);
            A_diag.SetRow(idx, cols, vals);

            // Off-diagonal block
            mfem::SparseMatrix A_offd, D_offd;
            HYPRE_BigInt *cmap_A, *cmap_D;
            laplacian_matrix->GetOffd(A_offd, cmap_A);
            deriv.GetOffd(D_offd, cmap_D);
            D_offd.GetRow(idx, cols, vals);
            A_offd.SetRow(idx, cols, vals);
        }
    }

    // STEP 2: Eliminate Dirichlet DOFs via OperatorHandle.
    // EliminateRowsCols zeroes Dirichlet rows/columns in laplacian_matrix
    // and stores the eliminated column entries in Ae_h (used each time
    // step by EliminateBC to correct the RHS for known boundary values).
    mfem::Array<int> ess_tdof_list;
    for (int idx : dirichlet_psi_nodes) ess_tdof_list.Append(idx);

    mfem::OperatorHandle A_h(laplacian_matrix, false);
    mfem::OperatorHandle Ae_h;
    Ae_h.EliminateRowsCols(A_h, ess_tdof_list);

    std::cout << "Setup: Created parallel Laplacian matrix with boundary conditions." << std::endl;

    // ================== Timing PHASE 2: FACTORIZATION ==================
    mfem::StopWatch factor_timer;
    factor_timer.Start();

    // Create parallel Hypre solver (inline for HypreParMatrix compatibility)
    mfem::Solver* linear_solver = nullptr;
    mfem::Solver* precond = nullptr;

    if (params.solver_type == "cg") {
        if (myid == 0) std::cout << "Using HyprePCG solver (parallel conjugate gradient)" << std::endl;
        mfem::HyprePCG* cg_solver = new mfem::HyprePCG(MPI_COMM_WORLD);
        cg_solver->SetTol(params.solver_rel_tol);
        cg_solver->SetMaxIter(params.solver_max_iter);
        cg_solver->SetPrintLevel(params.solver_print_level);
        cg_solver->SetOperator(*laplacian_matrix);
        linear_solver = cg_solver;
    }
    else if (params.solver_type == "gmres") {
        if (myid == 0) std::cout << "Using HypreGMRES solver (parallel GMRES)" << std::endl;
        mfem::HypreGMRES* gmres_solver = new mfem::HypreGMRES(MPI_COMM_WORLD);
        gmres_solver->SetTol(params.solver_rel_tol);
        gmres_solver->SetMaxIter(params.solver_max_iter);
        gmres_solver->SetPrintLevel(params.solver_print_level);
        gmres_solver->SetOperator(*laplacian_matrix);
        linear_solver = gmres_solver;
    }
    else if (params.solver_type == "hypre" || params.solver_type == "amg") {
        if (myid == 0) std::cout << "Using HypreBoomerAMG solver (algebraic multigrid)" << std::endl;
        mfem::HypreBoomerAMG* amg_solver = new mfem::HypreBoomerAMG(*laplacian_matrix);
        amg_solver->SetPrintLevel(params.solver_print_level);
        linear_solver = amg_solver;
    }
    else {
        MFEM_ABORT("Parallel solver type '" << params.solver_type << "' not supported. Use: cg, gmres, hypre");
    }

    factor_timer.Stop();
    perf_metrics.factorization_time = factor_timer.RealTime();
    std::cout << "Phase 2 - Solver setup: " << perf_metrics.factorization_time << " s" << std::endl;

    // Initialize solution vectors
    mfem::Vector vorticity(num_nodes);
    mfem::Vector streamfunction(num_nodes);
    mfem::Vector vorticity_old(num_nodes); // for residual monitor
    mfem::Vector rhs(num_nodes);

    vorticity = 0.0;
    streamfunction = 0.0;

    // Derivative vectors (allocated once, reused throughout)
    // Streamfunction derivatives (for interior vorticity update)
    mfem::Vector dpsi_dy(num_nodes);      // ∂ψ/∂y
    mfem::Vector dpsi_dx(num_nodes);      // ∂ψ/∂x

    // CRITICAL FIX: Separate velocity vectors (for boundary conditions)
    mfem::Vector u_velocity(num_nodes);   // u = ∂ψ/∂y (with BCs applied)
    mfem::Vector v_velocity(num_nodes);   // v = -∂ψ/∂x (with BCs applied)

    // Velocity derivatives (for boundary vorticity calculation)
    mfem::Vector du_dy(num_nodes);        // ∂u/∂y
    mfem::Vector dv_dx(num_nodes);        // ∂v/∂x

    // Vorticity derivatives
    mfem::Vector domega_dx(num_nodes);
    mfem::Vector domega_dy(num_nodes);
    mfem::Vector d2omega_dx2(num_nodes);
    mfem::Vector d2omega_dy2(num_nodes);

    // Residual monitoring vectors
    mfem::Vector convection_term(num_nodes);
    mfem::Vector diffusion_term(num_nodes);

    // ====================================================================
    // PARAVIEW OUTPUT SETUP (Parallel)
    // ====================================================================
    // Create parallel finite element spaces for visualization
    mfem::ParFiniteElementSpace scalar_fes(mesh, &fec, 1);
    mfem::ParFiniteElementSpace vector_fes(mesh, &fec, mesh->Dimension());

    // Create parallel grid functions
    mfem::ParGridFunction vorticity_gf(&scalar_fes);
    mfem::ParGridFunction streamfunction_gf(&scalar_fes);
    mfem::ParGridFunction velocity_gf(&vector_fes);

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
    double current_dt = params.dt; //fixed dt
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
// PART 2: TIME-STEPPING LOOP
// ====================================================================
    // ================== Timing PHASE 3: SIMULATION  ==================
    mfem::StopWatch main_loop_timer;
    main_loop_timer.Start();

    // Main simulation loop (implementing Algorithm from Bourantas et al. 2019)

    for (int time_step = 1; time_step <= num_timesteps; ++time_step) {
        double current_time = time_step * current_dt;
        perf_metrics.total_iterations = time_step;

        vorticity_old = vorticity;

        // ================================================================
        // STEP 1: COMPUTE ALL DERIVATIVES (reused by all subsequent steps)
        // ================================================================
        dy_matrix.Mult(streamfunction, dpsi_dy);      // ∂ψ/∂y = u
        dx_matrix.Mult(streamfunction, dpsi_dx);      // ∂ψ/∂x = -v (will negate later)
        dx_matrix.Mult(vorticity, domega_dx);         // ∂ω/∂x
        dy_matrix.Mult(vorticity, domega_dy);         // ∂ω/∂y
        dxx_matrix.Mult(vorticity, d2omega_dx2);      // ∂²ω/∂x²
        dyy_matrix.Mult(vorticity, d2omega_dy2);      // ∂²ω/∂y²

        // ================================================================
        // STEP 2: UPDATE VORTICITY (using derivatives computed above) (interior only)
        // ================================================================
        // Use explicit Euler scheme (Eq. 11)
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)interior_nodes.size(); ++i) {
            int idx = interior_nodes[i];
            double convection = dpsi_dy[idx] * domega_dx[idx] - dpsi_dx[idx] * domega_dy[idx];
            double diffusion = (1.0 / params.reynolds_number) * (d2omega_dx2[idx] + d2omega_dy2[idx]);
            if (params.check_residuals) {
                convection_term[idx] = convection;
                diffusion_term[idx] = diffusion;
            }
            vorticity[idx] += current_dt * (diffusion - convection);
        }
#else
        for (int idx : interior_nodes) {
            double convection = dpsi_dy[idx] * domega_dx[idx] - dpsi_dx[idx] * domega_dy[idx];
            double diffusion = (1.0 / params.reynolds_number) * (d2omega_dx2[idx] + d2omega_dy2[idx]);
            if (params.check_residuals) {
                convection_term[idx] = convection;
                diffusion_term[idx] = diffusion;
            }
            vorticity[idx] += current_dt * (diffusion - convection);
        }
#endif


        // ================================================================
        // STEP 3: SOLVE STREAMFUNCTION POISSON EQUATION (Equation 12)
        // ================================================================
        // Solve -∇²ψ = ω instead ∇²ψ = -ω, so that RHS stays positive for linear solvers (cg)
        rhs = vorticity;

        // Apply boundary conditions to Poisson RHS:
        //   Neumann nodes: rhs = 0 (∂ψ/∂n = 0, natural outflow)
        //   Dirichlet nodes: EliminateBC sets rhs = psi_bc and corrects
        //   interior RHS for eliminated column contributions.
        for (const auto& [idx, axis] : neumann_psi_info) {
            rhs[idx] = 0.0;
        }
        A_h.EliminateBC(Ae_h, ess_tdof_list, psi_bc, rhs);

        linear_solver->Mult(rhs, streamfunction);

        // ================================================================
        // STEP 4: APPLY BOUNDARY CONDITIONS
        // ================================================================
        // Compute velocities from updated streamfunction:
        //   u = ∂ψ/∂y   (Eq. 3a)
        //   v = -∂ψ/∂x  (Eq. 3b)
        dy_matrix.Mult(streamfunction, u_velocity);  // u = ∂ψ/∂y
        dx_matrix.Mult(streamfunction, v_velocity);  // temp = ∂ψ/∂x
        v_velocity *= -1.0;                          // v = -∂ψ/∂x

#ifdef STREAMVORTI_WITH_ECL
        if (use_sdl_bcs) {
            // Apply boundary conditions from SDL using attribute-based lookup.
            // Each BC was assigned an MFEM boundary attribute at mesh setup
            // (sv_mesh_set_boundary_attribute_by_predicate); boundary_nodes_by_attr
            // maps attribute -> node list (built by IdentifyBoundaryNodesByAttributePar).
            //
            // Corner nodes appear in multiple attribute lists because MFEM
            // attributes are per-element (edge), not per-vertex. At corners,
            // two edges with different attributes share a vertex. The last BC
            // processed wins at these shared nodes — for lid-driven cavity this
            // is correct since all corners are no-slip regardless of order.
            //
            // BC type actions on velocity:
            //   no-slip:          u = 0, v = 0
            //   velocity (inlet): u = u_prescribed, v = v_prescribed
            //   pressure/outflow: leave computed values (natural outflow from ψ)
            //   slip:             leave computed values (ψ=const → zero normal vel)
            for (const auto& bc : sdl_boundaries) {
                auto it = boundary_nodes_by_attr.find(bc.attribute);
                if (it == boundary_nodes_by_attr.end()) continue;

                for (int vi : it->second) {
                    if (bc.type == "no-slip") {
                        u_velocity[vi] = 0.0;
                        v_velocity[vi] = 0.0;
                    } else if (bc.type == "velocity") {
                        if (bc.u_function) {
                            const double* vertex = mesh->GetVertex(tdof_to_vertex[vi]);
                            double x = vertex[0], y = vertex[1];
                            double z = (dim > 2) ? vertex[2] : 0.0;
                            u_velocity[vi] = bc.u_function->evaluateAt(x, y, z);
                        }
                        if (bc.v_function) {
                            const double* vertex = mesh->GetVertex(tdof_to_vertex[vi]);
                            double x = vertex[0], y = vertex[1];
                            double z = (dim > 2) ? vertex[2] : 0.0;
                            v_velocity[vi] = bc.v_function->evaluateAt(x, y, z);
                        }
                    } else if (bc.type == "pressure" || bc.type == "outflow") {
                        // Natural outflow: velocity from ψ derivatives stands
                    } else if (bc.type == "slip") {
                        // Free-slip: ψ=const gives zero normal velocity
                    } else {
                        MFEM_ABORT("Unknown boundary condition type '"
                                   << bc.type << "' for boundary '" << bc.name
                                   << "' at node " << vi);
                    }
                }
            }
        } else
#endif
        {
            // Fallback: hardcoded lid-driven cavity BCs
            for (int idx : bottom_nodes) { u_velocity[idx] = 0.0; v_velocity[idx] = 0.0; }  // u=0, v=0
            for (int idx : right_nodes)  { u_velocity[idx] = 0.0; v_velocity[idx] = 0.0; }  // u=0, v=0
            for (int idx : left_nodes)   { u_velocity[idx] = 0.0; v_velocity[idx] = 0.0; }  // u=0, v=0
            for (int idx : top_nodes)    { u_velocity[idx] = 1.0; v_velocity[idx] = 0.0; }  // u=1, v=0 (moving lid)
        }

        // Compute boundary vorticity: ω = ∂v/∂x - ∂u/∂y
        dy_matrix.Mult(u_velocity, du_dy);  // ∂u/∂y (u has BCs applied)
        dx_matrix.Mult(v_velocity, dv_dx);  // ∂v/∂x (v has BCs applied)

        // Set vorticity at ALL boundary nodes
        for (int idx : all_boundary_nodes) {
            vorticity[idx] = dv_dx[idx] - du_dy[idx];
        }

        // ================================================================
        // STEP 5: COMPUTE AND LOG RESIDUALS (if enabled)
        // ================================================================
        if (params.check_residuals && (time_step % params.residual_check_freq == 0)) {

            // Vorticity residual: ∂ω/∂t + K(ψ,ω) - L(ω) = 0
            double vort_sum_sq = 0.0, vort_max = 0.0;
            for (int idx : interior_nodes) {
                double domega_dt = (vorticity[idx] - vorticity_old[idx]) / current_dt;
                double conv = convection_term[idx];
                double diff = diffusion_term[idx];

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

                // Compute relative changes with proper threshold check
                double vort_norm = vorticity.Norml2();
                double stream_norm = streamfunction.Norml2();
                const double norm_threshold = 1e-12;

                vort_rel_change = (vort_norm > norm_threshold) ?
                                  diff_vort.Norml2() / vort_norm :
                                  diff_vort.Norml2();
                psi_rel_change = (stream_norm > norm_threshold) ?
                                 diff_stream.Norml2() / stream_norm :
                                 diff_stream.Norml2();
            }

            // Store in history
            residual_history.Add(time_step, current_time,
                                vort_rms, vort_max, psi_rms, psi_max,
                                vort_rel_change, psi_rel_change);
        }

        // ================================================================
        // STEP 6: ADAPTIVE TIMESTEP
        // Implement Gershgorin circle theorem estimation (Equation 42 from paper)
        // ================================================================
        // NOTE: Adaptive timestep disabled in parallel version due to HypreParMatrix
        // not supporting direct element access via GetI()/GetData()
        #ifndef MFEM_USE_MPI
        if (params.enable_adaptive_timestep &&
            time_step % params.gershgorin_frequency == 0) {
            // Compute Gershgorin bound: max_i { Σ_j (|L_ij| + |K_ij|) }
            double max_eigenvalue_bound = 0.0;

            for (int i = 0; i < num_nodes; ++i) {
                double row_sum = 0.0;

                // ============================================================
                // DIFFUSION OPERATOR: L = (1/Re) * (∂²/∂x² + ∂²/∂y²)
                // Accumulate |L_ij| = |(1/Re) * (dxx_matrix[i,j] + dyy_matrix[i,j])|
                // ============================================================

                // Contribution from ∂²/∂x² operator
                for (int j = dxx_matrix.GetI()[i]; j < dxx_matrix.GetI()[i + 1]; ++j) {
                    row_sum += std::abs((1.0 / params.reynolds_number) *
                                    dxx_matrix.GetData()[j]);
                }

                // Contribution from ∂²/∂y² operator
                for (int j = dyy_matrix.GetI()[i]; j < dyy_matrix.GetI()[i + 1]; ++j) {
                    row_sum += std::abs((1.0 / params.reynolds_number) *
                                    dyy_matrix.GetData()[j]);
                }

                // ============================================================
                // CONVECTION OPERATOR: K = u*(∂/∂x) + v*(∂/∂y)
                // Accumulate |K_ij| = |u_i * dx_matrix[i,j] + v_i * dy_matrix[i,j]|
                // ============================================================

                // Contribution from u*(∂/∂x) term
                for (int j = dx_matrix.GetI()[i]; j < dx_matrix.GetI()[i + 1]; ++j) {
                    row_sum += std::abs(dpsi_dy[i] * dx_matrix.GetData()[j]);
                }

                // Contribution from v*(∂/∂y) term (v = -∂ψ/∂x)
                for (int j = dy_matrix.GetI()[i]; j < dy_matrix.GetI()[i + 1]; ++j) {
                    row_sum += std::abs(-dpsi_dx[i] * dy_matrix.GetData()[j]);
                }

                // Track maximum row sum across all nodes
                max_eigenvalue_bound = std::max(max_eigenvalue_bound, row_sum);
            }

            // ============================================================
            // CRITICAL TIMESTEP: dt ≤ 2/|λ_max|  (Equation 40)
            // ============================================================
            double dt_critical = (max_eigenvalue_bound > 0.0) ?
                                2.0 / max_eigenvalue_bound : 1e-4;

            // Apply safety factor and adjust timestep if necessary
            if (current_dt > dt_critical / params.gershgorin_safety_factor) {
                current_dt = dt_critical / params.gershgorin_safety_factor;
                num_timesteps = static_cast<int>(params.final_time / current_dt);

                std::cout << "Adaptive timestep adjustment:" << std::endl;
                std::cout << "  |λ_max| ≈ " << std::scientific << std::setprecision(4)
                        << max_eigenvalue_bound << std::endl;
                std::cout << "  dt_critical = " << dt_critical << std::endl;
                std::cout << "  dt_new = " << current_dt
                        << " (safety factor: " << params.gershgorin_safety_factor << ")"
                        << std::endl;
                std::cout << "  Updated total steps = " << num_timesteps << std::endl;
            }
        }
        #endif // MFEM_USE_MPI

        // ================================================================
        // STEP 7: OUTPUT SOLUTION PERIODICALLY (if enabled)
        // ================================================================
        if (save_solutions && (time_step % params.output_frequency == 0)) {
            SaveSolutionToFile(vorticity, streamfunction, u_velocity, v_velocity,
                            params.output_prefix, time_step, dat_dir, false);
        }

        // ================================================================
        // STEP 8: PARAVIEW OUTPUT (if enabled)
        // ================================================================
        if (paraview_output && (time_step % params.paraview_output_frequency == 0)) {

            for (int i = 0; i < num_nodes; ++i) {
                vorticity_gf[i] = vorticity[i];
                streamfunction_gf[i] = streamfunction[i];
                velocity_gf[i] = u_velocity[i];                    // u
                velocity_gf[i + num_nodes] = v_velocity[i];      // v = -∂ψ/∂x
            }

            paraview_dc.SetTime(current_time);
            paraview_dc.SetCycle(time_step);
            paraview_dc.Save();
        }

        // ================================================================
        // STEP 9: CHECK STEADY STATE
        // ================================================================
        if (time_step >= 100 && (time_step % params.steady_state_freq == 0)) {

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

                // Compute relative changes with proper threshold check
                double vort_norm = vorticity.Norml2();
                double stream_norm = streamfunction.Norml2();
                const double norm_threshold = 1e-12;

                double vort_change = (vort_norm > norm_threshold) ?
                                     diff_vort.Norml2() / vort_norm :
                                     diff_vort.Norml2();
                double stream_change = (stream_norm > norm_threshold) ?
                                       diff_stream.Norml2() / stream_norm :
                                       diff_stream.Norml2();
                double max_change = std::max(vort_change, stream_change);

                // Report progress
                std::cout << "Step " << time_step << " (t=" << std::fixed
                        << std::setprecision(6) << current_time
                        << "): Temporal ‖Δ‖/‖·‖ = " << std::scientific << max_change;

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
                    std::cout << " (not converged)" << std::endl;
                    steady_consecutive_passes = 0; // Reset counter
                }

                // Update previous values for next check
                steady_prev_vorticity = vorticity;
                steady_prev_streamfunction = streamfunction;
            }
        }

        // ====================================================================
        // SIMPLE REAL-TIME MONITORING
        // ====================================================================

        // Ranges of Vorticity and Streamfunction results
        if (time_step % 1000 == 0) {
            std::cout << "\nStep " << time_step
                    << " | t=" << std::fixed << std::setprecision(2) << current_time
                    << " | ω: [" << std::setprecision(3) << vorticity.Min()
                    << ", " << vorticity.Max() << "]"
                    << " | ψ: [" << streamfunction.Min()
                    << ", " << streamfunction.Max() << "]" << std::endl;
        }

    } // End of time-stepping loop

    main_loop_timer.Stop();
    double total_solve_time = main_loop_timer.RealTime();

    std::cout << "\nSimulation completed successfully in "
            << main_loop_timer.RealTime() << " seconds." << std::endl;

// ====================================================================
// PART 3: POST-PROCESSING
// ====================================================================

    // Extract line probe data for validation
#ifdef STREAMVORTI_WITH_ECL
    if (!sdl_probes.empty()) {
        // Use SDL-configured probes
        for (const auto& probe : sdl_probes) {
            std::string filename = dat_dir + "/" + params.output_prefix
                                   + "_" + probe.name + ".dat";
            ExtractCenterline(u_velocity, v_velocity, &scalar_fes,
                              filename, probe.axis, probe.position);
        }
    } else
#endif
    {
        // Default probes for CLI mode (cavity at x=0.5, y=0.5)
        ExtractCenterline(u_velocity, v_velocity, &scalar_fes,
                          dat_dir + "/" + params.output_prefix + "_u_centerline_x0.5.dat",
                          'x', 0.5);
        ExtractCenterline(u_velocity, v_velocity, &scalar_fes,
                          dat_dir + "/" + params.output_prefix + "_v_centerline_y0.5.dat",
                          'y', 0.5);
    }

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
    std::cout << "  Total iterations(time steps) completed: " << perf_metrics.total_iterations << "\n";
    std::cout << "  Average time per iteration: " << total_solve_time / perf_metrics.total_iterations << " s\n\n";


    // Final solution statistics - compute GLOBAL max/min across all ranks
    double local_vort_max = vorticity.Max();
    double local_vort_min = vorticity.Min();
    double local_psi_max = streamfunction.Max();
    double local_psi_min = streamfunction.Min();

    double global_vort_max, global_vort_min;
    double global_psi_max, global_psi_min;

    MPI_Reduce(&local_vort_max, &global_vort_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_vort_min, &global_vort_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_psi_max, &global_psi_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_psi_min, &global_psi_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        std::cout << "Final simulation statistics (GLOBAL):" << std::endl;
        std::cout << "  - Maximum vorticity: " << global_vort_max << std::endl;
        std::cout << "  - Minimum vorticity: " << global_vort_min << std::endl;
        std::cout << "  - Maximum streamfunction: " << global_psi_max << std::endl;
        std::cout << "  - Minimum streamfunction: " << global_psi_min << std::endl;
    }

#ifdef STREAMVORTI_WITH_ECL
    // Shutdown ECL if it was initialized
    if (sdl_mode) {
        StreamVorti::Lisp::Runtime::shutdown();
    }
#endif

    if (myid == 0) {
        std::cout << "main: success!" << std::endl;
    }

    // Clean up all Hypre/MPI objects before MPI_Finalize (Hypre requires this)
    delete derivs;
    delete laplacian_matrix;
    delete linear_solver;
    delete mesh;

    MPI_Finalize();
    return EXIT_SUCCESS;
}

// ============================================================================
// FUNCTION IMPLEMENTATIONS
// ============================================================================

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

StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int num_neighbors) {
    MFEM_ABORT("Serial InitialiseDCPSE not available in parallel build. "
               "Use InitialiseParDCPSE() which creates ParDcpse2d with "
               "distributed ghost exchange.");
    return nullptr;
}

#ifdef MFEM_USE_MPI
StreamVorti::ParDcpse2d* InitialiseParDCPSE(mfem::ParGridFunction& gf, int dim, int num_neighbors) {
    std::cout << "InitialiseParDCPSE: Initialising parallel DC PSE derivatives." << std::endl;
    mfem::StopWatch timer;
    timer.Start();

    if (dim != 2) {
        MFEM_ABORT("InitialiseParDCPSE: Only 2D supported currently. Got dim=" << dim);
    }

    StreamVorti::ParDcpse2d* derivs = new StreamVorti::ParDcpse2d(gf, num_neighbors);

    std::cout << "InitialiseParDCPSE: Execution time for parallel DCPSE initialisation: "
              << timer.RealTime() << " s" << std::endl;

    timer.Clear();
    derivs->Update();  // Includes ghost exchange via ParSupportDomain::Update()
    std::cout << "InitialiseParDCPSE: Execution time for parallel DCPSE computation (includes ghost exchange): "
              << timer.RealTime() << " s" << std::endl;

    return derivs;
}
#endif

void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                           int dim, bool save_d, bool save_dd,
                           const std::string& dat_dir) {
    MFEM_ABORT("Serial SaveDerivativeMatrices not available in parallel build. "
               "Parallel DCPSE derivatives use HypreParMatrix; parallel export "
               "not yet implemented.");
}

void IdentifyBoundaryNodes(mfem::Mesh* mesh,
                                          std::vector<int>& bottom_nodes,
                                          std::vector<int>& right_nodes,
                                          std::vector<int>& top_nodes,
                                          std::vector<int>& left_nodes,
                                          std::vector<int>& interior_nodes,
                                          int num_true_verts) {
    MFEM_ABORT("Serial IdentifyBoundaryNodes not available in parallel build. "
               "Use IdentifyBoundaryNodesPar() which uses GetLocalTDofNumber() "
               "for correct DOF ownership.");
}

void SaveSolutionToFile(const mfem::Vector& vorticity,
                       const mfem::Vector& streamfunction,
                       const mfem::Vector& u_velocity,
                       const mfem::Vector& v_velocity,
                       const std::string& filename,
                       int timestep,
                       const std::string& dat_dir,
                       bool is_final)
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

SolverPackage* CreateLinearSolver(const std::string& solver_type,
                                   const mfem::SparseMatrix& A,
                                   const SimulationParams& params)
{
    MFEM_ABORT("Serial CreateLinearSolver not available in parallel build. "
               "Parallel solvers (HyprePCG, HypreGMRES, HypreBoomerAMG) are "
               "created inline with HypreParMatrix.");
    return nullptr;
}

std::string GetCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

void IdentifyBoundaryNodesByAttribute(
    mfem::Mesh* mesh,
    std::map<int, std::vector<int>>& boundary_nodes,
    std::vector<int>& interior_nodes,
    int num_true_verts)
{
    MFEM_ABORT("Serial IdentifyBoundaryNodesByAttribute not available in parallel build. "
               "Use IdentifyBoundaryNodesByAttributePar() which uses "
               "GetLocalTDofNumber() for correct DOF ownership.");
}

void ExtractCenterline(const mfem::Vector& u_velocity,
                       const mfem::Vector& v_velocity,
                       mfem::ParFiniteElementSpace* pfes,
                       const std::string& filename,
                       char axis,
                       double position,
                       double tol)
{
    MPI_Comm comm = pfes->GetComm();
    int myrank, nranks;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nranks);

    mfem::Mesh* mesh = pfes->GetMesh();
    int nldofs = pfes->GetNDofs();  // local DOFs (owned + shared)

    // Collect owned nodes along the line: iterate LOCAL DOFs, skip non-owned.
    // For P1 H1 elements: local DOF ldof = local vertex ldof.
    // True DOF index = pfes->GetLocalTDofNumber(ldof).
    // u_velocity[tdof] is the correct velocity at the owned vertex.
    std::vector<double> local_vary, local_u, local_v;
    for (int ldof = 0; ldof < nldofs; ++ldof) {
        int tdof = pfes->GetLocalTDofNumber(ldof);
        if (tdof < 0) continue;  // shared vertex owned by another rank

        const double* vertex = mesh->GetVertex(ldof);
        double fixed_val = (axis == 'x') ? vertex[0] : vertex[1];
        double vary_val  = (axis == 'x') ? vertex[1] : vertex[0];

        if (std::abs(fixed_val - position) < tol) {
            local_vary.push_back(vary_val);
            local_u.push_back(u_velocity[tdof]);
            local_v.push_back(v_velocity[tdof]);
        }
    }

    // Gather all centerline data to rank 0 via MPI
    int local_n = (int)local_vary.size();
    std::vector<int> counts(nranks, 0), displs(nranks, 0);
    MPI_Gather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

    int total_n = 0;
    if (myrank == 0) {
        for (int r = 0; r < nranks; ++r) { displs[r] = total_n; total_n += counts[r]; }
    }

    std::vector<double> all_vary(total_n), all_u(total_n), all_v(total_n);
    MPI_Gatherv(local_vary.data(), local_n, MPI_DOUBLE,
                all_vary.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
    MPI_Gatherv(local_u.data(), local_n, MPI_DOUBLE,
                all_u.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);
    MPI_Gatherv(local_v.data(), local_n, MPI_DOUBLE,
                all_v.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    // Only rank 0 writes the file
    if (myrank == 0) {
        // Sort by varying coordinate
        std::vector<int> order(total_n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return all_vary[a] < all_vary[b]; });

        std::ofstream out(filename);
        out << "# StreamVorti Centerline Data\n";
        out << "# Axis: " << axis << " = " << position << "\n";
        out << (axis == 'x' ? "# y  u  v\n" : "# x  u  v\n");
        out.precision(8);
        for (int idx : order) {
            out << all_vary[idx] << " " << all_u[idx] << " " << all_v[idx] << "\n";
        }
        std::cout << "Saved centerline data to: " << filename
                  << " (" << total_n << " points)" << std::endl;
    }
}

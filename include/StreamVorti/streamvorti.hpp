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
 *      Weizheng(Will) LI
 */



/*!
   \file stream_vorti.hpp
   \brief StreamVorti software header file. Should be included in a project to use StreamVorti.
   \author Konstantinos A. Mountris
   \date 12/01/2018
   \author Weizheng Li
   \date 17/10/2025
*/

#ifndef STREAMVORTI_STREAM_VORTI_SIM_HPP_
#define STREAMVORTI_STREAM_VORTI_SIM_HPP_

// Collecting StreamVorti modules' header files.

#include "StreamVorti/approximants/dcpse.hpp"
#include "StreamVorti/approximants/dcpse_2d.hpp"
#include "StreamVorti/approximants/dcpse_3d.hpp"
#include "StreamVorti/support_domain/support_domain.hpp"

// Standard library includes for struct definitions
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Forward declarations
namespace mfem {
    class Mesh;
    class GridFunction;
    class SparseMatrix;
    class Vector;
    class Solver;
}

/*!
 * \namespace StreamVorti Classes collection for implementation of
 *            the Strong-Form Meshless Stream Function - Vorticity formulation.
 */

// ============================================================================
// SIMULATION CONFIGURATION STRUCTURES
// ============================================================================

/**
 * @brief Simulation parameters structure containing all configuration settings
 *
 * Encapsulates physical parameters, solver options, output settings, and
 * convergence criteria for the stream function-vorticity solver.
 */
struct SimulationParams {
    // Physical parameters
    double final_time = 60.0;              ///< Simulation end time [seconds]
    double dt = 1e-3;                      ///< Timestep [s] (CFL-limited, typical: 1e-4 to 1e-3)
    int reynolds_number = 1000;            ///< Reynolds number (range: 100-5000)
    int num_neighbors = 25;                ///< DC PSE stencil size (25-40 recommended)

    // Output file settings
    std::string output_prefix = "mfem_square10x10";  ///< Base name for output files
    std::string output_extension = ".dat";           ///< File extension for output files

    // Solution output frequencies
    int output_frequency = 100;            ///< Save solution every N steps (0=disable)
    int paraview_output_frequency = 1000;  ///< Save ParaView VTU every N steps (0=disable)

    // Adaptive timestepping (Gershgorin criterion)
    int gershgorin_frequency = 1000;       ///< Check Gershgorin stability every N steps
    bool enable_adaptive_timestep = false; ///< Enable adaptive dt (Gershgorin bound)
    double gershgorin_safety_factor = 5.0; ///< Safety factor for dt_critical (typ: 2-10)

    // Linear solver configuration
    std::string solver_type = "umfpack";   ///< Solver: umfpack/klu/cg/gmres/minres/bicgstab
    int solver_max_iter = 500;             ///< Maximum iterations for iterative solvers
    double solver_rel_tol = 1e-8;          ///< Relative tolerance (typical: 1e-6 to 1e-10)
    int solver_print_level = 0;            ///< Verbosity (0=quiet, 1=summary, 2=verbose)

    // Steady-state convergence detection
    double steady_state_tol = 1e-12;       ///< L2 relative change threshold (typical: 1e-5)
    int steady_state_freq = 100;           ///< Check steady state every N timesteps
    int steady_state_checks = 3;           ///< Consecutive passes required to declare convergence

    // Residual monitoring and diagnostics
    bool check_residuals = false;          ///< Enable residual computation and logging
    int residual_check_freq = 100;         ///< Compute residuals every N timesteps
};

/**
 * @brief Structure to store residual history during time-stepping
 *
 * Tracks vorticity and streamfunction residuals over time for convergence
 * monitoring and diagnostics. Supports CSV export and summary statistics.
 */
struct ResidualHistory {
    std::vector<int> timesteps;                   ///< Timestep numbers when residuals computed
    std::vector<double> times;                    ///< Physical times [s] at residual checks
    std::vector<double> vorticity_rms;            ///< RMS residual of vorticity equation
    std::vector<double> vorticity_max;            ///< Max absolute residual of vorticity
    std::vector<double> streamfunction_rms;       ///< RMS residual of Poisson (∇²ψ+ω=0)
    std::vector<double> streamfunction_max;       ///< Max absolute residual of Poisson
    std::vector<double> vorticity_change;         ///< Relative L2 change in vorticity
    std::vector<double> streamfunction_change;    ///< Relative L2 change in streamfunction

    /**
     * @brief Add a residual data point to the history
     * @param step Timestep number
     * @param time Physical time
     * @param vort_rms RMS residual of vorticity equation
     * @param vort_max Maximum absolute residual of vorticity
     * @param psi_rms RMS residual of streamfunction equation
     * @param psi_max Maximum absolute residual of streamfunction
     * @param vort_rel_change Relative change in vorticity (for steady-state monitoring)
     * @param psi_rel_change Relative change in streamfunction (for steady-state monitoring)
     */
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

    /**
     * @brief Save residual history to CSV file
     * @param filename Output file path
     */
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

    /**
     * @brief Print summary statistics of residual history
     *
     * Displays initial/final residuals and detects stagnant convergence
     */
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

/**
 * @brief Structure to hold performance timing metrics
 *
 * Tracks computational costs of different phases: derivative computation,
 * matrix factorization, and per-iteration solve time.
 */
struct PerformanceMetrics {
    double derivative_time = 0.0;          ///< DC PSE operator time [s] (Phase 1)
    double factorization_time = 0.0;       ///< Matrix factorization time [s] (Phase 2)
    double solution_time_per_iter = 0.0;   ///< Avg time per timestep [s] (Phase 3)
    int grid_size = 0;                     ///< Grid resolution (nx = ny for square grids)
    int total_iterations = 0;              ///< Total number of timesteps completed

    /**
     * @brief Print formatted performance metrics table
     * @param metrics The performance metrics to display
     */
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
        std::string grid_str = std::to_string(metrics.grid_size) + " x " +
                               std::to_string(metrics.grid_size);
        std::cout << std::setw(10) << grid_str
                  << std::setw(17) << std::fixed << std::setprecision(4)
                  << metrics.derivative_time << " s"
                  << std::setw(17) << metrics.factorization_time << " s"
                  << std::setw(17) << metrics.solution_time_per_iter << " s\n";
        std::cout << "========================================================================\n\n";
    }
};

/**
 * @brief Struct for linear solver and optional preconditioner
 *
 * Manages lifetime of solver and preconditioner pointers together.
 * Ensures proper cleanup and prevents memory leaks.
 * Preconditioner is nullptr for direct solvers (UMFPACK, KLU).
 */
struct SolverPackage {
    mfem::Solver* solver = nullptr;         ///< Linear solver pointer (owned, deleted in destructor)
    mfem::Solver* preconditioner = nullptr; ///< Preconditioner pointer (nullptr for direct solvers)

    /**
     * @brief Destructor: automatically cleans up solver and preconditioner
     *
     * Ensures proper cleanup even if exceptions occur.
     */
    ~SolverPackage() {
        delete solver;
        delete preconditioner;
    }
};

// ============================================================================
// FUNCTION DECLARATIONS
// ============================================================================

/**
 * @brief Create a new Cartesian mesh or load from file
 *
 * Generates uniform quadrilateral (2D) or hexahedral (3D) mesh if no file
 * is specified. Otherwise loads mesh from MFEM format file.
 *
 * @param mesh_file Path to mesh file (empty string to generate new mesh)
 * @param dim Spatial dimension (2 or 3)
 * @param nx Number of elements in x-direction
 * @param ny Number of elements in y-direction
 * @param nz Number of elements in z-direction
 * @param sx Physical size in x-direction
 * @param sy Physical size in y-direction
 * @param sz Physical size in z-direction
 * @param save_mesh Whether to save generated mesh to file
 * @return Pointer to created/loaded mesh
 */
mfem::Mesh* CreateOrLoadMesh(const char* mesh_file, int dim, int nx, int ny, int nz,
                            double sx, double sy, double sz, bool save_mesh);

/**
 * @brief Initialize DC PSE derivative operators
 *
 * Creates 2D or 3D DC PSE object, performs neighbor search, and computes
 * derivative matrices. Timing information is printed to console.
 *
 * @param gf Grid function containing nodal coordinates
 * @param dim Spatial dimension (2 or 3)
 * @param num_neighbors Number of neighbors for DC PSE stencil
 * @return Pointer to DC PSE derivative object
 */
StreamVorti::Dcpse* InitialiseDCPSE(mfem::GridFunction& gf, int dim, int num_neighbors);

/**
 * @brief Save DC PSE derivative matrices to file
 *
 * Exports sparse matrices for first derivatives (dx, dy, dz) and second
 * derivatives (dxx, dxy, dxz, dyy, dyz, dzz) based on dimension and flags.
 *
 * @param derivs DC PSE derivative object
 * @param params Simulation parameters (for output naming)
 * @param dim Spatial dimension
 * @param save_d Save first derivatives (gradient)
 * @param save_dd Save second derivatives (Hessian)
 * @param dat_dir Output directory path
 */
void SaveDerivativeMatrices(StreamVorti::Dcpse* derivs, const SimulationParams& params,
                            int dim, bool save_d, bool save_dd,
                            const std::string& dat_dir);

/**
 * @brief Identify boundary and interior nodes for lid-driven cavity
 *
 * Classifies mesh vertices into boundary segments (bottom, right, top, left)
 * and interior nodes. Corner nodes are included in horizontal boundaries (top/bottom).
 *
 * @param mesh The mesh to analyze
 * @param bottom_nodes Output: nodes on bottom boundary (y=0)
 * @param right_nodes Output: nodes on right boundary (x=1, excluding corners)
 * @param top_nodes Output: nodes on top boundary (y=1)
 * @param left_nodes Output: nodes on left boundary (x=0, excluding corners)
 * @param interior_nodes Output: nodes not on any boundary
 */
void IdentifyBoundaryNodes(mfem::Mesh* mesh, std::vector<int>& bottom_nodes,
                            std::vector<int>& right_nodes, std::vector<int>& top_nodes,
                            std::vector<int>& left_nodes, std::vector<int>& interior_nodes);

/**
 * @brief Create and configure linear solver based on type
 *
 * Supports direct solvers (UMFPACK, KLU) and iterative solvers (CG, GMRES,
 * MINRES, BiCGSTAB, HYPRE AMG). Iterative solvers use Gauss-Seidel
 * preconditioner except MINRES (no preconditioner for SPD systems).
 *
 * @param solver_type Solver name: "umfpack", "klu", "cg", "gmres", "minres", "bicgstab", "hypre"
 * @param A Coefficient matrix (streamfunction Laplacian)
 * @param params Simulation parameters (tolerances, max iterations, print level)
 * @return Pointer to solver package (solver + optional preconditioner)
 *
 * @note For this Poisson equation, the matrix is SPD so CG is optimal.
 *       MINRES is designed for indefinite systems but works.
 *       HypreBoomerAMG solver requires HypreParMatrix.
 */
SolverPackage* CreateLinearSolver(const std::string& solver_type, const mfem::SparseMatrix& A,
                                   const SimulationParams& params);

/**
 * @brief Save solution fields to DAT files
 *
 * Writes vorticity, streamfunction, and velocity components to separate
 * files with high precision (12 digits). Filename suffix indicates whether
 * this is final solution or intermediate output.
 *
 * @param vorticity Vorticity field
 * @param streamfunction Streamfunction field
 * @param u_velocity Horizontal velocity (u = ∂ψ/∂y)
 * @param v_velocity Vertical velocity (v = -∂ψ/∂x)
 * @param filename Base filename prefix
 * @param timestep Current timestep number
 * @param dat_dir Output directory
 * @param is_final True for final solution, false for intermediate
 */
void SaveSolutionToFile(const mfem::Vector& vorticity, const mfem::Vector& streamfunction,
                       const mfem::Vector& u_velocity, const mfem::Vector& v_velocity,
                       const std::string& filename, int timestep,
                       const std::string& dat_dir, bool is_final);

/**
 * @brief Get current date and time as formatted string
 *
 * Uses system clock to generate timestamp for logging and output filenames.
 *
 * @return String in format "YYYY-MM-DD HH:MM:SS"
 */
std::string GetCurrentDateTime();

/**
 * @brief General-purpose boundary node identification using MFEM boundary attributes
 *
 * Instead of hardcoding x=0, x=1, y=0, y=1, this uses the boundary attributes
 * assigned by the mesh generator or SDL file. Works with any geometry.
 *
 * @param mesh The MFEM mesh
 * @param boundary_nodes Output: Map from boundary attribute to list of node indices
 * @param interior_nodes Output: List of interior node indices
 */
void IdentifyBoundaryNodesByAttribute(
    mfem::Mesh* mesh,
    std::map<int, std::vector<int>>& boundary_nodes,
    std::vector<int>& interior_nodes);

/**
 * @brief Extract velocity along a centerline for validation
 *
 * Used to compare simulation results against Ghia et al. (1982) or
 * Erturk & Corke (2005) benchmark data.
 *
 * @param u_velocity U-velocity field
 * @param v_velocity V-velocity field
 * @param mesh The MFEM mesh
 * @param filename Output file path
 * @param axis Fixed coordinate ('x' or 'y')
 * @param position Value of the fixed coordinate (e.g., 0.5 for x=0.5)
 * @param tol Tolerance for "on the line" detection
 */
void ExtractCenterline(const mfem::Vector& u_velocity,
                       const mfem::Vector& v_velocity,
                       mfem::Mesh* mesh,
                       const std::string& filename,
                       char axis,
                       double position,
                       double tol = 0.01);


#endif //STREAMVORTI_STREAM_VORTI_SIM_HPP_

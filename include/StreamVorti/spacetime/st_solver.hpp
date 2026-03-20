/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Main space-time Navier-Stokes solver class.
 *
 * Implements the time-slab marching algorithm:
 *   1. Build space-time slab mesh from spatial mesh
 *   2. Set up FE spaces on the slab
 *   3. Solve nonlinear NS system on the slab
 *   4. Extract solution at top of slab (new initial condition)
 *   5. Optionally adapt mesh based on error indicators
 *   6. Advance to next slab
 *
 * Supports both serial and parallel execution.
 */

#ifndef STREAMVORTI_SPACETIME_SOLVER_HPP
#define STREAMVORTI_SPACETIME_SOLVER_HPP

#include "st_config.hpp"
#include "st_mesh.hpp"
#include "st_integrators.hpp"
#include "mfem.hpp"
#include <memory>
#include <vector>

namespace StreamVorti {
namespace SpaceTime {

/// Main space-time Navier-Stokes solver.
///
/// Solves incompressible Navier-Stokes equations using a time-slab
/// approach with stabilized space-time finite elements.
///
/// Usage:
///   SolverConfig config;
///   config.reynolds = 100;
///   config.mesh_file = "cylinder.mesh";
///   ...
///   STNavierStokesSolver solver(config);
///   solver.Initialize();
///   solver.Solve();
class STNavierStokesSolver
{
public:
    /// Construct solver with given configuration
    STNavierStokesSolver(const SolverConfig &config);

    /// Construct parallel solver
    STNavierStokesSolver(const SolverConfig &config, MPI_Comm comm);

    ~STNavierStokesSolver();

    /// Initialize the solver (load mesh, set up spaces, etc.)
    void Initialize();

    /// Run the full simulation (all time slabs)
    void Solve();

    /// Solve a single time slab [t_n, t_{n+1}]
    /// Returns the solution at t_{n+1} as initial condition for next slab
    void SolveTimeSlab(double t_start, double t_end,
                       const mfem::Vector &u_initial);

    /// Set body force coefficient
    void SetBodyForce(mfem::VectorCoefficient *f) { body_force_ = f; }

    /// Set initial condition for velocity
    void SetInitialVelocity(mfem::VectorCoefficient &u0);

    /// Set initial condition for pressure
    void SetInitialPressure(mfem::Coefficient &p0);

    /// Get the current velocity solution
    const mfem::GridFunction &GetVelocity() const { return *velocity_gf_; }

    /// Get the current pressure solution
    const mfem::GridFunction &GetPressure() const { return *pressure_gf_; }

    /// Get the spatial mesh
    mfem::Mesh *GetSpatialMesh() { return spatial_mesh_; }

    /// Get current simulation time
    double GetTime() const { return current_time_; }

    /// Get the solver configuration
    const SolverConfig &GetConfig() const { return config_; }

private:
    // ---- Initialization helpers ----
    void LoadSpatialMesh();
    void SetupFESpaces();
    void SetupSolver();
    void SetupOutput();

    // ---- Time-slab operations ----

    /// Build the space-time mesh for a slab
    void BuildSpaceTimeSlab(double t_start, double t_end);

    /// Set up FE spaces on the space-time slab
    void SetupSlabFESpaces();

    /// Assemble the space-time NS system on the slab
    void AssembleSlabSystem();

    /// Solve the nonlinear system on the slab
    void SolveSlabSystem();

    /// Extract solution at top of slab (t = t_{n+1})
    void ExtractTopSolution(mfem::Vector &u_top, mfem::Vector &p_top);

    /// Project initial condition onto bottom of slab
    void ProjectInitialCondition(const mfem::Vector &u_prev);

    // ---- Adaptive refinement ----

    /// Compute error indicators on the space-time slab
    void ComputeErrorIndicators(mfem::Vector &indicators);

    /// Adapt the spatial mesh based on error indicators
    bool AdaptSpatialMesh(const mfem::Vector &indicators);

    /// Compute element-wise residual error indicator
    void ComputeResidualIndicator(mfem::Vector &indicators);

    /// Compute velocity gradient error indicator
    void ComputeVelocityGradientIndicator(mfem::Vector &indicators);

    /// Compute vorticity-based error indicator
    void ComputeVorticityIndicator(mfem::Vector &indicators);

    // ---- Output ----

    /// Save the current solution to ParaView files
    void SaveOutput(int slab_index, double time);

    /// Save the adapted mesh
    void SaveMesh(int slab_index);

    // ---- Members ----

    SolverConfig config_;
    bool is_parallel_;
    MPI_Comm comm_;
    int myid_;
    int num_procs_;

    // Spatial mesh and FE
    mfem::Mesh *spatial_mesh_;
    mfem::FiniteElementCollection *vel_fec_;
    mfem::FiniteElementCollection *pres_fec_;
    mfem::FiniteElementSpace *vel_fes_;
    mfem::FiniteElementSpace *pres_fes_;

    // Space-time slab
    std::unique_ptr<SpaceTimeMesh> st_mesh_;
    mfem::FiniteElementCollection *st_vel_fec_;
    mfem::FiniteElementCollection *st_pres_fec_;
    mfem::FiniteElementSpace *st_vel_fes_;
    mfem::FiniteElementSpace *st_pres_fes_;

    // Solution fields (on spatial mesh, for output)
    mfem::GridFunction *velocity_gf_;
    mfem::GridFunction *pressure_gf_;
    mfem::GridFunction *vorticity_gf_;

    // Space-time solution (on slab)
    mfem::BlockVector *st_solution_;

    // NS operator and solver
    std::unique_ptr<STNavierStokesOperator> ns_operator_;
    mfem::NewtonSolver *newton_solver_;
    mfem::Solver *linear_solver_;
    mfem::Solver *preconditioner_;

    // Body force
    mfem::VectorCoefficient *body_force_;

    // Initial conditions
    mfem::Vector u_initial_;
    mfem::Vector p_initial_;

    // Output
    mfem::ParaViewDataCollection *paraview_dc_;

    // State
    double current_time_;
    int slab_count_;

#ifdef MFEM_USE_MPI
    // Parallel spatial mesh and FE
    mfem::ParMesh *par_spatial_mesh_;
    mfem::ParFiniteElementSpace *par_vel_fes_;
    mfem::ParFiniteElementSpace *par_pres_fes_;
#endif
};

} // namespace SpaceTime
} // namespace StreamVorti

#endif // STREAMVORTI_SPACETIME_SOLVER_HPP

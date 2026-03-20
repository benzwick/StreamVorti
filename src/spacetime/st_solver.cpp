/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Main space-time solver implementation.
 */

#include <StreamVorti/spacetime/st_solver.hpp>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cassert>

namespace StreamVorti {
namespace SpaceTime {

// =========================================================================
// Constructor (serial)
// =========================================================================

STNavierStokesSolver::STNavierStokesSolver(const SolverConfig &config)
    : config_(config),
      is_parallel_(false),
      comm_(MPI_COMM_WORLD),
      myid_(0),
      num_procs_(1),
      spatial_mesh_(nullptr),
      vel_fec_(nullptr),
      pres_fec_(nullptr),
      vel_fes_(nullptr),
      pres_fes_(nullptr),
      st_vel_fec_(nullptr),
      st_pres_fec_(nullptr),
      st_vel_fes_(nullptr),
      st_pres_fes_(nullptr),
      velocity_gf_(nullptr),
      pressure_gf_(nullptr),
      vorticity_gf_(nullptr),
      st_solution_(nullptr),
      newton_solver_(nullptr),
      linear_solver_(nullptr),
      preconditioner_(nullptr),
      body_force_(nullptr),
      paraview_dc_(nullptr),
      current_time_(config.t_start),
      slab_count_(0)
#ifdef MFEM_USE_MPI
      ,
      par_spatial_mesh_(nullptr),
      par_vel_fes_(nullptr),
      par_pres_fes_(nullptr)
#endif
{
    config_.ComputeDerivedQuantities();
}

// =========================================================================
// Constructor (parallel)
// =========================================================================

STNavierStokesSolver::STNavierStokesSolver(const SolverConfig &config,
                                           MPI_Comm comm)
    : config_(config),
      is_parallel_(true),
      comm_(comm),
      spatial_mesh_(nullptr),
      vel_fec_(nullptr),
      pres_fec_(nullptr),
      vel_fes_(nullptr),
      pres_fes_(nullptr),
      st_vel_fec_(nullptr),
      st_pres_fec_(nullptr),
      st_vel_fes_(nullptr),
      st_pres_fes_(nullptr),
      velocity_gf_(nullptr),
      pressure_gf_(nullptr),
      vorticity_gf_(nullptr),
      st_solution_(nullptr),
      newton_solver_(nullptr),
      linear_solver_(nullptr),
      preconditioner_(nullptr),
      body_force_(nullptr),
      paraview_dc_(nullptr),
      current_time_(config.t_start),
      slab_count_(0)
#ifdef MFEM_USE_MPI
      ,
      par_spatial_mesh_(nullptr),
      par_vel_fes_(nullptr),
      par_pres_fes_(nullptr)
#endif
{
    MPI_Comm_rank(comm_, &myid_);
    MPI_Comm_size(comm_, &num_procs_);
    config_.ComputeDerivedQuantities();
}

// =========================================================================
// Destructor
// =========================================================================

STNavierStokesSolver::~STNavierStokesSolver()
{
    delete paraview_dc_;
    delete st_solution_;
    delete newton_solver_;
    delete linear_solver_;
    delete preconditioner_;
    delete vorticity_gf_;
    delete pressure_gf_;
    delete velocity_gf_;
    delete st_pres_fes_;
    delete st_vel_fes_;
    delete st_pres_fec_;
    delete st_vel_fec_;
    delete pres_fes_;
    delete vel_fes_;
    delete pres_fec_;
    delete vel_fec_;
#ifdef MFEM_USE_MPI
    delete par_pres_fes_;
    delete par_vel_fes_;
    delete par_spatial_mesh_;
#endif
    delete spatial_mesh_;
}

// =========================================================================
// Initialize
// =========================================================================

void STNavierStokesSolver::Initialize()
{
    if (myid_ == 0)
    {
        std::cout << "========================================" << std::endl;
        std::cout << "StreamVorti Space-Time Navier-Stokes Solver" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Formulation: ";
        switch (config_.formulation)
        {
        case Formulation::ST_SUPG_PSPG:
            std::cout << "ST-SUPG/PSPG"; break;
        case Formulation::ST_VMS:
            std::cout << "ST-VMS"; break;
        case Formulation::ST_GLS:
            std::cout << "ST-GLS"; break;
        }
        std::cout << std::endl;
        std::cout << "Element type: "
                  << (config_.element_type == ElementType::Simplex
                      ? "Simplex" : "Tensor-Product")
                  << std::endl;
        std::cout << "Spatial dim:  " << config_.spatial_dim << std::endl;
        std::cout << "Re:           " << config_.reynolds << std::endl;
        std::cout << "Velocity order: " << config_.velocity_order << std::endl;
        std::cout << "Pressure order: " << config_.pressure_order << std::endl;
        std::cout << "dt_slab:      " << config_.dt_slab << std::endl;
        std::cout << "t_final:      " << config_.t_final << std::endl;
        std::cout << "Stabilization: SUPG="
                  << config_.stabilization.enable_supg
                  << " PSPG=" << config_.stabilization.enable_pspg
                  << " LSIC=" << config_.stabilization.enable_lsic
                  << std::endl;
        std::cout << "AMR:          "
                  << (config_.amr.strategy != AdaptiveStrategy::None
                      ? "ENABLED" : "DISABLED")
                  << std::endl;
        std::cout << "Parallel:     "
                  << (is_parallel_ ? "YES" : "NO");
        if (is_parallel_)
            std::cout << " (" << num_procs_ << " processes)";
        std::cout << std::endl;
        std::cout << "========================================" << std::endl;
    }

    LoadSpatialMesh();
    SetupFESpaces();
    SetupSolver();
    SetupOutput();
}

// =========================================================================
// LoadSpatialMesh
// =========================================================================

void STNavierStokesSolver::LoadSpatialMesh()
{
    if (config_.mesh_file.empty())
    {
        if (myid_ == 0)
        {
            std::cerr << "Error: No mesh file specified." << std::endl;
        }
        return;
    }

    // Load spatial mesh
    spatial_mesh_ = new mfem::Mesh(config_.mesh_file.c_str(), 1, 1);
    assert(spatial_mesh_->Dimension() == config_.spatial_dim);

    // Apply serial refinement
    for (int i = 0; i < config_.serial_refine; ++i)
    {
        spatial_mesh_->UniformRefinement();
    }

#ifdef MFEM_USE_MPI
    if (is_parallel_)
    {
        par_spatial_mesh_ = new mfem::ParMesh(comm_, *spatial_mesh_);

        // Apply parallel refinement
        for (int i = 0; i < config_.parallel_refine; ++i)
        {
            par_spatial_mesh_->UniformRefinement();
        }

        if (myid_ == 0)
        {
            std::cout << "Spatial mesh: "
                      << par_spatial_mesh_->GetGlobalNE() << " elements, "
                      << par_spatial_mesh_->GetNV() << " vertices (local)"
                      << std::endl;
        }
    }
    else
#endif
    {
        if (myid_ == 0)
        {
            std::cout << "Spatial mesh: "
                      << spatial_mesh_->GetNE() << " elements, "
                      << spatial_mesh_->GetNV() << " vertices"
                      << std::endl;
        }
    }
}

// =========================================================================
// SetupFESpaces
// =========================================================================

void STNavierStokesSolver::SetupFESpaces()
{
    int sdim = config_.spatial_dim;

    // FE collections for the spatial mesh (used for output/visualization)
    vel_fec_ = new mfem::H1_FECollection(config_.velocity_order, sdim);
    pres_fec_ = new mfem::H1_FECollection(config_.pressure_order, sdim);

#ifdef MFEM_USE_MPI
    if (is_parallel_ && par_spatial_mesh_)
    {
        par_vel_fes_ = new mfem::ParFiniteElementSpace(
            par_spatial_mesh_, vel_fec_, sdim);
        par_pres_fes_ = new mfem::ParFiniteElementSpace(
            par_spatial_mesh_, pres_fec_, 1);

        velocity_gf_ = new mfem::GridFunction(par_vel_fes_);
        pressure_gf_ = new mfem::GridFunction(par_pres_fes_);
    }
    else
#endif
    {
        vel_fes_ = new mfem::FiniteElementSpace(
            spatial_mesh_, vel_fec_, sdim);
        pres_fes_ = new mfem::FiniteElementSpace(
            spatial_mesh_, pres_fec_, 1);

        velocity_gf_ = new mfem::GridFunction(vel_fes_);
        pressure_gf_ = new mfem::GridFunction(pres_fes_);
    }

    *velocity_gf_ = 0.0;
    *pressure_gf_ = 0.0;

    // Vorticity for output (scalar in 2D, vector in 3D)
    if (sdim == 2)
    {
        mfem::FiniteElementSpace *vort_fes =
            new mfem::FiniteElementSpace(spatial_mesh_, pres_fec_, 1);
        vorticity_gf_ = new mfem::GridFunction(vort_fes);
        // Note: vort_fes ownership is tricky; in production code
        // we'd manage this more carefully
    }

    if (myid_ == 0)
    {
        std::cout << "Velocity DOFs: " << velocity_gf_->Size() << std::endl;
        std::cout << "Pressure DOFs: " << pressure_gf_->Size() << std::endl;
    }
}

// =========================================================================
// SetupSolver
// =========================================================================

void STNavierStokesSolver::SetupSolver()
{
    // The solver is set up per-slab since the space-time FE spaces
    // change with each slab (especially with AMR)
}

// =========================================================================
// SetupOutput
// =========================================================================

void STNavierStokesSolver::SetupOutput()
{
    mfem::Mesh *output_mesh = spatial_mesh_;
#ifdef MFEM_USE_MPI
    if (is_parallel_ && par_spatial_mesh_)
    {
        output_mesh = par_spatial_mesh_;
    }
#endif

    paraview_dc_ = new mfem::ParaViewDataCollection(
        config_.output.prefix, output_mesh);
    paraview_dc_->SetPrefixPath(config_.output.directory);
    paraview_dc_->SetDataFormat(mfem::VTKFormat::BINARY);
    paraview_dc_->SetHighOrderOutput(true);
    paraview_dc_->SetLevelsOfDetail(config_.velocity_order);

    if (config_.output.save_velocity)
    {
        paraview_dc_->RegisterField("Velocity", velocity_gf_);
    }
    if (config_.output.save_pressure)
    {
        paraview_dc_->RegisterField("Pressure", pressure_gf_);
    }
    if (config_.output.save_vorticity && vorticity_gf_)
    {
        paraview_dc_->RegisterField("Vorticity", vorticity_gf_);
    }
}

// =========================================================================
// SetInitialVelocity
// =========================================================================

void STNavierStokesSolver::SetInitialVelocity(mfem::VectorCoefficient &u0)
{
    velocity_gf_->ProjectCoefficient(u0);
    u_initial_.SetSize(velocity_gf_->Size());
    u_initial_ = *velocity_gf_;
}

// =========================================================================
// SetInitialPressure
// =========================================================================

void STNavierStokesSolver::SetInitialPressure(mfem::Coefficient &p0)
{
    pressure_gf_->ProjectCoefficient(p0);
    p_initial_.SetSize(pressure_gf_->Size());
    p_initial_ = *pressure_gf_;
}

// =========================================================================
// Solve (main time-slab marching loop)
// =========================================================================

void STNavierStokesSolver::Solve()
{
    if (myid_ == 0)
    {
        std::cout << "\nStarting space-time solve..." << std::endl;
    }

    // Initialize from initial conditions
    if (u_initial_.Size() == 0)
    {
        u_initial_.SetSize(velocity_gf_->Size());
        u_initial_ = 0.0;
    }

    // Save initial condition
    SaveOutput(0, config_.t_start);

    current_time_ = config_.t_start;
    slab_count_ = 0;

    // Time-slab marching loop
    while (current_time_ < config_.t_final - 1e-12)
    {
        double t_slab_end = std::min(current_time_ + config_.dt_slab,
                                     config_.t_final);

        slab_count_++;

        if (myid_ == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "\n--- Time slab " << slab_count_
                      << ": [" << current_time_ << ", "
                      << t_slab_end << "] ---" << std::endl;
        }

        // Solve this time slab
        SolveTimeSlab(current_time_, t_slab_end, u_initial_);

        // Update time
        current_time_ = t_slab_end;

        // Save output
        if (slab_count_ % config_.output.save_frequency == 0)
        {
            SaveOutput(slab_count_, current_time_);
        }

        // Adaptive mesh refinement
        if (config_.amr.strategy != AdaptiveStrategy::None &&
            slab_count_ % config_.amr.amr_frequency == 0)
        {
            mfem::Vector indicators;
            ComputeErrorIndicators(indicators);

            bool mesh_changed = AdaptSpatialMesh(indicators);
            if (mesh_changed && myid_ == 0)
            {
                std::cout << "  Mesh adapted. New element count: "
                          << spatial_mesh_->GetNE() << std::endl;
            }
        }
    }

    if (myid_ == 0)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Simulation complete." << std::endl;
        std::cout << "  Total time slabs: " << slab_count_ << std::endl;
        std::cout << "  Final time: " << current_time_ << std::endl;
        std::cout << "========================================" << std::endl;
    }
}

// =========================================================================
// SolveTimeSlab
// =========================================================================

void STNavierStokesSolver::SolveTimeSlab(double t_start, double t_end,
                                         const mfem::Vector &u_initial)
{
    // 1. Build space-time slab mesh
    BuildSpaceTimeSlab(t_start, t_end);

    // 2. Set up FE spaces on the slab
    SetupSlabFESpaces();

    // 3. Assemble and solve the nonlinear system
    AssembleSlabSystem();
    SolveSlabSystem();

    // 4. Extract solution at top of slab
    mfem::Vector u_top, p_top;
    ExtractTopSolution(u_top, p_top);

    // 5. Update spatial solution fields for output
    // Copy the extracted top-of-slab velocity back to the spatial grid function
    if (u_top.Size() == velocity_gf_->Size())
    {
        *velocity_gf_ = u_top;
    }
    if (p_top.Size() == pressure_gf_->Size())
    {
        *pressure_gf_ = p_top;
    }

    // 6. Compute vorticity for output (2D only)
    if (config_.spatial_dim == 2 && vorticity_gf_)
    {
        // ω = ∂v/∂x - ∂u/∂y
        // Use MFEM's CurlGridFunctionCoefficient or manual computation
        velocity_gf_->ComputeCurl2D(*vorticity_gf_);
    }

    // 7. Update initial condition for next slab
    u_initial_ = u_top;

    // Clean up slab-specific data
    delete st_solution_;
    st_solution_ = nullptr;
    ns_operator_.reset();
    delete st_pres_fes_; st_pres_fes_ = nullptr;
    delete st_vel_fes_;  st_vel_fes_ = nullptr;
    delete st_pres_fec_; st_pres_fec_ = nullptr;
    delete st_vel_fec_;  st_vel_fec_ = nullptr;
}

// =========================================================================
// BuildSpaceTimeSlab
// =========================================================================

void STNavierStokesSolver::BuildSpaceTimeSlab(double t_start, double t_end)
{
#ifdef MFEM_USE_MPI
    if (is_parallel_ && par_spatial_mesh_)
    {
        st_mesh_ = std::make_unique<SpaceTimeMesh>(
            *par_spatial_mesh_, t_start, t_end,
            config_.element_type, config_.n_time_steps_per_slab);
    }
    else
#endif
    {
        st_mesh_ = std::make_unique<SpaceTimeMesh>(
            *spatial_mesh_, t_start, t_end,
            config_.element_type, config_.n_time_steps_per_slab);
    }

    if (myid_ == 0)
    {
        std::cout << "  ST mesh: "
                  << st_mesh_->GetMesh()->GetNE() << " elements, "
                  << st_mesh_->GetMesh()->GetNV() << " vertices"
                  << std::endl;
    }
}

// =========================================================================
// SetupSlabFESpaces
// =========================================================================

void STNavierStokesSolver::SetupSlabFESpaces()
{
    int st_dim = config_.spatial_dim + 1;

    // For the space-time slab, we use standard H1 elements
    // in the (d+1)-dimensional domain
    st_vel_fec_ = new mfem::H1_FECollection(config_.velocity_order, st_dim);
    st_pres_fec_ = new mfem::H1_FECollection(config_.pressure_order, st_dim);

    mfem::Mesh *st_m = st_mesh_->GetMesh();

    // Velocity: vector-valued with spatial_dim components
    st_vel_fes_ = new mfem::FiniteElementSpace(
        st_m, st_vel_fec_, config_.spatial_dim);
    st_pres_fes_ = new mfem::FiniteElementSpace(
        st_m, st_pres_fec_, 1);

    if (myid_ == 0)
    {
        std::cout << "  ST velocity DOFs: " << st_vel_fes_->GetTrueVSize()
                  << std::endl;
        std::cout << "  ST pressure DOFs: " << st_pres_fes_->GetTrueVSize()
                  << std::endl;
    }
}

// =========================================================================
// AssembleSlabSystem
// =========================================================================

void STNavierStokesSolver::AssembleSlabSystem()
{
    ns_operator_ = std::make_unique<STNavierStokesOperator>(
        *st_vel_fes_, *st_pres_fes_, config_);

    if (body_force_)
    {
        ns_operator_->SetBodyForce(body_force_);
    }

    if (u_initial_.Size() > 0)
    {
        ns_operator_->SetPreviousSolution(u_initial_);
    }

    // Set up block structure for the solution
    const mfem::Array<int> &offsets = ns_operator_->GetBlockOffsets();
    st_solution_ = new mfem::BlockVector(offsets);
    *st_solution_ = 0.0;

    // Project initial condition onto the space-time slab
    ProjectInitialCondition(u_initial_);
}

// =========================================================================
// SolveSlabSystem
// =========================================================================

void STNavierStokesSolver::SolveSlabSystem()
{
    // Set up linear solver for the Jacobian system
    if (config_.linear_solver.type == "gmres")
    {
        mfem::GMRESSolver *gmres = new mfem::GMRESSolver();
        gmres->SetRelTol(config_.linear_solver.rel_tol);
        gmres->SetMaxIter(config_.linear_solver.max_iterations);
        gmres->SetKDim(config_.linear_solver.kdim);
        gmres->SetPrintLevel(config_.linear_solver.print_level);
        linear_solver_ = gmres;
    }
    else if (config_.linear_solver.type == "umfpack")
    {
#ifdef MFEM_USE_SUITESPARSE
        linear_solver_ = new mfem::UMFPackSolver();
#else
        // Fall back to GMRES
        mfem::GMRESSolver *gmres = new mfem::GMRESSolver();
        gmres->SetRelTol(1e-10);
        gmres->SetMaxIter(500);
        linear_solver_ = gmres;
#endif
    }
    else
    {
        // Default to GMRES
        mfem::GMRESSolver *gmres = new mfem::GMRESSolver();
        gmres->SetRelTol(config_.linear_solver.rel_tol);
        gmres->SetMaxIter(config_.linear_solver.max_iterations);
        linear_solver_ = gmres;
    }

    // Set up Newton solver
    newton_solver_ = new mfem::NewtonSolver();
    newton_solver_->SetSolver(*linear_solver_);
    newton_solver_->SetOperator(*ns_operator_);
    newton_solver_->SetRelTol(config_.nonlinear_solver.rel_tol);
    newton_solver_->SetAbsTol(config_.nonlinear_solver.abs_tol);
    newton_solver_->SetMaxIter(config_.nonlinear_solver.max_iterations);
    newton_solver_->SetPrintLevel(config_.nonlinear_solver.print_level);

    // Initial guess
    mfem::Vector zero(st_solution_->Size());
    zero = 0.0;

    // Solve
    newton_solver_->Mult(zero, *st_solution_);

    if (myid_ == 0)
    {
        if (newton_solver_->GetConverged())
        {
            std::cout << "  Newton converged in "
                      << newton_solver_->GetNumIterations()
                      << " iterations." << std::endl;
        }
        else
        {
            std::cout << "  WARNING: Newton did NOT converge after "
                      << newton_solver_->GetNumIterations()
                      << " iterations." << std::endl;
        }
    }

    // Clean up per-slab solver objects
    delete newton_solver_;
    newton_solver_ = nullptr;
    delete linear_solver_;
    linear_solver_ = nullptr;
}

// =========================================================================
// ProjectInitialCondition
// =========================================================================

void STNavierStokesSolver::ProjectInitialCondition(
    const mfem::Vector &u_prev)
{
    // Project the spatial velocity from the previous slab onto the
    // bottom face of the current space-time slab.
    //
    // For dG in time: this sets up the jump term [u]_{t_n} = u_n^+ - u_n^-
    // For cG in time: this sets Dirichlet BC at the bottom face

    if (u_prev.Size() == 0) return;

    // The bottom face DOFs correspond to the spatial mesh DOFs
    // at the first time layer. We need to identify which ST DOFs
    // map to the bottom face and set their values.

    mfem::BlockVector &sol = *st_solution_;
    mfem::Vector &u_sol = sol.GetBlock(0);

    // For simplicity, project the initial condition as a constant-in-time
    // field throughout the slab as the initial guess for Newton
    int n_spatial_vel_dofs = u_prev.Size();
    int n_st_vel_dofs = u_sol.Size();

    if (n_spatial_vel_dofs > 0 && n_st_vel_dofs > 0)
    {
        // Replicate spatial solution across time layers
        int n_time_layers = st_mesh_->GetNumTimeLayers();
        int dofs_per_layer = n_st_vel_dofs / n_time_layers;

        if (dofs_per_layer > 0 && dofs_per_layer <= n_spatial_vel_dofs)
        {
            for (int tl = 0; tl < n_time_layers; ++tl)
            {
                for (int i = 0; i < dofs_per_layer; ++i)
                {
                    if (i < n_spatial_vel_dofs)
                    {
                        u_sol(tl * dofs_per_layer + i) = u_prev(i);
                    }
                }
            }
        }
    }
}

// =========================================================================
// ExtractTopSolution
// =========================================================================

void STNavierStokesSolver::ExtractTopSolution(mfem::Vector &u_top,
                                               mfem::Vector &p_top)
{
    // Extract the solution at the top face (t = t_{n+1}) of the
    // space-time slab. This becomes the initial condition for the
    // next time slab.

    const mfem::BlockVector &sol = *st_solution_;
    const mfem::Vector &u_sol = sol.GetBlock(0);
    const mfem::Vector &p_sol = sol.GetBlock(1);

    int n_time_layers = st_mesh_->GetNumTimeLayers();
    int n_st_vel_dofs = u_sol.Size();
    int n_st_pres_dofs = p_sol.Size();

    int vel_dofs_per_layer = n_st_vel_dofs / n_time_layers;
    int pres_dofs_per_layer = n_st_pres_dofs / n_time_layers;

    // Extract top layer
    int top_layer = n_time_layers - 1;

    u_top.SetSize(vel_dofs_per_layer);
    for (int i = 0; i < vel_dofs_per_layer; ++i)
    {
        u_top(i) = u_sol(top_layer * vel_dofs_per_layer + i);
    }

    p_top.SetSize(pres_dofs_per_layer);
    for (int i = 0; i < pres_dofs_per_layer; ++i)
    {
        p_top(i) = p_sol(top_layer * pres_dofs_per_layer + i);
    }
}

// =========================================================================
// Adaptive mesh refinement
// =========================================================================

void STNavierStokesSolver::ComputeErrorIndicators(mfem::Vector &indicators)
{
    switch (config_.amr.indicator)
    {
    case ErrorIndicator::Residual:
        ComputeResidualIndicator(indicators);
        break;
    case ErrorIndicator::VelocityGradient:
        ComputeVelocityGradientIndicator(indicators);
        break;
    case ErrorIndicator::VorticityBased:
        ComputeVorticityIndicator(indicators);
        break;
    case ErrorIndicator::Kelly:
        ComputeVelocityGradientIndicator(indicators); // fallback
        break;
    case ErrorIndicator::ZZ_Recovery:
        ComputeVelocityGradientIndicator(indicators); // fallback
        break;
    }
}

void STNavierStokesSolver::ComputeResidualIndicator(mfem::Vector &indicators)
{
    int ne = spatial_mesh_->GetNE();
    indicators.SetSize(ne);

    // Element-wise residual error indicator
    // η_K² = h_K² ||R_M||²_K + h_K ||[∂u/∂n]||²_∂K
    //
    // where R_M = ρ(∂u/∂t + u·∇u) + ∇p - μ∇²u - f

    for (int i = 0; i < ne; ++i)
    {
        mfem::ElementTransformation *Trans =
            spatial_mesh_->GetElementTransformation(i);
        double h_K = spatial_mesh_->GetElementSize(i);

        // Compute velocity gradient norm on this element
        double grad_u_norm = 0.0;
        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(spatial_mesh_->GetElementGeometry(i),
                               2 * config_.velocity_order);

        for (int q = 0; q < ir.GetNPoints(); ++q)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);
            Trans->SetIntPoint(&ip);
            double w = ip.weight * Trans->Weight();

            // Evaluate velocity gradient
            mfem::DenseMatrix grad_u;
            velocity_gf_->GetVectorGradient(*Trans, grad_u);

            // Frobenius norm of gradient
            double norm_sq = 0.0;
            for (int r = 0; r < grad_u.Height(); ++r)
            {
                for (int c = 0; c < grad_u.Width(); ++c)
                {
                    norm_sq += grad_u(r, c) * grad_u(r, c);
                }
            }
            grad_u_norm += w * norm_sq;
        }

        indicators(i) = h_K * std::sqrt(grad_u_norm);
    }
}

void STNavierStokesSolver::ComputeVelocityGradientIndicator(
    mfem::Vector &indicators)
{
    int ne = spatial_mesh_->GetNE();
    indicators.SetSize(ne);

    for (int i = 0; i < ne; ++i)
    {
        mfem::ElementTransformation *Trans =
            spatial_mesh_->GetElementTransformation(i);

        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(spatial_mesh_->GetElementGeometry(i),
                               2 * config_.velocity_order);

        double indicator = 0.0;
        for (int q = 0; q < ir.GetNPoints(); ++q)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);
            Trans->SetIntPoint(&ip);
            double w = ip.weight * Trans->Weight();

            mfem::DenseMatrix grad_u;
            velocity_gf_->GetVectorGradient(*Trans, grad_u);

            double norm_sq = 0.0;
            for (int r = 0; r < grad_u.Height(); ++r)
            {
                for (int c = 0; c < grad_u.Width(); ++c)
                {
                    norm_sq += grad_u(r, c) * grad_u(r, c);
                }
            }
            indicator += w * norm_sq;
        }
        indicators(i) = std::sqrt(indicator);
    }
}

void STNavierStokesSolver::ComputeVorticityIndicator(
    mfem::Vector &indicators)
{
    int ne = spatial_mesh_->GetNE();
    indicators.SetSize(ne);

    for (int i = 0; i < ne; ++i)
    {
        mfem::ElementTransformation *Trans =
            spatial_mesh_->GetElementTransformation(i);

        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(spatial_mesh_->GetElementGeometry(i),
                               2 * config_.velocity_order);

        double indicator = 0.0;
        for (int q = 0; q < ir.GetNPoints(); ++q)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(q);
            Trans->SetIntPoint(&ip);
            double w = ip.weight * Trans->Weight();

            // Compute vorticity: ω = ∂v/∂x - ∂u/∂y (2D)
            mfem::DenseMatrix grad_u;
            velocity_gf_->GetVectorGradient(*Trans, grad_u);

            double vort = 0.0;
            if (config_.spatial_dim == 2)
            {
                vort = grad_u(1, 0) - grad_u(0, 1); // ∂v/∂x - ∂u/∂y
            }
            indicator += w * vort * vort;
        }
        indicators(i) = std::sqrt(indicator);
    }
}

bool STNavierStokesSolver::AdaptSpatialMesh(const mfem::Vector &indicators)
{
    if (config_.amr.strategy == AdaptiveStrategy::None)
        return false;

    int ne = spatial_mesh_->GetNE();
    if (ne >= config_.amr.max_elements)
        return false;

    // Determine refinement threshold using Dorfler marking
    double max_indicator = indicators.Max();
    double threshold = config_.amr.refine_fraction * max_indicator;

    mfem::Array<int> elements_to_refine;
    for (int i = 0; i < ne; ++i)
    {
        if (indicators(i) > threshold)
        {
            elements_to_refine.Append(i);
        }
    }

    if (elements_to_refine.Size() == 0)
        return false;

    // Check refinement level limit
    // (MFEM tracks this internally for NCMesh)
    spatial_mesh_->GeneralRefinement(elements_to_refine);

    // Rebuild FE spaces on the refined mesh
    // (The old grid functions are invalid after refinement)
    delete vel_fes_;
    delete pres_fes_;
    vel_fes_ = new mfem::FiniteElementSpace(
        spatial_mesh_, vel_fec_, config_.spatial_dim);
    pres_fes_ = new mfem::FiniteElementSpace(
        spatial_mesh_, pres_fec_, 1);

    // Create new grid functions and interpolate
    mfem::GridFunction *new_vel = new mfem::GridFunction(vel_fes_);
    mfem::GridFunction *new_pres = new mfem::GridFunction(pres_fes_);

    // Project current solution onto refined mesh
    // (MFEM handles the interpolation)
    new_vel->ProjectGridFunction(*velocity_gf_);
    new_pres->ProjectGridFunction(*pressure_gf_);

    delete velocity_gf_;
    delete pressure_gf_;
    velocity_gf_ = new_vel;
    pressure_gf_ = new_pres;

    // Update initial condition size
    u_initial_.SetSize(velocity_gf_->Size());
    u_initial_ = *velocity_gf_;

    // Update output
    SetupOutput();

    return true;
}

// =========================================================================
// Output
// =========================================================================

void STNavierStokesSolver::SaveOutput(int slab_index, double time)
{
    if (!paraview_dc_) return;

    paraview_dc_->SetCycle(slab_index);
    paraview_dc_->SetTime(time);
    paraview_dc_->Save();

    if (myid_ == 0)
    {
        std::cout << "  Output saved at t = " << time << std::endl;
    }
}

void STNavierStokesSolver::SaveMesh(int slab_index)
{
    if (!config_.output.save_mesh) return;

    std::string filename = config_.output.directory + "/" +
                          config_.output.prefix + "_mesh_" +
                          std::to_string(slab_index) + ".mesh";
    std::ofstream ofs(filename);
    spatial_mesh_->Print(ofs);
}

} // namespace SpaceTime
} // namespace StreamVorti

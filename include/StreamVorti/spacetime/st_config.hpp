/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
 *
 * Space-Time configuration and data structures.
 *
 * Based on the ST-VMS (Space-Time Variational Multiscale) formulation
 * of Behr, Tezduyar, and collaborators:
 *
 * References:
 *   [1] Tezduyar, Behr, Mittal, Liou (1992). "A new strategy for finite
 *       element computations involving moving boundaries and interfaces -
 *       The deforming-spatial-domain/stabilized-space-time procedure."
 *       Comput. Methods Appl. Mech. Engrg., 94(3), 339-351.
 *
 *   [2] Tezduyar (2003). "Computation of moving boundaries and interfaces
 *       and stabilization parameters." Int. J. Numer. Meth. Fluids, 43, 555-575.
 *
 *   [3] Karyofylli, Frings, Elgeti, Behr (2018). "Simplex space-time meshes
 *       in two-phase flow simulations." Int. J. Numer. Meth. Fluids, 86, 218-230.
 *
 *   [4] von Danwitz, Karyofylli, Hosters, Behr (2019). "Simplex space-time
 *       meshes in compressible flow simulations." Int. J. Numer. Meth. Fluids,
 *       91, 29-42.
 *
 *   [5] Behr (2001). "Stabilized space-time finite elements for
 *       unsteady fluid flows." Comput. Methods Appl. Mech. Engrg.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef STREAMVORTI_SPACETIME_CONFIG_HPP
#define STREAMVORTI_SPACETIME_CONFIG_HPP

#include <string>
#include <vector>
#include <functional>

namespace mfem { class Mesh; }

namespace StreamVorti {
namespace SpaceTime {

/// Space-time formulation type
enum class Formulation
{
    ST_SUPG_PSPG,   ///< Classical SUPG/PSPG stabilized (Tezduyar & Behr)
    ST_VMS,          ///< Variational Multiscale (Bazilevs, Calo, et al.)
    ST_GLS           ///< Galerkin/Least-Squares
};

/// Element type for the space-time slab
enum class ElementType
{
    Simplex,         ///< Tetrahedra (2D+t) or pentatopes (3D+t)
    TensorProduct    ///< Hexahedra (2D+t) or tesseracts (3D+t)
};

/// Temporal discretization within the slab
enum class TemporalScheme
{
    TimeContinuous,     ///< Continuous Galerkin in time (cG)
    TimeDiscontinuous   ///< Discontinuous Galerkin in time (dG)
};

/// Adaptive refinement strategy
enum class AdaptiveStrategy
{
    None,              ///< No adaptive refinement
    Spatial,           ///< Refine only in spatial dimensions
    Temporal,          ///< Refine only in temporal direction
    SpaceTime,         ///< Coupled space-time refinement
    Anisotropic        ///< Direction-dependent refinement
};

/// Error indicator type for AMR
enum class ErrorIndicator
{
    VelocityGradient,  ///< ||grad(u)||_L2 element-wise
    VorticityBased,    ///< Based on local vorticity magnitude
    Residual,          ///< Element residual of NS equations
    ZZ_Recovery,       ///< Zienkiewicz-Zhu recovery-based
    Kelly              ///< Kelly error estimator (face jumps)
};

/// Boundary condition type
enum class BCType
{
    Dirichlet,         ///< Prescribed velocity
    Neumann,           ///< Prescribed traction
    NoSlip,            ///< Zero velocity (u = 0)
    Slip,              ///< Free-slip (u·n = 0)
    Inflow,            ///< Prescribed inflow velocity profile
    Outflow,           ///< Stress-free outflow
    Pressure,          ///< Prescribed pressure
    Periodic           ///< Periodic boundary
};

/// Velocity function type: (x, y, [z], t) -> (vx, vy, [vz])
using VelocityFunction = std::function<void(const double *x, double t,
                                            double *v)>;

/// Boundary condition specification
struct BoundaryCondition
{
    int attribute;           ///< MFEM boundary attribute
    BCType type;             ///< Type of BC
    VelocityFunction vel_fn; ///< Velocity function (for Dirichlet/Inflow)
    double pressure;         ///< Pressure value (for Pressure BC)
    std::string name;        ///< Human-readable name
};

/// Stabilization parameters configuration
struct StabilizationParams
{
    bool enable_supg = true;   ///< SUPG stabilization
    bool enable_pspg = true;   ///< PSPG stabilization
    bool enable_lsic = true;   ///< LSIC (least-squares incompressibility)
    double c_1 = 4.0;          ///< Coefficient for tau_SUGN1
    double c_2 = 2.0;          ///< Coefficient for tau_SUGN2
    double c_3 = 4.0;          ///< Coefficient for tau_SUGN3 (LSIC)
};

/// Nonlinear solver parameters
struct NonlinearSolverParams
{
    int max_iterations = 30;
    double rel_tol = 1e-8;
    double abs_tol = 1e-12;
    int print_level = 1;
};

/// Linear solver parameters
struct LinearSolverParams
{
    std::string type = "gmres";  ///< gmres, cg, minres, umfpack, hypre
    int max_iterations = 500;
    double rel_tol = 1e-10;
    int kdim = 200;              ///< GMRES restart dimension
    int print_level = 0;
    std::string preconditioner = "ilu"; ///< ilu, amg, jacobi, gs
};

/// Adaptive mesh refinement parameters
struct AMRParams
{
    AdaptiveStrategy strategy = AdaptiveStrategy::None;
    ErrorIndicator indicator = ErrorIndicator::Residual;
    double refine_fraction = 0.3;     ///< Fraction of elements to refine
    double coarsen_fraction = 0.05;   ///< Fraction of elements to coarsen
    int max_refinement_level = 5;     ///< Maximum refinement depth
    int min_elements = 100;           ///< Minimum elements (prevent coarsening)
    int max_elements = 500000;        ///< Maximum elements (prevent explosion)
    int amr_frequency = 1;            ///< Re-mesh every N time slabs
    double temporal_error_tol = 1e-3; ///< Tolerance for temporal refinement
};

/// Output configuration
struct OutputParams
{
    std::string directory = "ParaView";
    std::string prefix = "spacetime_ns";
    int save_frequency = 1;        ///< Save every N time slabs
    bool save_velocity = true;
    bool save_pressure = true;
    bool save_vorticity = true;
    bool save_mesh = false;        ///< Save adapted mesh
    bool save_error = false;       ///< Save error indicator field
    std::string format = "vtk";    ///< vtk, mfem
};

/// Complete solver configuration
struct SolverConfig
{
    // Problem definition
    int spatial_dim = 2;                ///< Spatial dimension (2 or 3)
    double reynolds = 100.0;            ///< Reynolds number
    double density = 1.0;               ///< Fluid density
    double viscosity = 0.01;            ///< Kinematic viscosity (= 1/Re if nondim)

    // Space-time discretization
    Formulation formulation = Formulation::ST_VMS;
    ElementType element_type = ElementType::Simplex;
    TemporalScheme temporal_scheme = TemporalScheme::TimeDiscontinuous;
    int velocity_order = 2;             ///< Polynomial order for velocity
    int pressure_order = 1;             ///< Polynomial order for pressure

    // Time parameters
    double t_start = 0.0;
    double t_final = 10.0;
    double dt_slab = 0.1;              ///< Time-slab thickness
    int n_time_steps_per_slab = 1;     ///< Sub-divisions within slab

    // Mesh
    std::string mesh_file;              ///< Spatial mesh file
    mfem::Mesh *spatial_mesh = nullptr; ///< Pre-built mesh (e.g. from SDL); not owned
    int serial_refine = 0;              ///< Uniform refinements before partitioning
    int parallel_refine = 0;            ///< Uniform refinements after partitioning

    // Sub-configurations
    StabilizationParams stabilization;
    NonlinearSolverParams nonlinear_solver;
    LinearSolverParams linear_solver;
    AMRParams amr;
    OutputParams output;

    // Boundary conditions
    std::vector<BoundaryCondition> boundary_conditions;

    // Derived quantities
    void ComputeDerivedQuantities()
    {
        if (reynolds > 0.0)
            viscosity = 1.0 / reynolds;
    }
};

} // namespace SpaceTime
} // namespace StreamVorti

#endif // STREAMVORTI_SPACETIME_CONFIG_HPP

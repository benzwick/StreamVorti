/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2020-2026 Benjamin F. Zwick
 *
 * Main driver for the space-time velocity-pressure Navier-Stokes solver.
 *
 * This program solves the incompressible Navier-Stokes equations using
 * stabilized space-time finite elements (ST-VMS / SUPG-PSPG-LSIC) with
 * adaptive mesh refinement in both space and time.
 *
 * The solver is based on the methods of Tezduyar, Behr, and collaborators:
 *
 *   - DSD/SST (Deforming-Spatial-Domain / Stabilized Space-Time)
 *   - ST-VMS (Space-Time Variational Multiscale)
 *   - SUPG/PSPG/LSIC stabilization
 *
 * For 2D spatial problems, the space-time domain is 3D (x, y, t),
 * discretized with tetrahedra (simplex) or hexahedra (tensor-product).
 *
 * Usage:
 *   CLI mode:
 *     ./SpaceTimeNS -m mesh.mesh -Re 100 -dt 0.1 -tf 10.0 -o 2 -pv
 *
 *   SDL mode (requires ECL):
 *     ./SpaceTimeNS -f demo/st_cylinder.lisp -lp lisp -pv
 *
 * Contributors:
 *     Benjamin F. ZWICK
 */

#include <StreamVorti/spacetime/st_config.hpp>
#include <StreamVorti/spacetime/st_solver.hpp>

#ifdef STREAMVORTI_WITH_ECL
#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_loader.hpp>
#endif

#include "mfem.hpp"

#include <iostream>
#include <string>
#include <cstdlib>
#include <functional>

using namespace StreamVorti::SpaceTime;

// =========================================================================
// Boundary condition functions for standard benchmark problems
// =========================================================================

namespace BenchmarkBC {

/// Parabolic inflow profile for channel flow
/// u(y) = U_max * 4 * y * (H - y) / H^2
struct ParabolicInflow
{
    double U_max;    ///< Maximum centerline velocity
    double H;        ///< Channel height
    int flow_dir;    ///< Flow direction (0=x, 1=y)
    int normal_dir;  ///< Wall-normal direction

    void operator()(const double *x, double t, double *v) const
    {
        double y = x[normal_dir];
        double u_profile = U_max * 4.0 * y * (H - y) / (H * H);

        for (int d = 0; d < 3; ++d) v[d] = 0.0;
        v[flow_dir] = u_profile;
    }
};

/// Uniform inflow
struct UniformInflow
{
    double U;
    int flow_dir;

    void operator()(const double *x, double t, double *v) const
    {
        for (int d = 0; d < 3; ++d) v[d] = 0.0;
        v[flow_dir] = U;
    }
};

/// No-slip wall (zero velocity)
void NoSlip(const double *x, double t, double *v)
{
    v[0] = 0.0;
    v[1] = 0.0;
}

/// Lid-driven cavity top wall
void LidVelocity(const double *x, double t, double *v)
{
    v[0] = 1.0;
    v[1] = 0.0;
}

} // namespace BenchmarkBC

// =========================================================================
// VelocityCoefficient wrapper for boundary conditions
// =========================================================================

class VelocityBCCoefficient : public mfem::VectorCoefficient
{
public:
    VelocityBCCoefficient(int dim,
                          std::function<void(const double *, double, double *)> fn)
        : mfem::VectorCoefficient(dim), fn_(fn) {}

    void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override
    {
        double x[3] = {0.0, 0.0, 0.0};
        mfem::Vector transip(x, T.GetSpaceDim());
        T.Transform(ip, transip);
        fn_(x, 0.0, V.GetData());
    }

private:
    std::function<void(const double *, double, double *)> fn_;
};

// =========================================================================
// Parse command-line arguments
// =========================================================================

SolverConfig ParseCommandLine(int argc, char *argv[])
{
    SolverConfig config;

    // Default values
    const char *mesh_file = "";
    const char *sdl_file = "";
    const char *lisp_path = "";
    int dim = 2;
    double Re = 100.0;
    double dt_slab = 0.1;
    double t_final = 10.0;
    int vel_order = 2;
    int pres_order = 1;
    int n_time_steps = 1;
    int serial_ref = 0;
    int parallel_ref = 0;
    int element_type = 0;       // 0 = simplex, 1 = tensor-product
    int formulation = 1;        // 0 = SUPG/PSPG, 1 = VMS, 2 = GLS
    int temporal_scheme = 1;    // 0 = cG, 1 = dG
    bool enable_supg = true;
    bool enable_pspg = true;
    bool enable_lsic = true;
    bool enable_amr = false;
    int amr_indicator = 0;      // 0 = residual, 1 = grad, 2 = vort
    double amr_refine_frac = 0.3;
    int amr_max_level = 5;
    int amr_max_elems = 500000;
    int amr_freq = 1;
    bool paraview = false;
    const char *output_dir = "ParaView";
    const char *output_prefix = "spacetime_ns";
    int save_freq = 1;
    int newton_max_iter = 30;
    double newton_rel_tol = 1e-8;
    double newton_abs_tol = 1e-12;
    const char *linear_solver = "gmres";
    int linear_max_iter = 500;
    double linear_rel_tol = 1e-10;
    int gmres_kdim = 200;
    int print_level = 1;

    mfem::OptionsParser args(argc, argv);

    // Problem definition
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Spatial mesh file.");
    args.AddOption(&sdl_file, "-f", "--sdl-file",
                   "SDL simulation definition file (requires ECL).");
    args.AddOption(&lisp_path, "-lp", "--lisp-path",
                   "Path to SDL Lisp source files.");
    args.AddOption(&dim, "-dim", "--dimension",
                   "Spatial dimension (2 or 3).");
    args.AddOption(&Re, "-Re", "--reynolds",
                   "Reynolds number.");

    // Discretization
    args.AddOption(&vel_order, "-ov", "--velocity-order",
                   "Velocity finite element order.");
    args.AddOption(&pres_order, "-op", "--pressure-order",
                   "Pressure finite element order.");
    args.AddOption(&element_type, "-et", "--element-type",
                   "Element type: 0=simplex (tet), 1=tensor-product (hex).");
    args.AddOption(&formulation, "-form", "--formulation",
                   "Formulation: 0=SUPG/PSPG, 1=VMS, 2=GLS.");
    args.AddOption(&temporal_scheme, "-ts", "--temporal-scheme",
                   "Temporal scheme: 0=cG, 1=dG.");

    // Time parameters
    args.AddOption(&dt_slab, "-dt", "--dt-slab",
                   "Time slab thickness.");
    args.AddOption(&t_final, "-tf", "--t-final",
                   "Final simulation time.");
    args.AddOption(&n_time_steps, "-nts", "--time-steps-per-slab",
                   "Number of time steps per slab.");

    // Mesh refinement
    args.AddOption(&serial_ref, "-rs", "--serial-refine",
                   "Number of serial uniform refinements.");
    args.AddOption(&parallel_ref, "-rp", "--parallel-refine",
                   "Number of parallel uniform refinements.");

    // Stabilization
    args.AddOption(&enable_supg, "-supg", "--enable-supg",
                   "-no-supg", "--disable-supg",
                   "Enable SUPG stabilization.");
    args.AddOption(&enable_pspg, "-pspg", "--enable-pspg",
                   "-no-pspg", "--disable-pspg",
                   "Enable PSPG stabilization.");
    args.AddOption(&enable_lsic, "-lsic", "--enable-lsic",
                   "-no-lsic", "--disable-lsic",
                   "Enable LSIC stabilization.");

    // Adaptive mesh refinement
    args.AddOption(&enable_amr, "-amr", "--enable-amr",
                   "-no-amr", "--disable-amr",
                   "Enable adaptive mesh refinement.");
    args.AddOption(&amr_indicator, "-amri", "--amr-indicator",
                   "AMR error indicator: 0=residual, 1=grad, 2=vorticity.");
    args.AddOption(&amr_refine_frac, "-amrf", "--amr-refine-fraction",
                   "Fraction of elements to refine.");
    args.AddOption(&amr_max_level, "-amrl", "--amr-max-level",
                   "Maximum AMR refinement level.");
    args.AddOption(&amr_max_elems, "-amre", "--amr-max-elements",
                   "Maximum number of elements.");
    args.AddOption(&amr_freq, "-amrn", "--amr-frequency",
                   "Re-mesh every N time slabs.");

    // Nonlinear solver
    args.AddOption(&newton_max_iter, "-nmi", "--newton-max-iter",
                   "Maximum Newton iterations.");
    args.AddOption(&newton_rel_tol, "-nrt", "--newton-rel-tol",
                   "Newton relative tolerance.");
    args.AddOption(&newton_abs_tol, "-nat", "--newton-abs-tol",
                   "Newton absolute tolerance.");

    // Linear solver
    args.AddOption(&linear_solver, "-ls", "--linear-solver",
                   "Linear solver: gmres, cg, umfpack, hypre.");
    args.AddOption(&linear_max_iter, "-lmi", "--linear-max-iter",
                   "Maximum linear solver iterations.");
    args.AddOption(&linear_rel_tol, "-lrt", "--linear-rel-tol",
                   "Linear solver relative tolerance.");
    args.AddOption(&gmres_kdim, "-gk", "--gmres-kdim",
                   "GMRES restart dimension.");

    // Output
    args.AddOption(&paraview, "-pv", "--paraview",
                   "-no-pv", "--no-paraview",
                   "Enable ParaView output.");
    args.AddOption(&output_dir, "-od", "--output-dir",
                   "Output directory.");
    args.AddOption(&output_prefix, "-op_name", "--output-prefix",
                   "Output filename prefix.");
    args.AddOption(&save_freq, "-sf", "--save-frequency",
                   "Save every N time slabs.");
    args.AddOption(&print_level, "-pl", "--print-level",
                   "Print level (0=quiet, 1=normal, 2=verbose).");

    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        exit(1);
    }
    args.PrintOptions(std::cout);

    // Fill config from parsed options
    config.spatial_dim = dim;
    config.reynolds = Re;
    config.density = 1.0;
    config.mesh_file = mesh_file;

    config.velocity_order = vel_order;
    config.pressure_order = pres_order;
    config.element_type = (element_type == 0) ? ElementType::Simplex
                                              : ElementType::TensorProduct;
    switch (formulation)
    {
    case 0: config.formulation = Formulation::ST_SUPG_PSPG; break;
    case 1: config.formulation = Formulation::ST_VMS; break;
    case 2: config.formulation = Formulation::ST_GLS; break;
    }
    config.temporal_scheme = (temporal_scheme == 0)
                                 ? TemporalScheme::TimeContinuous
                                 : TemporalScheme::TimeDiscontinuous;

    config.dt_slab = dt_slab;
    config.t_final = t_final;
    config.n_time_steps_per_slab = n_time_steps;
    config.serial_refine = serial_ref;
    config.parallel_refine = parallel_ref;

    config.stabilization.enable_supg = enable_supg;
    config.stabilization.enable_pspg = enable_pspg;
    config.stabilization.enable_lsic = enable_lsic;

    if (enable_amr)
    {
        config.amr.strategy = AdaptiveStrategy::SpaceTime;
        switch (amr_indicator)
        {
        case 0: config.amr.indicator = ErrorIndicator::Residual; break;
        case 1: config.amr.indicator = ErrorIndicator::VelocityGradient; break;
        case 2: config.amr.indicator = ErrorIndicator::VorticityBased; break;
        }
        config.amr.refine_fraction = amr_refine_frac;
        config.amr.max_refinement_level = amr_max_level;
        config.amr.max_elements = amr_max_elems;
        config.amr.amr_frequency = amr_freq;
    }

    config.nonlinear_solver.max_iterations = newton_max_iter;
    config.nonlinear_solver.rel_tol = newton_rel_tol;
    config.nonlinear_solver.abs_tol = newton_abs_tol;
    config.nonlinear_solver.print_level = print_level;

    config.linear_solver.type = linear_solver;
    config.linear_solver.max_iterations = linear_max_iter;
    config.linear_solver.rel_tol = linear_rel_tol;
    config.linear_solver.kdim = gmres_kdim;
    config.linear_solver.print_level = (print_level > 1) ? 1 : 0;

    if (paraview)
    {
        config.output.directory = output_dir;
        config.output.prefix = output_prefix;
        config.output.save_frequency = save_freq;
        config.output.save_velocity = true;
        config.output.save_pressure = true;
        config.output.save_vorticity = true;
    }

    return config;
}

// =========================================================================
// Setup boundary conditions for known benchmark problems
// =========================================================================

void SetupCylinderBCs(SolverConfig &config, double H, double U_max)
{
    // Cylinder in channel benchmark (Schafer & Turek, 1996)
    //
    // Boundary attributes (from Gmsh):
    //   1 = inlet
    //   2 = outlet
    //   3 = top wall
    //   4 = bottom wall
    //   5 = cylinder surface

    BenchmarkBC::ParabolicInflow inflow{U_max, H, 0, 1};

    // Inlet: parabolic velocity profile
    BoundaryCondition bc_inlet;
    bc_inlet.attribute = 1;
    bc_inlet.type = BCType::Inflow;
    bc_inlet.vel_fn = [inflow](const double *x, double t, double *v)
    {
        inflow(x, t, v);
    };
    bc_inlet.name = "inlet";
    config.boundary_conditions.push_back(bc_inlet);

    // Outlet: stress-free
    BoundaryCondition bc_outlet;
    bc_outlet.attribute = 2;
    bc_outlet.type = BCType::Outflow;
    bc_outlet.name = "outlet";
    config.boundary_conditions.push_back(bc_outlet);

    // Top wall: no-slip
    BoundaryCondition bc_top;
    bc_top.attribute = 3;
    bc_top.type = BCType::NoSlip;
    bc_top.name = "top";
    config.boundary_conditions.push_back(bc_top);

    // Bottom wall: no-slip
    BoundaryCondition bc_bottom;
    bc_bottom.attribute = 4;
    bc_bottom.type = BCType::NoSlip;
    bc_bottom.name = "bottom";
    config.boundary_conditions.push_back(bc_bottom);

    // Cylinder: no-slip
    BoundaryCondition bc_cyl;
    bc_cyl.attribute = 5;
    bc_cyl.type = BCType::NoSlip;
    bc_cyl.name = "cylinder";
    config.boundary_conditions.push_back(bc_cyl);
}

void SetupCavityBCs(SolverConfig &config)
{
    // Lid-driven cavity
    // Attributes: 1=bottom, 2=right, 3=top(lid), 4=left

    BoundaryCondition bc_bottom;
    bc_bottom.attribute = 1;
    bc_bottom.type = BCType::NoSlip;
    bc_bottom.name = "bottom";
    config.boundary_conditions.push_back(bc_bottom);

    BoundaryCondition bc_right;
    bc_right.attribute = 2;
    bc_right.type = BCType::NoSlip;
    bc_right.name = "right";
    config.boundary_conditions.push_back(bc_right);

    BoundaryCondition bc_lid;
    bc_lid.attribute = 3;
    bc_lid.type = BCType::Dirichlet;
    bc_lid.vel_fn = BenchmarkBC::LidVelocity;
    bc_lid.name = "lid";
    config.boundary_conditions.push_back(bc_lid);

    BoundaryCondition bc_left;
    bc_left.attribute = 4;
    bc_left.type = BCType::NoSlip;
    bc_left.name = "left";
    config.boundary_conditions.push_back(bc_left);
}

// =========================================================================
// MAIN
// =========================================================================

int main(int argc, char *argv[])
{
    // Initialize MPI if available
    mfem::Mpi::Init(argc, argv);
    int myid = mfem::Mpi::WorldRank();
    int num_procs = mfem::Mpi::WorldSize();
    bool parallel = (num_procs > 1);

    // Parse command line
    SolverConfig config = ParseCommandLine(argc, argv);

    // Detect benchmark problem from mesh file name
    std::string mesh_name = config.mesh_file;
    bool is_cylinder = (mesh_name.find("cylinder") != std::string::npos);
    bool is_cavity = (mesh_name.find("cavity") != std::string::npos);

    if (is_cylinder)
    {
        double H = 0.41;      // Channel height (DFG benchmark)
        double U_max = 1.5;   // Max inflow velocity
        SetupCylinderBCs(config, H, U_max);
    }
    else if (is_cavity)
    {
        SetupCavityBCs(config);
    }

    // Create and run solver
    if (parallel)
    {
        STNavierStokesSolver solver(config, MPI_COMM_WORLD);
        solver.Initialize();
        solver.Solve();
    }
    else
    {
        STNavierStokesSolver solver(config);
        solver.Initialize();
        solver.Solve();
    }

    return EXIT_SUCCESS;
}

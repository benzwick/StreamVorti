/*
 * StreamVorti - Adaptive Space-Time Navier-Stokes Solver
 * Copyright (C) 2026 Benjamin F. Zwick
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
 * Problem definitions (boundary conditions, physics, mesh) are specified
 * via SDL files — there are no hardcoded benchmark setups.
 *
 * Usage:
 *   SDL mode (requires ECL):
 *     ./SpaceTimeNS -f demo/st_cavity.lisp -lp lisp -pv
 *     ./SpaceTimeNS -f demo/st_cylinder.lisp -lp lisp -pv
 *
 *   CLI overrides (combined with SDL file):
 *     ./SpaceTimeNS -f demo/st_cavity.lisp -lp lisp -Re 1000 -dt 0.05 -pv
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
// SDL to SolverConfig conversion
// =========================================================================

#ifdef STREAMVORTI_WITH_ECL

/// Map SDL BC type string to SpaceTime::BCType enum
static BCType SDLTypeToBC(const std::string &type)
{
    if (type == "no-slip")   return BCType::NoSlip;
    if (type == "velocity")  return BCType::Dirichlet;
    if (type == "inflow")    return BCType::Inflow;
    if (type == "outflow")   return BCType::Outflow;
    if (type == "pressure")  return BCType::Pressure;
    if (type == "slip")      return BCType::Slip;
    if (type == "neumann")   return BCType::Neumann;
    if (type == "periodic")  return BCType::Periodic;

    std::cerr << "Warning: unknown SDL BC type '" << type
              << "', defaulting to Dirichlet." << std::endl;
    return BCType::Dirichlet;
}

/// Build a VelocityFunction from SDL boundary condition's Lisp functions.
/// Returns a callable (x, t) -> v or nullptr if no velocity functions exist.
static VelocityFunction MakeVelocityFn(
    const StreamVorti::Lisp::BoundaryCondition &sdl_bc)
{
    // Capture raw pointers; the SDL config must outlive the solver run.
    auto *u_fn = sdl_bc.u_function.get();
    auto *v_fn = sdl_bc.v_function.get();
    auto *w_fn = sdl_bc.w_function.get();
    auto *scalar_fn = sdl_bc.function.get();

    if (u_fn || v_fn || w_fn)
    {
        return [u_fn, v_fn, w_fn](const double *x, double t, double *v)
        {
            v[0] = u_fn ? u_fn->evaluateAt(x[0], x[1], 0.0) : 0.0;
            v[1] = v_fn ? v_fn->evaluateAt(x[0], x[1], 0.0) : 0.0;
            v[2] = w_fn ? w_fn->evaluateAt(x[0], x[1], 0.0) : 0.0;
        };
    }
    if (scalar_fn)
    {
        return [scalar_fn](const double *x, double t, double *v)
        {
            v[0] = scalar_fn->evaluateAt(x[0], x[1], 0.0);
            v[1] = 0.0;
            v[2] = 0.0;
        };
    }
    return nullptr;
}

/// Populate SolverConfig from a loaded SDL SimulationConfig
static void LoadSDLConfig(const StreamVorti::Lisp::SimulationConfig &sdl,
                          SolverConfig &config)
{
    // Physics
    config.spatial_dim = sdl.dimension;
    config.reynolds = sdl.physics.reynolds;
    config.density = sdl.physics.density;

    // Mesh (SDL loader builds the mesh; pass it to the solver)
    if (sdl.mesh)
    {
        config.spatial_mesh = sdl.mesh.get();
    }

    // Boundary conditions
    config.boundary_conditions.clear();
    for (const auto &sdl_bc : sdl.boundaries)
    {
        BoundaryCondition bc;
        bc.attribute = sdl_bc.attribute;
        bc.type = SDLTypeToBC(sdl_bc.type);
        bc.name = sdl_bc.name;
        bc.vel_fn = MakeVelocityFn(sdl_bc);
        bc.pressure = 0.0;
        config.boundary_conditions.push_back(std::move(bc));
    }
}

#endif // STREAMVORTI_WITH_ECL

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

SolverConfig ParseCommandLine(int argc, char *argv[],
                              std::string &sdl_file_out,
                              std::string &lisp_path_out)
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

    sdl_file_out = sdl_file;
    lisp_path_out = lisp_path;

    return config;
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
    std::string sdl_file, lisp_path;
    SolverConfig config = ParseCommandLine(argc, argv, sdl_file, lisp_path);

    // Load problem definition from SDL file.
    // All boundary conditions, physics parameters, and mesh are defined
    // in SDL — there are no hardcoded problem setups.
#ifdef STREAMVORTI_WITH_ECL
    // Keep SDL config alive for the duration of the solver run, since
    // VelocityFunction lambdas capture pointers to its LispFunctions.
    StreamVorti::Lisp::SimulationConfig sdl_config;

    if (!sdl_file.empty())
    {
        StreamVorti::Lisp::Runtime::init(lisp_path);
        sdl_config = StreamVorti::Lisp::Loader::load(sdl_file);
        LoadSDLConfig(sdl_config, config);

        if (myid == 0)
        {
            std::cout << "Loaded SDL: " << sdl_config.name << std::endl;
            std::cout << "  Boundaries: "
                      << config.boundary_conditions.size() << std::endl;
        }
    }
#endif

    if (config.boundary_conditions.empty())
    {
        if (myid == 0)
        {
            std::cerr << "Error: no boundary conditions defined.\n"
                      << "Use an SDL file to define the problem:\n"
                      << "  ./SpaceTimeNS -f demo/st_cavity.lisp -lp lisp -pv\n"
                      << "  ./SpaceTimeNS -f demo/st_cylinder.lisp -lp lisp -pv\n";
#ifndef STREAMVORTI_WITH_ECL
            std::cerr << "\nNote: rebuild with -DSTREAMVORTI_WITH_ECL=ON "
                      << "to enable SDL support." << std::endl;
#endif
        }
        return EXIT_FAILURE;
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

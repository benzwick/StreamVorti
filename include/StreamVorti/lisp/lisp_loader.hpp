/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
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
 */

/**
 * @file lisp_loader.hpp
 * @brief SDL (Simulation Definition Language) file loader
 *
 * Provides functionality to load and execute .lisp SDL files,
 * extracting simulation configuration for StreamVorti.
 */

#ifndef STREAMVORTI_LISP_LISP_LOADER_HPP_
#define STREAMVORTI_LISP_LISP_LOADER_HPP_

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <functional>

// Include lisp_bridge.hpp for LispFunction definition
// (required because BoundaryCondition uses unique_ptr<LispFunction>
// which needs complete type for destructor)
#include "StreamVorti/lisp/lisp_bridge.hpp"

// Forward declarations
namespace mfem {
    class Mesh;
    class GridFunction;
}

namespace StreamVorti {
namespace Lisp {

/**
 * @struct BoundaryCondition
 * @brief Represents a boundary condition loaded from SDL
 */
struct BoundaryCondition {
    std::string name;           ///< Region name
    int attribute;              ///< Mesh boundary attribute
    std::string type;           ///< "velocity", "pressure", etc.

    /// For scalar BCs (pressure, temperature, or legacy single-function velocity)
    std::unique_ptr<LispFunction> function;

    /// For vector velocity BCs (separate u, v, w components)
    std::unique_ptr<LispFunction> u_function;  ///< u-velocity (x-direction)
    std::unique_ptr<LispFunction> v_function;  ///< v-velocity (y-direction)
    std::unique_ptr<LispFunction> w_function;  ///< w-velocity (z-direction, 3D only)

    BoundaryCondition() = default;
    BoundaryCondition(BoundaryCondition&&) = default;
    BoundaryCondition& operator=(BoundaryCondition&&) = default;
};

/**
 * @struct DCPSEParams
 * @brief DCPSE discretization parameters from SDL
 */
struct DCPSEParams {
    int num_neighbors = 25;     ///< Number of support domain neighbors
    double cutoff_radius = 30;  ///< Cutoff radius multiplier
    double support_radius = 5;  ///< Support radius multiplier
};

/**
 * @struct SolverParams
 * @brief Time integration and solver parameters from SDL
 */
struct SolverParams {
    std::string timestepping = "explicit-euler";
    double dt = 0.001;
    double end_time = 1.0;
    double tolerance = 1e-6;
    int max_iterations = 10000;
};

/**
 * @struct PhysicsParams
 * @brief Physics parameters from SDL
 */
struct PhysicsParams {
    std::string type = "incompressible-navier-stokes";
    std::string formulation = "stream-vorticity";
    double reynolds = 100.0;
    double density = 1.0;
    double viscosity = 0.01;
};

/**
 * @struct OutputParams
 * @brief Output configuration from SDL
 */
struct OutputParams {
    std::string format = "vtk";
    double interval = 0.1;
    std::string directory = "results/";
    std::vector<std::string> fields;
};

/**
 * @struct SimulationConfig
 * @brief Complete simulation configuration loaded from SDL
 */
struct SimulationConfig {
    std::string name;
    int version = 1;
    int dimension = 2;

    // Mesh (owned)
    std::unique_ptr<mfem::Mesh> mesh;

    // Boundary conditions
    std::vector<BoundaryCondition> boundaries;

    // Parameters
    PhysicsParams physics;
    DCPSEParams dcpse;
    SolverParams solver;
    OutputParams output;

    SimulationConfig() = default;
    SimulationConfig(SimulationConfig&&) = default;
    SimulationConfig& operator=(SimulationConfig&&) = default;
};

/**
 * @class Loader
 * @brief Loads SDL files and extracts simulation configuration
 *
 * Usage:
 * @code
 * auto config = StreamVorti::Lisp::Loader::load("simulation.lisp");
 * // Use config.mesh, config.boundaries, config.dcpse, etc.
 * @endcode
 */
class Loader {
public:
    /**
     * @brief Load simulation from SDL file
     *
     * @param path Path to .lisp SDL file
     * @return SimulationConfig populated from file
     * @throws EclException on parse or evaluation error
     */
    static SimulationConfig load(const std::string& path);

    /**
     * @brief Load simulation from SDL string
     *
     * @param sdl_content SDL content as string
     * @return SimulationConfig populated from content
     * @throws EclException on parse or evaluation error
     */
    static SimulationConfig loadFromString(const std::string& sdl_content);

    /**
     * @brief Get a Lisp function by name from the loaded SDL
     *
     * @param name Function name (in current package)
     * @return LispFunction wrapper, or invalid if not found
     */
    static LispFunction getFunction(const std::string& name);

    /**
     * @brief Extract mesh from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return MFEM Mesh pointer (caller owns)
     */
    static std::unique_ptr<mfem::Mesh> extractMesh(EclObject sim_obj);

    /**
     * @brief Extract boundary conditions from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return Vector of boundary conditions
     */
    static std::vector<BoundaryCondition> extractBoundaries(EclObject sim_obj);

    /**
     * @brief Extract DCPSE parameters from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return DCPSEParams struct
     */
    static DCPSEParams extractDCPSEParams(EclObject sim_obj);

    /**
     * @brief Extract solver parameters from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return SolverParams struct
     */
    static SolverParams extractSolverParams(EclObject sim_obj);

    /**
     * @brief Extract physics parameters from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return PhysicsParams struct
     */
    static PhysicsParams extractPhysicsParams(EclObject sim_obj);

    /**
     * @brief Extract output parameters from loaded SDL
     *
     * @param sim_obj The Lisp simulation object
     * @return OutputParams struct
     */
    static OutputParams extractOutputParams(EclObject sim_obj);

private:
    /**
     * @brief Helper to get property from Lisp plist
     */
    static EclObject getProperty(EclObject plist, const std::string& key);

    /**
     * @brief Helper to get integer property
     */
    static int getIntProperty(EclObject plist, const std::string& key, int default_value);

    /**
     * @brief Helper to get double property
     */
    static double getDoubleProperty(EclObject plist, const std::string& key, double default_value);

    /**
     * @brief Helper to get string property
     */
    static std::string getStringProperty(EclObject plist, const std::string& key,
                                          const std::string& default_value);
};

/**
 * @class SimulationRunner
 * @brief Runs a simulation from loaded SDL configuration
 */
class SimulationRunner {
public:
    /**
     * @brief Create runner from configuration
     */
    explicit SimulationRunner(SimulationConfig&& config);

    /**
     * @brief Initialize the simulation (create DCPSE, apply BCs)
     */
    void initialize();

    /**
     * @brief Run the simulation
     */
    void run();

    /**
     * @brief Get current simulation time
     */
    double currentTime() const { return current_time_; }

    /**
     * @brief Get reference to configuration
     */
    const SimulationConfig& config() const { return config_; }

private:
    SimulationConfig config_;
    double current_time_ = 0.0;
    int step_ = 0;
};

} // namespace Lisp
} // namespace StreamVorti

#endif // STREAMVORTI_LISP_LISP_LOADER_HPP_

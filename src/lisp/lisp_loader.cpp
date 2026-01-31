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
 * @file lisp_loader.cpp
 * @brief Implementation of SDL file loader
 */

// IMPORTANT: Include MFEM FIRST before ECL to avoid macro conflicts
// ECL defines macros (like Ct) that conflict with MFEM member names
#include "mfem.hpp"

// Now include ECL
#include <ecl/ecl.h>

// Include our headers after MFEM and ECL
#include "StreamVorti/lisp/lisp_loader.hpp"
#include "StreamVorti/lisp/ecl_runtime.hpp"
#include "StreamVorti/lisp/lisp_bridge.hpp"
#include "StreamVorti/lisp/lisp_mesh.hpp"

#include <iostream>
#include <sstream>
#include <filesystem>

namespace StreamVorti {
namespace Lisp {

// Helper functions for type conversion
static inline cl_object to_cl(EclObject x) { return reinterpret_cast<cl_object>(x); }
static inline EclObject to_ecl(cl_object x) { return reinterpret_cast<EclObject>(x); }

SimulationConfig Loader::load(const std::string& path)
{
    // Ensure runtime is initialized
    if (!Runtime::isInitialized()) {
        throw EclException("ECL runtime not initialized. Call Runtime::init() first.");
    }

    // Check file exists
    if (!std::filesystem::exists(path)) {
        throw EclException("SDL file not found: " + path);
    }

    // Load and evaluate the SDL file
    Runtime::loadFile(path);

    // Get the simulation object from Lisp
    // The SDL file should define *current-simulation* or return a simulation object
    EclObject sim_obj = Runtime::eval("(sdl:get-current-simulation)");

    if (Bridge::isNil(sim_obj)) {
        throw EclException("No simulation defined in SDL file: " + path);
    }

    // Extract configuration from the Lisp object
    SimulationConfig config;

    // Get basic properties
    config.name = getStringProperty(sim_obj, "name", "unnamed");
    config.version = getIntProperty(sim_obj, "version", 1);
    config.dimension = getIntProperty(sim_obj, "dimension", 2);

    // Extract mesh
    config.mesh = extractMesh(sim_obj);

    // Extract boundary conditions
    config.boundaries = extractBoundaries(sim_obj);

    // Extract parameters
    config.physics = extractPhysicsParams(sim_obj);
    config.dcpse = extractDCPSEParams(sim_obj);
    config.solver = extractSolverParams(sim_obj);
    config.output = extractOutputParams(sim_obj);

    return config;
}

SimulationConfig Loader::loadFromString(const std::string& sdl_content)
{
    if (!Runtime::isInitialized()) {
        throw EclException("ECL runtime not initialized. Call Runtime::init() first.");
    }

    // Evaluate the SDL content directly
    std::ostringstream eval_expr;
    eval_expr << "(progn " << sdl_content << ")";

    Runtime::eval(eval_expr.str());

    // Get the simulation object
    EclObject sim_obj = Runtime::eval("(sdl:get-current-simulation)");

    if (Bridge::isNil(sim_obj)) {
        throw EclException("No simulation defined in SDL content");
    }

    SimulationConfig config;
    config.name = getStringProperty(sim_obj, "name", "unnamed");
    config.version = getIntProperty(sim_obj, "version", 1);
    config.dimension = getIntProperty(sim_obj, "dimension", 2);
    config.mesh = extractMesh(sim_obj);
    config.boundaries = extractBoundaries(sim_obj);
    config.physics = extractPhysicsParams(sim_obj);
    config.dcpse = extractDCPSEParams(sim_obj);
    config.solver = extractSolverParams(sim_obj);
    config.output = extractOutputParams(sim_obj);

    return config;
}

LispFunction Loader::getFunction(const std::string& name)
{
    return LispFunction(name, "SDL");
}

std::unique_ptr<mfem::Mesh> Loader::extractMesh(EclObject sim_obj)
{
    // Get mesh specification from simulation object
    EclObject mesh_spec = getProperty(sim_obj, "mesh");

    if (Bridge::isNil(mesh_spec)) {
        throw EclException("No mesh specification in simulation");
    }

    // Check if it's a generated mesh or loaded mesh
    EclObject mesh_type = getProperty(mesh_spec, "type");
    std::string type_str = Bridge::toCppString(mesh_type);

    if (type_str == "generated" || type_str == "generate") {
        // Generate mesh from parameters
        int dim = getIntProperty(sim_obj, "dimension", 2);
        EclObject divisions = getProperty(mesh_spec, "divisions");

        std::vector<int> divs = Bridge::toIntVector(divisions);

        if (dim == 2) {
            int nx = divs.size() > 0 ? divs[0] : 10;
            int ny = divs.size() > 1 ? divs[1] : 10;
            double sx = getDoubleProperty(mesh_spec, "size-x", 1.0);
            double sy = getDoubleProperty(mesh_spec, "size-y", 1.0);

            std::string elem_type_str = getStringProperty(mesh_spec, "element-type", "quad");
            int elem_type = (elem_type_str == "tri" || elem_type_str == "triangle") ? 2 : 3;

            return MeshWrapper::makeCartesian2D(nx, ny, elem_type, sx, sy);
        }
        else if (dim == 3) {
            int nx = divs.size() > 0 ? divs[0] : 10;
            int ny = divs.size() > 1 ? divs[1] : 10;
            int nz = divs.size() > 2 ? divs[2] : 10;
            double sx = getDoubleProperty(mesh_spec, "size-x", 1.0);
            double sy = getDoubleProperty(mesh_spec, "size-y", 1.0);
            double sz = getDoubleProperty(mesh_spec, "size-z", 1.0);

            std::string elem_type_str = getStringProperty(mesh_spec, "element-type", "hex");
            int elem_type = (elem_type_str == "tet" || elem_type_str == "tetrahedron") ? 4 : 5;

            return MeshWrapper::makeCartesian3D(nx, ny, nz, elem_type, sx, sy, sz);
        }
    }
    else if (type_str == "loaded" || type_str == "load") {
        // Load mesh from file
        std::string mesh_path = getStringProperty(mesh_spec, "path", "");
        if (mesh_path.empty()) {
            throw EclException("Mesh load path not specified");
        }
        return MeshWrapper::loadFromFile(mesh_path);
    }

    throw EclException("Unknown mesh type: " + type_str);
}

std::vector<BoundaryCondition> Loader::extractBoundaries(EclObject sim_obj)
{
    std::vector<BoundaryCondition> boundaries;

    EclObject bc_list = getProperty(sim_obj, "boundaries");
    if (Bridge::isNil(bc_list)) {
        return boundaries;
    }

    // Iterate over boundary condition list using cl_cdr
    cl_object cl_bc_list = to_cl(bc_list);
    while (cl_bc_list != ECL_NIL && !Bridge::isNil(to_ecl(cl_consp(cl_bc_list)))) {
        EclObject bc_spec = to_ecl(cl_car(cl_bc_list));

        BoundaryCondition bc;
        bc.name = getStringProperty(bc_spec, "name", "");
        bc.attribute = getIntProperty(bc_spec, "attribute", 1);
        bc.type = getStringProperty(bc_spec, "type", "velocity");

        // Get the boundary function(s)
        if (bc.type == "velocity") {
            // Try to get separate u and v functions for vector velocity BCs
            EclObject u_func = getProperty(bc_spec, "u-function");
            EclObject v_func = getProperty(bc_spec, "v-function");
            EclObject w_func = getProperty(bc_spec, "w-function");

            if (!Bridge::isNil(u_func) && Bridge::isFunction(u_func)) {
                bc.u_function = std::make_unique<LispFunction>(u_func);
            }
            if (!Bridge::isNil(v_func) && Bridge::isFunction(v_func)) {
                bc.v_function = std::make_unique<LispFunction>(v_func);
            }
            if (!Bridge::isNil(w_func) && Bridge::isFunction(w_func)) {
                bc.w_function = std::make_unique<LispFunction>(w_func);
            }

            // Fallback: if single function provided, use it for u-component (legacy support)
            EclObject func = getProperty(bc_spec, "function");
            if (!Bridge::isNil(func) && Bridge::isFunction(func)) {
                bc.function = std::make_unique<LispFunction>(func);
                // If no separate u-function, use the legacy function for u
                if (!bc.u_function) {
                    EclObject func_copy = getProperty(bc_spec, "function");
                    if (!Bridge::isNil(func_copy) && Bridge::isFunction(func_copy)) {
                        bc.u_function = std::make_unique<LispFunction>(func_copy);
                    }
                }
            }
        } else {
            // For non-velocity BCs (pressure, temperature), use single function
            EclObject func = getProperty(bc_spec, "function");
            if (!Bridge::isNil(func) && Bridge::isFunction(func)) {
                bc.function = std::make_unique<LispFunction>(func);
            }
        }

        boundaries.push_back(std::move(bc));

        // Move to next in list
        cl_bc_list = cl_cdr(cl_bc_list);
    }

    return boundaries;
}

DCPSEParams Loader::extractDCPSEParams(EclObject sim_obj)
{
    DCPSEParams params;

    EclObject dcpse_spec = getProperty(sim_obj, "discretization");
    if (Bridge::isNil(dcpse_spec)) {
        return params;
    }

    params.num_neighbors = getIntProperty(dcpse_spec, "num-neighbors", 25);
    params.cutoff_radius = getDoubleProperty(dcpse_spec, "cutoff-radius", 30.0);
    params.support_radius = getDoubleProperty(dcpse_spec, "support-radius", 5.0);

    return params;
}

SolverParams Loader::extractSolverParams(EclObject sim_obj)
{
    SolverParams params;

    EclObject solver_spec = getProperty(sim_obj, "solver");
    if (Bridge::isNil(solver_spec)) {
        return params;
    }

    params.timestepping = getStringProperty(solver_spec, "timestepping", "explicit-euler");
    params.dt = getDoubleProperty(solver_spec, "dt", 0.001);
    params.end_time = getDoubleProperty(solver_spec, "end-time", 1.0);
    params.tolerance = getDoubleProperty(solver_spec, "tolerance", 1e-6);
    params.max_iterations = getIntProperty(solver_spec, "max-iterations", 10000);

    return params;
}

PhysicsParams Loader::extractPhysicsParams(EclObject sim_obj)
{
    PhysicsParams params;

    EclObject physics_spec = getProperty(sim_obj, "physics");
    if (Bridge::isNil(physics_spec)) {
        return params;
    }

    params.type = getStringProperty(physics_spec, "type", "incompressible-navier-stokes");
    params.formulation = getStringProperty(physics_spec, "formulation", "stream-vorticity");
    params.reynolds = getDoubleProperty(physics_spec, "reynolds", 100.0);
    params.density = getDoubleProperty(physics_spec, "density", 1.0);
    params.viscosity = getDoubleProperty(physics_spec, "viscosity", 0.01);

    return params;
}

OutputParams Loader::extractOutputParams(EclObject sim_obj)
{
    OutputParams params;

    EclObject output_spec = getProperty(sim_obj, "output");
    if (Bridge::isNil(output_spec)) {
        return params;
    }

    params.format = getStringProperty(output_spec, "format", "vtk");
    params.interval = getDoubleProperty(output_spec, "interval", 0.1);
    params.directory = getStringProperty(output_spec, "directory", "results/");

    // Get fields list
    EclObject fields = getProperty(output_spec, "fields");
    if (!Bridge::isNil(fields) && Bridge::isList(fields)) {
        cl_object cl_fields = to_cl(fields);
        while (cl_fields != ECL_NIL && !Bridge::isNil(to_ecl(cl_consp(cl_fields)))) {
            EclObject field = to_ecl(cl_car(cl_fields));
            params.fields.push_back(Bridge::toCppString(field));
            cl_fields = cl_cdr(cl_fields);
        }
    }

    return params;
}

EclObject Loader::getProperty(EclObject plist, const std::string& key)
{
    // Use Lisp getf to get property from plist
    std::ostringstream expr;
    expr << "(getf '" << "(simulation-data) :" << key << ")";

    // For now, use a simpler approach - call the accessor function
    std::string accessor = "sdl:get-" + key;
    try {
        return Bridge::funcall(accessor, {plist});
    } catch (...) {
        return to_ecl(ECL_NIL);
    }
}

int Loader::getIntProperty(EclObject plist, const std::string& key, int default_value)
{
    EclObject val = getProperty(plist, key);
    if (Bridge::isNil(val) || !Bridge::isNumber(val)) {
        return default_value;
    }
    return Bridge::toCppInt(val);
}

double Loader::getDoubleProperty(EclObject plist, const std::string& key, double default_value)
{
    EclObject val = getProperty(plist, key);
    if (Bridge::isNil(val) || !Bridge::isNumber(val)) {
        return default_value;
    }
    return Bridge::toCppDouble(val);
}

std::string Loader::getStringProperty(EclObject plist, const std::string& key,
                                       const std::string& default_value)
{
    EclObject val = getProperty(plist, key);
    if (Bridge::isNil(val)) {
        return default_value;
    }
    return Bridge::toCppString(val);
}

// ==================== SimulationRunner ====================

SimulationRunner::SimulationRunner(SimulationConfig&& config)
    : config_(std::move(config))
{
}

void SimulationRunner::initialize()
{
    std::cout << "Initializing simulation: " << config_.name << std::endl;
    std::cout << "  Dimension: " << config_.dimension << std::endl;
    std::cout << "  Mesh vertices: " << config_.mesh->GetNV() << std::endl;
    std::cout << "  Boundaries: " << config_.boundaries.size() << std::endl;
    std::cout << "  DCPSE neighbors: " << config_.dcpse.num_neighbors << std::endl;

    // TODO: Create DCPSE operator, apply boundary conditions, etc.
}

void SimulationRunner::run()
{
    std::cout << "Running simulation..." << std::endl;
    std::cout << "  dt: " << config_.solver.dt << std::endl;
    std::cout << "  end_time: " << config_.solver.end_time << std::endl;

    // Time stepping loop (placeholder)
    while (current_time_ < config_.solver.end_time) {
        // TODO: Implement actual time stepping

        current_time_ += config_.solver.dt;
        step_++;

        // Output at intervals
        if (step_ % 100 == 0) {
            std::cout << "  Step " << step_ << ", t = " << current_time_ << std::endl;
        }
    }

    std::cout << "Simulation complete." << std::endl;
}

} // namespace Lisp
} // namespace StreamVorti

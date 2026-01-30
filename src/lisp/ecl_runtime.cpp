/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
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
 * @file ecl_runtime.cpp
 * @brief Implementation of ECL runtime management
 */

// Include ECL header
#include <ecl/ecl.h>

// Include our header after ECL
#include "StreamVorti/lisp/ecl_runtime.hpp"

#include <iostream>
#include <sstream>
#include <filesystem>

namespace StreamVorti {
namespace Lisp {

// Helper function for type conversion
static inline EclObject to_ecl(cl_object x) { return reinterpret_cast<EclObject>(x); }

// Static member definitions
bool Runtime::initialized_ = false;
std::string Runtime::last_error_;
std::string Runtime::lisp_path_;

void Runtime::init(const std::string& lisp_path)
{
    if (initialized_) {
        return; // Already initialized
    }

    // Store lisp path
    lisp_path_ = lisp_path;

    // Boot ECL with minimal arguments
    char* argv[] = { const_cast<char*>("streamvorti"), nullptr };
    int argc = 1;

    cl_boot(argc, argv);

    // Set up error handling
    // ECL uses conditions/restarts, we'll catch and convert to exceptions

    initialized_ = true;
    clearError();

    // If lisp_path provided, load the SDL packages
    if (!lisp_path_.empty()) {
        std::filesystem::path base_path(lisp_path_);

        // Load core SDL files in order
        std::vector<std::string> core_files = {
            "packages.lisp",
            "sdl-macros.lisp",
            "geometry.lisp",
            "mesh.lisp",
            "boundaries.lisp",
            "simulation.lisp"
        };

        for (const auto& file : core_files) {
            std::filesystem::path file_path = base_path / file;
            if (std::filesystem::exists(file_path)) {
                try {
                    loadFile(file_path.string());
                } catch (const EclException& e) {
                    std::cerr << "Warning: Failed to load " << file
                              << ": " << e.what() << std::endl;
                }
            }
        }
    }
}

void Runtime::shutdown()
{
    if (!initialized_) {
        return;
    }

    cl_shutdown();
    initialized_ = false;
}

bool Runtime::isInitialized()
{
    return initialized_;
}

EclObject Runtime::eval(const std::string& expr)
{
    if (!initialized_) {
        throw EclException("ECL runtime not initialized");
    }

    clearError();

    // Convert string to Lisp object and evaluate
    cl_object form = c_string_to_object(expr.c_str());

    // Use cl_safe_eval for error handling
    cl_object result = cl_safe_eval(form, Cnil, Cnil);

    // Check for NIL result which might indicate error
    // (actual error checking would need condition handling)
    if (result == OBJNULL) {
        last_error_ = "Evaluation returned null object";
        throw EclException(last_error_);
    }

    return to_ecl(result);
}

EclObject Runtime::safeEval(const std::string& expr)
{
    try {
        return eval(expr);
    } catch (const EclException& e) {
        last_error_ = e.what();
        return nullptr;
    }
}

void Runtime::loadFile(const std::string& path)
{
    if (!initialized_) {
        throw EclException("ECL runtime not initialized");
    }

    if (!std::filesystem::exists(path)) {
        throw EclException("File not found: " + path);
    }

    clearError();

    // Build load expression
    std::ostringstream load_expr;
    load_expr << "(load \"" << path << "\")";

    cl_object result = cl_safe_eval(
        c_string_to_object(load_expr.str().c_str()),
        Cnil, Cnil
    );

    if (result == Cnil) {
        last_error_ = "Failed to load file: " + path;
        throw EclException(last_error_);
    }
}

std::string Runtime::getLastError()
{
    return last_error_;
}

bool Runtime::hasError()
{
    return !last_error_.empty();
}

void Runtime::clearError()
{
    last_error_.clear();
}

const std::string& Runtime::getLispPath()
{
    return lisp_path_;
}

} // namespace Lisp
} // namespace StreamVorti

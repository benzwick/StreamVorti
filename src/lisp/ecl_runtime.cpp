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
        // Note: geometry.lisp must load before sdl-macros.lisp because
        // sdl-macros.lisp references streamvorti.geometry:shape
        std::vector<std::string> core_files = {
            "packages.lisp",
            "geometry.lisp",
            "mesh.lisp",
            "boundaries.lisp",
            "sdl-macros.lisp",
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

    // Wrap expression in handler-case to capture errors
    // Returns (values result nil) on success, (values nil error-message) on failure
    std::ostringstream wrapped_expr;
    wrapped_expr << "(handler-case "
                 << "  (values " << expr << " nil) "
                 << "  (error (c) "
                 << "    (values nil (format nil \"~A\" c))))";

    cl_object form = c_string_to_object(wrapped_expr.str().c_str());
    cl_object result = cl_safe_eval(form, Cnil, Cnil);

    if (result == OBJNULL) {
        last_error_ = "Evaluation returned null object";
        throw EclException(last_error_);
    }

    // Get the second value (error message if any)
    cl_object error_val = ecl_nth_value(ecl_process_env(), 1);

    if (error_val != Cnil && ECL_STRINGP(error_val)) {
        // There was an error - extract the message
        cl_object str = si_coerce_to_base_string(error_val);
        last_error_ = std::string(ecl_base_string_pointer_safe(str),
                                   ecl_length(str));
        throw EclException("Evaluation error: " + last_error_);
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

    // Build load expression wrapped in handler-case to capture errors
    // Returns T on success, error message string on failure
    std::ostringstream load_expr;
    load_expr << "(handler-case "
              << "  (progn (load \"" << path << "\") t) "
              << "  (error (c) "
              << "    (format nil \"~A\" c)))";

    cl_object result = cl_safe_eval(
        c_string_to_object(load_expr.str().c_str()),
        Cnil, Cnil
    );

    // Check if result is a string (error message) or T (success)
    if (result == OBJNULL) {
        last_error_ = "ECL evaluation returned null for: " + path;
        throw EclException(last_error_);
    }

    // If result is a string, it's an error message
    if (ECL_STRINGP(result)) {
        // Convert ECL string to C++ string
        cl_object str = si_coerce_to_base_string(result);
        last_error_ = std::string(ecl_base_string_pointer_safe(str),
                                   ecl_length(str));
        throw EclException("Error loading " + path + ": " + last_error_);
    }

    // If result is NIL (not T), something went wrong
    if (result == Cnil) {
        last_error_ = "Load returned NIL for: " + path;
        throw EclException(last_error_);
    }

    // Success: result should be T
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

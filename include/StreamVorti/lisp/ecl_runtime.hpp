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
 * @file ecl_runtime.hpp
 * @brief ECL (Embeddable Common Lisp) runtime management
 *
 * Provides initialization, shutdown, and evaluation functions for the
 * embedded Common Lisp runtime used by the SDL (Simulation Definition Language).
 */

#ifndef STREAMVORTI_LISP_ECL_RUNTIME_HPP_
#define STREAMVORTI_LISP_ECL_RUNTIME_HPP_

#include <string>
#include <vector>
#include <stdexcept>

// Forward declare ECL types - actual ECL header included only in .cpp files
// to avoid ECL macro conflicts with MFEM headers (e.g., ECL defines Ct macro
// which conflicts with mfem::Hybridization::Ct member variable)
// Use void* as opaque pointer - type safety enforced in implementation
typedef void* EclObject;

namespace StreamVorti {
namespace Lisp {

/**
 * @class EclException
 * @brief Exception class for ECL-related errors
 */
class EclException : public std::runtime_error {
public:
    explicit EclException(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @class Runtime
 * @brief Manages the ECL runtime lifecycle
 *
 * This class handles initialization and shutdown of the ECL runtime,
 * as well as providing methods for evaluating Lisp expressions.
 *
 * Usage:
 * @code
 * StreamVorti::Lisp::Runtime::init();
 * auto result = StreamVorti::Lisp::Runtime::eval("(+ 1 2)");
 * StreamVorti::Lisp::Runtime::shutdown();
 * @endcode
 */
class Runtime {
public:
    /**
     * @brief Initialize the ECL runtime
     *
     * Must be called before any Lisp operations. Safe to call multiple times.
     *
     * @param lisp_path Optional path to SDL Lisp source files
     * @throws EclException if initialization fails
     */
    static void init(const std::string& lisp_path = "");

    /**
     * @brief Shutdown the ECL runtime
     *
     * Should be called at program exit. Safe to call multiple times.
     */
    static void shutdown();

    /**
     * @brief Check if runtime is initialized
     * @return true if ECL has been initialized
     */
    static bool isInitialized();

    /**
     * @brief Evaluate a Lisp expression string
     *
     * @param expr The Lisp expression to evaluate
     * @return The result as an opaque ECL object pointer
     * @throws EclException on evaluation error
     */
    static EclObject eval(const std::string& expr);

    /**
     * @brief Safely evaluate a Lisp expression with error handling
     *
     * @param expr The Lisp expression to evaluate
     * @return The result as an opaque ECL object pointer, or nullptr on error
     */
    static EclObject safeEval(const std::string& expr);

    /**
     * @brief Load a Lisp file
     *
     * @param path Path to the .lisp file
     * @throws EclException if file cannot be loaded
     */
    static void loadFile(const std::string& path);

    /**
     * @brief Get the last error message
     * @return Error message string, empty if no error
     */
    static std::string getLastError();

    /**
     * @brief Check if there was an error in the last operation
     * @return true if there was an error
     */
    static bool hasError();

    /**
     * @brief Clear the error state
     */
    static void clearError();

    /**
     * @brief Get path to SDL Lisp source directory
     * @return Path string
     */
    static const std::string& getLispPath();

private:
    static bool initialized_;
    static std::string last_error_;
    static std::string lisp_path_;

    // Prevent instantiation
    Runtime() = delete;
    Runtime(const Runtime&) = delete;
    Runtime& operator=(const Runtime&) = delete;
};

} // namespace Lisp
} // namespace StreamVorti

#endif // STREAMVORTI_LISP_ECL_RUNTIME_HPP_

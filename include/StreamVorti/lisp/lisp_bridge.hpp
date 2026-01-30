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
 * @file lisp_bridge.hpp
 * @brief Type conversion bridge between C++ and Common Lisp
 *
 * Provides bidirectional type conversion between C++ types
 * (vectors, strings, numbers) and Lisp objects.
 */

#ifndef STREAMVORTI_LISP_LISP_BRIDGE_HPP_
#define STREAMVORTI_LISP_LISP_BRIDGE_HPP_

#include <string>
#include <vector>
#include <map>
#include <functional>

// Include ECL header for proper type definitions
#include <ecl/ecl.h>

// Forward declare MFEM types
namespace mfem {
    class Vector;
    class DenseMatrix;
    class SparseMatrix;
}

namespace StreamVorti {
namespace Lisp {

/**
 * @class Bridge
 * @brief Bidirectional type conversion between C++ and Lisp
 *
 * All methods are static. Requires ECL runtime to be initialized.
 */
class Bridge {
public:
    // ==================== C++ to Lisp ====================

    /**
     * @brief Convert double to Lisp double-float
     */
    static cl_object toDouble(double value);

    /**
     * @brief Convert int to Lisp integer
     */
    static cl_object toInt(int value);

    /**
     * @brief Convert int64 to Lisp integer
     */
    static cl_object toInt64(int64_t value);

    /**
     * @brief Convert bool to Lisp T/NIL
     */
    static cl_object toBool(bool value);

    /**
     * @brief Convert string to Lisp string
     */
    static cl_object toString(const std::string& value);

    /**
     * @brief Convert C++ vector to Lisp list
     */
    static cl_object toList(const std::vector<double>& vec);

    /**
     * @brief Convert C++ vector of ints to Lisp list
     */
    static cl_object toIntList(const std::vector<int>& vec);

    /**
     * @brief Convert MFEM Vector to Lisp list
     */
    static cl_object toList(const mfem::Vector& vec);

    /**
     * @brief Convert 2D array (row-major) to Lisp nested list
     */
    static cl_object toNestedList(const double* data, int rows, int cols);

    /**
     * @brief Create a Lisp keyword symbol
     */
    static cl_object toKeyword(const std::string& name);

    /**
     * @brief Create a Lisp symbol in current package
     */
    static cl_object toSymbol(const std::string& name);

    /**
     * @brief Create a Lisp cons cell
     */
    static cl_object cons(cl_object car, cl_object cdr);

    /**
     * @brief Create a Lisp list from multiple objects
     */
    static cl_object makeList(std::initializer_list<cl_object> items);

    // ==================== Lisp to C++ ====================

    /**
     * @brief Convert Lisp number to double
     */
    static double toCppDouble(cl_object obj);

    /**
     * @brief Convert Lisp integer to int
     */
    static int toCppInt(cl_object obj);

    /**
     * @brief Convert Lisp integer to int64
     */
    static int64_t toCppInt64(cl_object obj);

    /**
     * @brief Convert Lisp T/NIL to bool
     */
    static bool toCppBool(cl_object obj);

    /**
     * @brief Convert Lisp string to C++ string
     */
    static std::string toCppString(cl_object obj);

    /**
     * @brief Convert Lisp list to C++ vector of doubles
     */
    static std::vector<double> toVector(cl_object list);

    /**
     * @brief Convert Lisp list to C++ vector of ints
     */
    static std::vector<int> toIntVector(cl_object list);

    /**
     * @brief Get length of Lisp list
     */
    static size_t listLength(cl_object list);

    /**
     * @brief Get nth element of Lisp list (0-indexed)
     */
    static cl_object nth(cl_object list, size_t index);

    /**
     * @brief Get symbol name as string
     */
    static std::string symbolName(cl_object sym);

    // ==================== Type Checking ====================

    /**
     * @brief Check if object is NIL
     */
    static bool isNil(cl_object obj);

    /**
     * @brief Check if object is a number
     */
    static bool isNumber(cl_object obj);

    /**
     * @brief Check if object is a string
     */
    static bool isString(cl_object obj);

    /**
     * @brief Check if object is a list
     */
    static bool isList(cl_object obj);

    /**
     * @brief Check if object is a symbol
     */
    static bool isSymbol(cl_object obj);

    /**
     * @brief Check if object is a function
     */
    static bool isFunction(cl_object obj);

    // ==================== Function Calls ====================

    /**
     * @brief Call a Lisp function by name with arguments
     */
    static cl_object funcall(const std::string& func_name,
                             std::initializer_list<cl_object> args);

    /**
     * @brief Call a Lisp function object with arguments
     */
    static cl_object funcall(cl_object func,
                             std::initializer_list<cl_object> args);

    /**
     * @brief Apply a Lisp function to a list of arguments
     */
    static cl_object apply(cl_object func, cl_object args_list);

    /**
     * @brief Look up a symbol in a package
     */
    static cl_object findSymbol(const std::string& name,
                                const std::string& package = "SDL");

    /**
     * @brief Get the function value of a symbol
     */
    static cl_object symbolFunction(cl_object sym);

private:
    Bridge() = delete;
    Bridge(const Bridge&) = delete;
    Bridge& operator=(const Bridge&) = delete;
};

/**
 * @class LispFunction
 * @brief Wrapper for calling Lisp functions from C++
 *
 * Provides a convenient way to store and invoke Lisp functions.
 */
class LispFunction {
public:
    /**
     * @brief Create from a Lisp function object
     */
    explicit LispFunction(cl_object func);

    /**
     * @brief Create from a function name in a package
     */
    LispFunction(const std::string& name,
                 const std::string& package = "SDL");

    /**
     * @brief Check if function is valid
     */
    bool isValid() const;

    /**
     * @brief Call with no arguments
     */
    cl_object operator()() const;

    /**
     * @brief Call with arguments
     */
    cl_object operator()(std::initializer_list<cl_object> args) const;

    /**
     * @brief Evaluate boundary condition at point (convenience method)
     */
    double evaluateAt(double x, double y, double z = 0.0) const;

    /**
     * @brief Get the underlying Lisp object
     */
    cl_object get() const { return func_; }

private:
    cl_object func_;
};

} // namespace Lisp
} // namespace StreamVorti

#endif // STREAMVORTI_LISP_LISP_BRIDGE_HPP_

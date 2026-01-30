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

// Forward declare ECL types - actual ECL header included only in .cpp files
// to avoid ECL macro conflicts with MFEM headers (e.g., ECL defines Ct macro
// which conflicts with mfem::Hybridization::Ct member variable)
// Use void* as opaque pointer - type safety enforced in implementation
typedef void* EclObject;

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
    static EclObject toDouble(double value);

    /**
     * @brief Convert int to Lisp integer
     */
    static EclObject toInt(int value);

    /**
     * @brief Convert int64 to Lisp integer
     */
    static EclObject toInt64(int64_t value);

    /**
     * @brief Convert bool to Lisp T/NIL
     */
    static EclObject toBool(bool value);

    /**
     * @brief Convert string to Lisp string
     */
    static EclObject toString(const std::string& value);

    /**
     * @brief Convert C++ vector to Lisp list
     */
    static EclObject toList(const std::vector<double>& vec);

    /**
     * @brief Convert C++ vector of ints to Lisp list
     */
    static EclObject toIntList(const std::vector<int>& vec);

    /**
     * @brief Convert MFEM Vector to Lisp list
     */
    static EclObject toList(const mfem::Vector& vec);

    /**
     * @brief Convert 2D array (row-major) to Lisp nested list
     */
    static EclObject toNestedList(const double* data, int rows, int cols);

    /**
     * @brief Create a Lisp keyword symbol
     */
    static EclObject toKeyword(const std::string& name);

    /**
     * @brief Create a Lisp symbol in current package
     */
    static EclObject toSymbol(const std::string& name);

    /**
     * @brief Create a Lisp cons cell
     */
    static EclObject cons(EclObject car, EclObject cdr);

    /**
     * @brief Create a Lisp list from multiple objects
     */
    static EclObject makeList(std::initializer_list<EclObject> items);

    // ==================== Lisp to C++ ====================

    /**
     * @brief Convert Lisp number to double
     */
    static double toCppDouble(EclObject obj);

    /**
     * @brief Convert Lisp integer to int
     */
    static int toCppInt(EclObject obj);

    /**
     * @brief Convert Lisp integer to int64
     */
    static int64_t toCppInt64(EclObject obj);

    /**
     * @brief Convert Lisp T/NIL to bool
     */
    static bool toCppBool(EclObject obj);

    /**
     * @brief Convert Lisp string to C++ string
     */
    static std::string toCppString(EclObject obj);

    /**
     * @brief Convert Lisp list to C++ vector of doubles
     */
    static std::vector<double> toVector(EclObject list);

    /**
     * @brief Convert Lisp list to C++ vector of ints
     */
    static std::vector<int> toIntVector(EclObject list);

    /**
     * @brief Get length of Lisp list
     */
    static size_t listLength(EclObject list);

    /**
     * @brief Get nth element of Lisp list (0-indexed)
     */
    static EclObject nth(EclObject list, size_t index);

    /**
     * @brief Get symbol name as string
     */
    static std::string symbolName(EclObject sym);

    // ==================== Type Checking ====================

    /**
     * @brief Check if object is NIL
     */
    static bool isNil(EclObject obj);

    /**
     * @brief Check if object is a number
     */
    static bool isNumber(EclObject obj);

    /**
     * @brief Check if object is a string
     */
    static bool isString(EclObject obj);

    /**
     * @brief Check if object is a list
     */
    static bool isList(EclObject obj);

    /**
     * @brief Check if object is a symbol
     */
    static bool isSymbol(EclObject obj);

    /**
     * @brief Check if object is a function
     */
    static bool isFunction(EclObject obj);

    // ==================== Function Calls ====================

    /**
     * @brief Call a Lisp function by name with arguments
     */
    static EclObject funcall(const std::string& func_name,
                             std::initializer_list<EclObject> args);

    /**
     * @brief Call a Lisp function object with arguments
     */
    static EclObject funcall(EclObject func,
                             std::initializer_list<EclObject> args);

    /**
     * @brief Apply a Lisp function to a list of arguments
     */
    static EclObject apply(EclObject func, EclObject args_list);

    /**
     * @brief Look up a symbol in a package
     */
    static EclObject findSymbol(const std::string& name,
                                const std::string& package = "SDL");

    /**
     * @brief Get the function value of a symbol
     */
    static EclObject symbolFunction(EclObject sym);

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
    explicit LispFunction(EclObject func);

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
    EclObject operator()() const;

    /**
     * @brief Call with arguments
     */
    EclObject operator()(std::initializer_list<EclObject> args) const;

    /**
     * @brief Evaluate boundary condition at point (convenience method)
     */
    double evaluateAt(double x, double y, double z = 0.0) const;

    /**
     * @brief Get the underlying Lisp object
     */
    EclObject get() const { return func_; }

private:
    EclObject func_;
};

} // namespace Lisp
} // namespace StreamVorti

#endif // STREAMVORTI_LISP_LISP_BRIDGE_HPP_

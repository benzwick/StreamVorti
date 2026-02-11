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
 * @file lisp_bridge.cpp
 * @brief Implementation of C++/Lisp type conversion bridge
 */

// IMPORTANT: Include MFEM FIRST before ECL to avoid macro conflicts
// ECL defines macros (like Ct) that conflict with MFEM member names
#include "mfem.hpp"

// Now include ECL
#include <ecl/ecl.h>

// Include our headers after MFEM and ECL
#include "StreamVorti/lisp/lisp_bridge.hpp"
#include "StreamVorti/lisp/ecl_runtime.hpp"

#include <cctype>  // for std::toupper

namespace StreamVorti {
namespace Lisp {

// Helper functions for type conversion between EclObject (void*) and cl_object
static inline cl_object to_cl(EclObject x) { return reinterpret_cast<cl_object>(x); }
static inline EclObject to_ecl(cl_object x) { return reinterpret_cast<EclObject>(x); }

// ==================== C++ to Lisp ====================

EclObject Bridge::toDouble(double value)
{
    return to_ecl(ecl_make_double_float(value));
}

EclObject Bridge::toInt(int value)
{
    return to_ecl(ecl_make_fixnum(value));
}

EclObject Bridge::toInt64(int64_t value)
{
    return to_ecl(ecl_make_int64_t(value));
}

EclObject Bridge::toBool(bool value)
{
    return to_ecl(value ? ECL_T : ECL_NIL);
}

EclObject Bridge::toString(const std::string& value)
{
    return to_ecl(ecl_make_simple_base_string(value.c_str(), value.length()));
}

EclObject Bridge::toList(const std::vector<double>& vec)
{
    if (vec.empty()) {
        return to_ecl(ECL_NIL);
    }

    // Build list in reverse order, then reverse (more efficient)
    cl_object result = ECL_NIL;
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(to_cl(toDouble(*it)), result);
    }
    return to_ecl(result);
}

EclObject Bridge::toIntList(const std::vector<int>& vec)
{
    if (vec.empty()) {
        return to_ecl(ECL_NIL);
    }

    cl_object result = ECL_NIL;
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(to_cl(toInt(*it)), result);
    }
    return to_ecl(result);
}

EclObject Bridge::toList(const mfem::Vector& vec)
{
    if (vec.Size() == 0) {
        return to_ecl(ECL_NIL);
    }

    cl_object result = ECL_NIL;
    for (int i = vec.Size() - 1; i >= 0; --i) {
        result = ecl_cons(to_cl(toDouble(vec(i))), result);
    }
    return to_ecl(result);
}

EclObject Bridge::toNestedList(const double* data, int rows, int cols)
{
    cl_object result = ECL_NIL;

    // Build rows from last to first
    for (int i = rows - 1; i >= 0; --i) {
        cl_object row = ECL_NIL;
        for (int j = cols - 1; j >= 0; --j) {
            row = ecl_cons(to_cl(toDouble(data[i * cols + j])), row);
        }
        result = ecl_cons(row, result);
    }

    return to_ecl(result);
}

EclObject Bridge::toKeyword(const std::string& name)
{
    cl_object keyword_pkg = ecl_find_package("KEYWORD");
    cl_object name_str = ecl_make_simple_base_string(name.c_str(), name.length());
    int intern_flag;
    return to_ecl(ecl_intern(name_str, keyword_pkg, &intern_flag));
}

EclObject Bridge::toSymbol(const std::string& name)
{
    return to_ecl(ecl_read_from_cstring(name.c_str()));
}

EclObject Bridge::cons(EclObject car, EclObject cdr)
{
    return to_ecl(ecl_cons(to_cl(car), to_cl(cdr)));
}

EclObject Bridge::makeList(std::initializer_list<EclObject> items)
{
    cl_object result = ECL_NIL;
    // Build in reverse
    std::vector<EclObject> vec(items);
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(to_cl(*it), result);
    }
    return to_ecl(result);
}

// ==================== Lisp to C++ ====================

double Bridge::toCppDouble(EclObject obj)
{
    return ecl_to_double(to_cl(obj));
}

int Bridge::toCppInt(EclObject obj)
{
    return ecl_to_int(to_cl(obj));
}

int64_t Bridge::toCppInt64(EclObject obj)
{
    return ecl_to_int64_t(to_cl(obj));
}

bool Bridge::toCppBool(EclObject obj)
{
    return to_cl(obj) != ECL_NIL;
}

std::string Bridge::toCppString(EclObject obj)
{
    if (isNil(obj)) {
        return "";
    }

    cl_object cl_obj = to_cl(obj);

    if (isString(obj)) {
        cl_index length = ecl_length(cl_obj);
        std::string result;
        result.reserve(length);
        for (cl_index i = 0; i < length; ++i) {
            result.push_back(ecl_char(cl_obj, i));
        }
        return result;
    }

    // If it's a symbol, get its name
    if (isSymbol(obj)) {
        return toCppString(to_ecl(ecl_symbol_name(cl_obj)));
    }

    return "";
}

std::vector<double> Bridge::toVector(EclObject list)
{
    std::vector<double> result;

    if (isNil(list)) {
        return result;
    }

    cl_object cl_list = to_cl(list);
    result.reserve(listLength(list));

    while (cl_list != ECL_NIL && !isNil(to_ecl(cl_consp(cl_list)))) {
        result.push_back(toCppDouble(to_ecl(cl_car(cl_list))));
        cl_list = cl_cdr(cl_list);
    }

    return result;
}

std::vector<int> Bridge::toIntVector(EclObject list)
{
    std::vector<int> result;

    if (isNil(list)) {
        return result;
    }

    cl_object cl_list = to_cl(list);
    result.reserve(listLength(list));

    while (cl_list != ECL_NIL && !isNil(to_ecl(cl_consp(cl_list)))) {
        result.push_back(toCppInt(to_ecl(cl_car(cl_list))));
        cl_list = cl_cdr(cl_list);
    }

    return result;
}

size_t Bridge::listLength(EclObject list)
{
    if (isNil(list)) {
        return 0;
    }
    return ecl_length(to_cl(list));
}

EclObject Bridge::nth(EclObject list, size_t index)
{
    return to_ecl(ecl_nth(index, to_cl(list)));
}

std::string Bridge::symbolName(EclObject sym)
{
    if (!isSymbol(sym)) {
        return "";
    }
    return toCppString(to_ecl(ecl_symbol_name(to_cl(sym))));
}

// ==================== Type Checking ====================

bool Bridge::isNil(EclObject obj)
{
    cl_object cl_obj = to_cl(obj);
    return cl_obj == ECL_NIL || cl_obj == NULL;
}

bool Bridge::isNumber(EclObject obj)
{
    return !isNil(to_ecl(cl_numberp(to_cl(obj))));
}

bool Bridge::isString(EclObject obj)
{
    return !isNil(to_ecl(cl_stringp(to_cl(obj))));
}

bool Bridge::isList(EclObject obj)
{
    return !isNil(to_ecl(cl_listp(to_cl(obj))));
}

bool Bridge::isSymbol(EclObject obj)
{
    return !isNil(to_ecl(cl_symbolp(to_cl(obj))));
}

bool Bridge::isFunction(EclObject obj)
{
    return !isNil(to_ecl(cl_functionp(to_cl(obj))));
}

// ==================== Function Calls ====================

EclObject Bridge::funcall(const std::string& func_name,
                          std::initializer_list<EclObject> args)
{
    EclObject sym = findSymbol(func_name);
    if (isNil(sym)) {
        throw EclException("Function not found: " + func_name);
    }
    return funcall(symbolFunction(sym), args);
}

EclObject Bridge::funcall(EclObject func,
                          std::initializer_list<EclObject> args)
{
    size_t nargs = args.size();
    cl_object cl_func = to_cl(func);

    switch (nargs) {
        case 0:
            return to_ecl(::cl_funcall(1, cl_func));
        case 1: {
            auto it = args.begin();
            return to_ecl(::cl_funcall(2, cl_func, to_cl(*it)));
        }
        case 2: {
            auto it = args.begin();
            cl_object a1 = to_cl(*it++);
            cl_object a2 = to_cl(*it);
            return to_ecl(::cl_funcall(3, cl_func, a1, a2));
        }
        case 3: {
            auto it = args.begin();
            cl_object a1 = to_cl(*it++);
            cl_object a2 = to_cl(*it++);
            cl_object a3 = to_cl(*it);
            return to_ecl(::cl_funcall(4, cl_func, a1, a2, a3));
        }
        case 4: {
            auto it = args.begin();
            cl_object a1 = to_cl(*it++);
            cl_object a2 = to_cl(*it++);
            cl_object a3 = to_cl(*it++);
            cl_object a4 = to_cl(*it);
            return to_ecl(::cl_funcall(5, cl_func, a1, a2, a3, a4));
        }
        default: {
            // For more arguments, use apply with a list
            EclObject args_list = makeList(args);
            return apply(func, args_list);
        }
    }
}

EclObject Bridge::apply(EclObject func, EclObject args_list)
{
    return to_ecl(cl_apply(2, to_cl(func), to_cl(args_list)));
}

EclObject Bridge::findSymbol(const std::string& name,
                             const std::string& package)
{
    std::string actual_name = name;
    std::string actual_package = package;

    // Parse "package:name" format if present
    size_t colon_pos = name.find(':');
    if (colon_pos != std::string::npos) {
        actual_package = name.substr(0, colon_pos);
        actual_name = name.substr(colon_pos + 1);
    }

    // Convert to uppercase (Common Lisp default readtable behavior)
    for (char& c : actual_name) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
    for (char& c : actual_package) {
        c = std::toupper(static_cast<unsigned char>(c));
    }

    cl_object pkg = ecl_find_package(actual_package.c_str());
    if (pkg == ECL_NIL) {
        return to_ecl(ECL_NIL);
    }

    int intern_flag;
    cl_object name_str = ecl_make_simple_base_string(actual_name.c_str(),
                                                      actual_name.length());
    return to_ecl(ecl_intern(name_str, pkg, &intern_flag));
}

EclObject Bridge::symbolFunction(EclObject sym)
{
    return to_ecl(ecl_fdefinition(to_cl(sym)));
}

// ==================== LispFunction ====================

LispFunction::LispFunction(EclObject func)
    : func_(func)
{
    registerRoot();
}

LispFunction::LispFunction(const std::string& name,
                           const std::string& package)
{
    EclObject sym = Bridge::findSymbol(name, package);
    if (!Bridge::isNil(sym)) {
        // Check if symbol has a function binding before trying to get it
        // ecl_fdefinition will signal an error if function is undefined
        try {
            func_ = Bridge::symbolFunction(sym);
            // Verify it's actually a function
            if (!Bridge::isFunction(func_)) {
                func_ = nullptr;
            }
        } catch (...) {
            func_ = nullptr;
        }
    } else {
        func_ = nullptr;
    }
    registerRoot();
}

LispFunction::~LispFunction()
{
    deregisterRoot();
}

LispFunction::LispFunction(LispFunction&& other) noexcept
    : func_(other.func_)
{
    registerRoot();
    other.func_ = nullptr;
}

LispFunction& LispFunction::operator=(LispFunction&& other) noexcept
{
    if (this != &other) {
        deregisterRoot();
        func_ = other.func_;
        registerRoot();
        other.func_ = nullptr;
    }
    return *this;
}

void LispFunction::registerRoot()
{
    if (func_ != nullptr) {
        ecl_register_root(reinterpret_cast<cl_object*>(&func_));
    }
}

void LispFunction::deregisterRoot()
{
    // Set to nil so the GC root (which persists) doesn't keep a stale object alive
    func_ = to_ecl(ECL_NIL);
}

bool LispFunction::isValid() const
{
    return func_ != nullptr && !Bridge::isNil(func_) &&
           Bridge::isFunction(func_);
}

EclObject LispFunction::operator()() const
{
    if (!isValid()) {
        throw EclException("Invalid Lisp function");
    }
    return to_ecl(::cl_funcall(1, to_cl(func_)));
}

EclObject LispFunction::operator()(std::initializer_list<EclObject> args) const
{
    if (!isValid()) {
        throw EclException("Invalid Lisp function");
    }
    return Bridge::funcall(func_, args);
}

double LispFunction::evaluateAt(double x, double y, double z) const
{
    EclObject result = Bridge::funcall(func_, {
        Bridge::toDouble(x),
        Bridge::toDouble(y),
        Bridge::toDouble(z)
    });
    return Bridge::toCppDouble(result);
}

} // namespace Lisp
} // namespace StreamVorti

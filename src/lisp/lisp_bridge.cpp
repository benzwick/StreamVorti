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
 * @file lisp_bridge.cpp
 * @brief Implementation of C++/Lisp type conversion bridge
 */

#include "StreamVorti/lisp/lisp_bridge.hpp"
#include "StreamVorti/lisp/ecl_runtime.hpp"

#include <ecl/ecl.h>
#include "mfem.hpp"

namespace StreamVorti {
namespace Lisp {

// ==================== C++ to Lisp ====================

cl_object Bridge::toDouble(double value)
{
    return ecl_make_double_float(value);
}

cl_object Bridge::toInt(int value)
{
    return ecl_make_fixnum(value);
}

cl_object Bridge::toInt64(int64_t value)
{
    return ecl_make_int64_t(value);
}

cl_object Bridge::toBool(bool value)
{
    return value ? ECL_T : ECL_NIL;
}

cl_object Bridge::toString(const std::string& value)
{
    return ecl_make_simple_base_string(value.c_str(), value.length());
}

cl_object Bridge::toList(const std::vector<double>& vec)
{
    if (vec.empty()) {
        return ECL_NIL;
    }

    // Build list in reverse order, then reverse (more efficient)
    cl_object result = ECL_NIL;
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(toDouble(*it), result);
    }
    return result;
}

cl_object Bridge::toIntList(const std::vector<int>& vec)
{
    if (vec.empty()) {
        return ECL_NIL;
    }

    cl_object result = ECL_NIL;
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(toInt(*it), result);
    }
    return result;
}

cl_object Bridge::toList(const mfem::Vector& vec)
{
    if (vec.Size() == 0) {
        return ECL_NIL;
    }

    cl_object result = ECL_NIL;
    for (int i = vec.Size() - 1; i >= 0; --i) {
        result = ecl_cons(toDouble(vec(i)), result);
    }
    return result;
}

cl_object Bridge::toNestedList(const double* data, int rows, int cols)
{
    cl_object result = ECL_NIL;

    // Build rows from last to first
    for (int i = rows - 1; i >= 0; --i) {
        cl_object row = ECL_NIL;
        for (int j = cols - 1; j >= 0; --j) {
            row = ecl_cons(toDouble(data[i * cols + j]), row);
        }
        result = ecl_cons(row, result);
    }

    return result;
}

cl_object Bridge::toKeyword(const std::string& name)
{
    cl_object keyword_pkg = ecl_find_package("KEYWORD");
    return ecl_intern(name.c_str(), name.length(), keyword_pkg, NULL);
}

cl_object Bridge::toSymbol(const std::string& name)
{
    return ecl_read_from_cstring(name.c_str());
}

cl_object Bridge::cons(cl_object car, cl_object cdr)
{
    return ecl_cons(car, cdr);
}

cl_object Bridge::makeList(std::initializer_list<cl_object> items)
{
    cl_object result = ECL_NIL;
    // Build in reverse
    std::vector<cl_object> vec(items);
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        result = ecl_cons(*it, result);
    }
    return result;
}

// ==================== Lisp to C++ ====================

double Bridge::toCppDouble(cl_object obj)
{
    return ecl_to_double(obj);
}

int Bridge::toCppInt(cl_object obj)
{
    return ecl_to_int(obj);
}

int64_t Bridge::toCppInt64(cl_object obj)
{
    return ecl_to_int64_t(obj);
}

bool Bridge::toCppBool(cl_object obj)
{
    return obj != ECL_NIL;
}

std::string Bridge::toCppString(cl_object obj)
{
    if (isNil(obj)) {
        return "";
    }

    if (isString(obj)) {
        cl_index length = ecl_length(obj);
        std::string result;
        result.reserve(length);
        for (cl_index i = 0; i < length; ++i) {
            result.push_back(ecl_char(obj, i));
        }
        return result;
    }

    // If it's a symbol, get its name
    if (isSymbol(obj)) {
        return toCppString(ecl_symbol_name(obj));
    }

    return "";
}

std::vector<double> Bridge::toVector(cl_object list)
{
    std::vector<double> result;

    if (isNil(list)) {
        return result;
    }

    result.reserve(listLength(list));

    while (!isNil(list) && !isNil(cl_consp(list))) {
        result.push_back(toCppDouble(cl_car(list)));
        list = cl_cdr(list);
    }

    return result;
}

std::vector<int> Bridge::toIntVector(cl_object list)
{
    std::vector<int> result;

    if (isNil(list)) {
        return result;
    }

    result.reserve(listLength(list));

    while (!isNil(list) && !isNil(cl_consp(list))) {
        result.push_back(toCppInt(cl_car(list)));
        list = cl_cdr(list);
    }

    return result;
}

size_t Bridge::listLength(cl_object list)
{
    if (isNil(list)) {
        return 0;
    }
    return ecl_length(list);
}

cl_object Bridge::nth(cl_object list, size_t index)
{
    return ecl_nth(index, list);
}

std::string Bridge::symbolName(cl_object sym)
{
    if (!isSymbol(sym)) {
        return "";
    }
    return toCppString(ecl_symbol_name(sym));
}

// ==================== Type Checking ====================

bool Bridge::isNil(cl_object obj)
{
    return obj == ECL_NIL || obj == NULL;
}

bool Bridge::isNumber(cl_object obj)
{
    return !isNil(cl_numberp(obj));
}

bool Bridge::isString(cl_object obj)
{
    return !isNil(cl_stringp(obj));
}

bool Bridge::isList(cl_object obj)
{
    return !isNil(cl_listp(obj));
}

bool Bridge::isSymbol(cl_object obj)
{
    return !isNil(cl_symbolp(obj));
}

bool Bridge::isFunction(cl_object obj)
{
    return !isNil(cl_functionp(obj));
}

// ==================== Function Calls ====================

cl_object Bridge::funcall(const std::string& func_name,
                          std::initializer_list<cl_object> args)
{
    cl_object sym = findSymbol(func_name);
    if (isNil(sym)) {
        throw EclException("Function not found: " + func_name);
    }
    return funcall(symbolFunction(sym), args);
}

cl_object Bridge::funcall(cl_object func,
                          std::initializer_list<cl_object> args)
{
    size_t nargs = args.size();

    switch (nargs) {
        case 0:
            return cl_funcall(1, func);
        case 1: {
            auto it = args.begin();
            return cl_funcall(2, func, *it);
        }
        case 2: {
            auto it = args.begin();
            cl_object a1 = *it++;
            cl_object a2 = *it;
            return cl_funcall(3, func, a1, a2);
        }
        case 3: {
            auto it = args.begin();
            cl_object a1 = *it++;
            cl_object a2 = *it++;
            cl_object a3 = *it;
            return cl_funcall(4, func, a1, a2, a3);
        }
        case 4: {
            auto it = args.begin();
            cl_object a1 = *it++;
            cl_object a2 = *it++;
            cl_object a3 = *it++;
            cl_object a4 = *it;
            return cl_funcall(5, func, a1, a2, a3, a4);
        }
        default: {
            // For more arguments, use apply with a list
            cl_object args_list = makeList(args);
            return apply(func, args_list);
        }
    }
}

cl_object Bridge::apply(cl_object func, cl_object args_list)
{
    return cl_apply(2, func, args_list);
}

cl_object Bridge::findSymbol(const std::string& name,
                             const std::string& package)
{
    cl_object pkg = ecl_find_package(package.c_str());
    if (isNil(pkg)) {
        return ECL_NIL;
    }

    int intern_flag;
    cl_object name_str = ecl_make_simple_base_string(name.c_str(), name.length());
    return ecl_intern(name_str, pkg, &intern_flag);
}

cl_object Bridge::symbolFunction(cl_object sym)
{
    return ecl_fdefinition(sym);
}

// ==================== LispFunction ====================

LispFunction::LispFunction(cl_object func)
    : func_(func)
{
}

LispFunction::LispFunction(const std::string& name,
                           const std::string& package)
{
    cl_object sym = Bridge::findSymbol(name, package);
    if (!Bridge::isNil(sym)) {
        func_ = Bridge::symbolFunction(sym);
    } else {
        func_ = nullptr;
    }
}

bool LispFunction::isValid() const
{
    return func_ != nullptr && !Bridge::isNil(func_) &&
           Bridge::isFunction(func_);
}

cl_object LispFunction::operator()() const
{
    if (!isValid()) {
        throw EclException("Invalid Lisp function");
    }
    return cl_funcall(1, func_);
}

cl_object LispFunction::operator()(std::initializer_list<cl_object> args) const
{
    if (!isValid()) {
        throw EclException("Invalid Lisp function");
    }
    return Bridge::funcall(func_, args);
}

double LispFunction::evaluateAt(double x, double y, double z) const
{
    cl_object result = Bridge::funcall(func_, {
        Bridge::toDouble(x),
        Bridge::toDouble(y),
        Bridge::toDouble(z)
    });
    return Bridge::toCppDouble(result);
}

} // namespace Lisp
} // namespace StreamVorti

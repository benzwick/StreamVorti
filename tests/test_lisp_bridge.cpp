/**
 * @file test_lisp_bridge.cpp
 * @brief Unit tests for C++ <-> Lisp type conversion bridge
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

#ifdef STREAMVORTI_WITH_ECL

#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_bridge.hpp>

#include <ecl/ecl.h>
#include <vector>
#include <string>

// Helper macros to convert between EclObject and cl_object for ECL-specific operations
#define TO_CL(x) reinterpret_cast<cl_object>(x)
#define TO_ECL(x) reinterpret_cast<EclObject>(x)

class LispBridgeTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!StreamVorti::Lisp::Runtime::isInitialized()) {
            StreamVorti::Lisp::Runtime::init();
        }
    }
};

// ==================== C++ to Lisp Conversion ====================

TEST_F(LispBridgeTest, ToDoubleConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toDouble(3.14159);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNumber(obj));
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(obj), 3.14159, 1e-10);
}

TEST_F(LispBridgeTest, ToIntConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toInt(42);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNumber(obj));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(obj), 42);
}

TEST_F(LispBridgeTest, ToInt64Conversion) {
    int64_t large_num = 9223372036854775807LL;
    EclObject obj = StreamVorti::Lisp::Bridge::toInt64(large_num);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNumber(obj));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt64(obj), large_num);
}

TEST_F(LispBridgeTest, ToBoolTrueConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toBool(true);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::toCppBool(obj));
}

TEST_F(LispBridgeTest, ToBoolFalseConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toBool(false);
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::toCppBool(obj));
}

TEST_F(LispBridgeTest, ToStringConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toString("hello world");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isString(obj));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppString(obj), "hello world");
}

TEST_F(LispBridgeTest, ToStringEmptyConversion) {
    EclObject obj = StreamVorti::Lisp::Bridge::toString("");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isString(obj));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppString(obj), "");
}

TEST_F(LispBridgeTest, ToListDoubleConversion) {
    std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0};
    EclObject obj = StreamVorti::Lisp::Bridge::toList(vec);

    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isList(obj));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::listLength(obj), 5);

    auto result = StreamVorti::Lisp::Bridge::toVector(obj);
    EXPECT_EQ(result.size(), 5);
    for (size_t i = 0; i < vec.size(); ++i) {
        EXPECT_DOUBLE_EQ(result[i], vec[i]);
    }
}

TEST_F(LispBridgeTest, ToListEmptyConversion) {
    std::vector<double> vec;
    EclObject obj = StreamVorti::Lisp::Bridge::toList(vec);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNil(obj));
}

TEST_F(LispBridgeTest, ToIntListConversion) {
    std::vector<int> vec = {10, 20, 30};
    EclObject obj = StreamVorti::Lisp::Bridge::toIntList(vec);

    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isList(obj));
    auto result = StreamVorti::Lisp::Bridge::toIntVector(obj);
    EXPECT_EQ(result.size(), 3);
    EXPECT_EQ(result[0], 10);
    EXPECT_EQ(result[1], 20);
    EXPECT_EQ(result[2], 30);
}

// ==================== Lisp to C++ Conversion ====================

TEST_F(LispBridgeTest, FromLispInt) {
    EclObject obj = StreamVorti::Lisp::Runtime::eval("123");
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(obj), 123);
}

TEST_F(LispBridgeTest, FromLispDouble) {
    // Use d0 suffix to get double-float precision in Common Lisp
    EclObject obj = StreamVorti::Lisp::Runtime::eval("2.718281828d0");
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(obj), 2.718281828, 1e-8);
}

TEST_F(LispBridgeTest, FromLispString) {
    EclObject obj = StreamVorti::Lisp::Runtime::eval("\"test string\"");
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppString(obj), "test string");
}

TEST_F(LispBridgeTest, FromLispList) {
    EclObject obj = StreamVorti::Lisp::Runtime::eval("'(1.0 2.0 3.0)");
    auto vec = StreamVorti::Lisp::Bridge::toVector(obj);

    EXPECT_EQ(vec.size(), 3);
    EXPECT_DOUBLE_EQ(vec[0], 1.0);
    EXPECT_DOUBLE_EQ(vec[1], 2.0);
    EXPECT_DOUBLE_EQ(vec[2], 3.0);
}

// ==================== Type Checking ====================

TEST_F(LispBridgeTest, IsNil) {
    EclObject nil = StreamVorti::Lisp::Runtime::eval("nil");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNil(nil));

    EclObject not_nil = StreamVorti::Lisp::Runtime::eval("t");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::isNil(not_nil));
}

TEST_F(LispBridgeTest, IsNumber) {
    EclObject num = StreamVorti::Lisp::Runtime::eval("42");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isNumber(num));

    EclObject str = StreamVorti::Lisp::Runtime::eval("\"hello\"");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::isNumber(str));
}

TEST_F(LispBridgeTest, IsString) {
    EclObject str = StreamVorti::Lisp::Runtime::eval("\"hello\"");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isString(str));

    EclObject num = StreamVorti::Lisp::Runtime::eval("42");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::isString(num));
}

TEST_F(LispBridgeTest, IsList) {
    EclObject list = StreamVorti::Lisp::Runtime::eval("'(1 2 3)");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isList(list));

    EclObject num = StreamVorti::Lisp::Runtime::eval("42");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::isList(num));
}

TEST_F(LispBridgeTest, IsSymbol) {
    EclObject sym = StreamVorti::Lisp::Runtime::eval("'my-symbol");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isSymbol(sym));

    EclObject num = StreamVorti::Lisp::Runtime::eval("42");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::isSymbol(num));
}

// ==================== List Operations ====================

TEST_F(LispBridgeTest, ListLength) {
    EclObject list = StreamVorti::Lisp::Runtime::eval("'(a b c d e)");
    EXPECT_EQ(StreamVorti::Lisp::Bridge::listLength(list), 5);
}

TEST_F(LispBridgeTest, ListNth) {
    EclObject list = StreamVorti::Lisp::Runtime::eval("'(10 20 30 40)");

    EclObject first = StreamVorti::Lisp::Bridge::nth(list, 0);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(first), 10);

    EclObject third = StreamVorti::Lisp::Bridge::nth(list, 2);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(third), 30);
}

// ==================== Cons and MakeList ====================

TEST_F(LispBridgeTest, Cons) {
    EclObject a = StreamVorti::Lisp::Bridge::toInt(1);
    EclObject b = StreamVorti::Lisp::Bridge::toInt(2);
    EclObject ecl_nil = TO_ECL(ECL_NIL);  // unused, for reference

    // Create (1 . 2)
    EclObject pair = StreamVorti::Lisp::Bridge::cons(a, b);
    // Use TO_CL/TO_ECL macros to convert between EclObject and cl_object for ECL macros
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(TO_ECL(ECL_CONS_CAR(TO_CL(pair)))), 1);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(TO_ECL(ECL_CONS_CDR(TO_CL(pair)))), 2);
}

TEST_F(LispBridgeTest, MakeList) {
    EclObject list = StreamVorti::Lisp::Bridge::makeList({
        StreamVorti::Lisp::Bridge::toInt(1),
        StreamVorti::Lisp::Bridge::toInt(2),
        StreamVorti::Lisp::Bridge::toInt(3)
    });

    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isList(list));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::listLength(list), 3);

    auto vec = StreamVorti::Lisp::Bridge::toIntVector(list);
    EXPECT_EQ(vec[0], 1);
    EXPECT_EQ(vec[1], 2);
    EXPECT_EQ(vec[2], 3);
}

// ==================== LispFunction Tests ====================

TEST_F(LispBridgeTest, LispFunctionCall) {
    // Define a simple function
    StreamVorti::Lisp::Runtime::eval("(defun test-double (x) (* 2 x))");

    StreamVorti::Lisp::LispFunction func("TEST-DOUBLE", "CL-USER");
    EXPECT_TRUE(func.isValid());

    EclObject result = func({StreamVorti::Lisp::Bridge::toInt(21)});
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 42);
}

TEST_F(LispBridgeTest, LispFunctionEvaluateAt) {
    // Define a 2D function f(x,y) = x + y
    StreamVorti::Lisp::Runtime::eval(
        "(defun test-sum-2d (x y z) (+ x y))");

    StreamVorti::Lisp::LispFunction func("TEST-SUM-2D", "CL-USER");
    EXPECT_TRUE(func.isValid());

    double result = func.evaluateAt(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(result, 7.0);
}

#endif // STREAMVORTI_WITH_ECL

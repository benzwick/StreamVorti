/**
 * @file test_ecl_runtime.cpp
 * @brief Unit tests for ECL runtime integration
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

#ifdef STREAMVORTI_WITH_ECL

#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_bridge.hpp>

#include <ecl/ecl.h>

class EclRuntimeTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!StreamVorti::Lisp::Runtime::isInitialized()) {
            StreamVorti::Lisp::Runtime::init();
        }
    }

    void TearDown() override {
        // Don't shutdown ECL between tests - it can only be initialized once
    }
};

// Basic initialization tests
TEST_F(EclRuntimeTest, InitializationSucceeds) {
    EXPECT_TRUE(StreamVorti::Lisp::Runtime::isInitialized());
}

TEST_F(EclRuntimeTest, MultipleInitCalls) {
    // Should be safe to call init multiple times
    EXPECT_NO_THROW(StreamVorti::Lisp::Runtime::init());
    EXPECT_TRUE(StreamVorti::Lisp::Runtime::isInitialized());
}

// Basic evaluation tests
TEST_F(EclRuntimeTest, EvalSimpleArithmetic) {
    auto result = StreamVorti::Lisp::Runtime::eval("(+ 1 2)");
    EXPECT_NE(result, nullptr);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 3);
}

TEST_F(EclRuntimeTest, EvalMultiplication) {
    auto result = StreamVorti::Lisp::Runtime::eval("(* 6 7)");
    EXPECT_NE(result, nullptr);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 42);
}

TEST_F(EclRuntimeTest, EvalFloatingPoint) {
    auto result = StreamVorti::Lisp::Runtime::eval("(+ 1.5 2.5)");
    EXPECT_NE(result, nullptr);
    EXPECT_DOUBLE_EQ(StreamVorti::Lisp::Bridge::toCppDouble(result), 4.0);
}

TEST_F(EclRuntimeTest, EvalNestedExpressions) {
    auto result = StreamVorti::Lisp::Runtime::eval("(* (+ 2 3) (- 10 4))");
    EXPECT_NE(result, nullptr);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 30);
}

// Safe evaluation tests
TEST_F(EclRuntimeTest, SafeEvalSuccess) {
    auto result = StreamVorti::Lisp::Runtime::safeEval("(+ 1 1)");
    EXPECT_NE(result, nullptr);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 2);
}

// String evaluation
TEST_F(EclRuntimeTest, EvalString) {
    auto result = StreamVorti::Lisp::Runtime::eval("\"hello\"");
    EXPECT_NE(result, nullptr);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isString(result));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppString(result), "hello");
}

// List evaluation
TEST_F(EclRuntimeTest, EvalList) {
    auto result = StreamVorti::Lisp::Runtime::eval("'(1 2 3)");
    EXPECT_NE(result, nullptr);
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::isList(result));
    EXPECT_EQ(StreamVorti::Lisp::Bridge::listLength(result), 3);
}

// Boolean evaluation
TEST_F(EclRuntimeTest, EvalBooleanTrue) {
    auto result = StreamVorti::Lisp::Runtime::eval("t");
    EXPECT_TRUE(StreamVorti::Lisp::Bridge::toCppBool(result));
}

TEST_F(EclRuntimeTest, EvalBooleanFalse) {
    auto result = StreamVorti::Lisp::Runtime::eval("nil");
    EXPECT_FALSE(StreamVorti::Lisp::Bridge::toCppBool(result));
}

// Lambda evaluation
TEST_F(EclRuntimeTest, EvalLambda) {
    // Define a function
    StreamVorti::Lisp::Runtime::eval("(defun test-add (a b) (+ a b))");

    // Call the function
    auto result = StreamVorti::Lisp::Runtime::eval("(test-add 3 4)");
    EXPECT_NE(result, nullptr);
    EXPECT_EQ(StreamVorti::Lisp::Bridge::toCppInt(result), 7);
}

TEST_F(EclRuntimeTest, EvalComplexFunction) {
    // Define a parabolic profile function
    StreamVorti::Lisp::Runtime::eval(
        "(defun parabolic (y) (* 4.0 y (- 1.0 y)))");

    // Test at y = 0.5 (should give maximum = 1.0)
    auto result = StreamVorti::Lisp::Runtime::eval("(parabolic 0.5)");
    EXPECT_NE(result, nullptr);
    EXPECT_DOUBLE_EQ(StreamVorti::Lisp::Bridge::toCppDouble(result), 1.0);

    // Test at y = 0 (should give 0)
    result = StreamVorti::Lisp::Runtime::eval("(parabolic 0.0)");
    EXPECT_DOUBLE_EQ(StreamVorti::Lisp::Bridge::toCppDouble(result), 0.0);

    // Test at y = 1 (should give 0)
    result = StreamVorti::Lisp::Runtime::eval("(parabolic 1.0)");
    EXPECT_DOUBLE_EQ(StreamVorti::Lisp::Bridge::toCppDouble(result), 0.0);
}

// Error handling
TEST_F(EclRuntimeTest, ClearError) {
    StreamVorti::Lisp::Runtime::clearError();
    EXPECT_FALSE(StreamVorti::Lisp::Runtime::hasError());
    EXPECT_TRUE(StreamVorti::Lisp::Runtime::getLastError().empty());
}

// Math functions
TEST_F(EclRuntimeTest, EvalMathSin) {
    auto result = StreamVorti::Lisp::Runtime::eval("(sin 0.0)");
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(result), 0.0, 1e-10);
}

TEST_F(EclRuntimeTest, EvalMathCos) {
    auto result = StreamVorti::Lisp::Runtime::eval("(cos 0.0)");
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(result), 1.0, 1e-10);
}

TEST_F(EclRuntimeTest, EvalMathExp) {
    auto result = StreamVorti::Lisp::Runtime::eval("(exp 0.0)");
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(result), 1.0, 1e-10);
}

TEST_F(EclRuntimeTest, EvalMathSqrt) {
    auto result = StreamVorti::Lisp::Runtime::eval("(sqrt 4.0)");
    EXPECT_NEAR(StreamVorti::Lisp::Bridge::toCppDouble(result), 2.0, 1e-10);
}

#endif // STREAMVORTI_WITH_ECL

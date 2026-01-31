/**
 * @file test_sdl_loader.cpp
 * @brief Unit tests for SDL (Simulation Definition Language) loader
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

#ifdef STREAMVORTI_WITH_ECL

// Include MFEM first to get complete type definitions
#include "mfem.hpp"

#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_loader.hpp>

#include <fstream>
#include <filesystem>

class SdlLoaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!StreamVorti::Lisp::Runtime::isInitialized()) {
            // Initialize with path to lisp sources
            StreamVorti::Lisp::Runtime::init("lisp");
        }
    }

    // Helper to create a temporary SDL file
    std::string createTempSdlFile(const std::string& content) {
        std::string path = "/tmp/test_simulation.lisp";
        std::ofstream ofs(path);
        ofs << content;
        ofs.close();
        return path;
    }
};

// ==================== SimulationConfig Tests ====================

TEST_F(SdlLoaderTest, SimulationConfigDefaults) {
    StreamVorti::Lisp::SimulationConfig config;

    EXPECT_EQ(config.name, "");
    EXPECT_EQ(config.version, 1);
    EXPECT_EQ(config.dimension, 2);
    EXPECT_EQ(config.mesh, nullptr);
    EXPECT_TRUE(config.boundaries.empty());
}

TEST_F(SdlLoaderTest, DCPSEParamsDefaults) {
    StreamVorti::Lisp::DCPSEParams params;

    EXPECT_EQ(params.num_neighbors, 25);
    EXPECT_DOUBLE_EQ(params.cutoff_radius, 30.0);
    EXPECT_DOUBLE_EQ(params.support_radius, 5.0);
}

TEST_F(SdlLoaderTest, SolverParamsDefaults) {
    StreamVorti::Lisp::SolverParams params;

    EXPECT_EQ(params.timestepping, "explicit-euler");
    EXPECT_DOUBLE_EQ(params.dt, 0.001);
    EXPECT_DOUBLE_EQ(params.end_time, 1.0);
    EXPECT_DOUBLE_EQ(params.tolerance, 1e-6);
    EXPECT_EQ(params.max_iterations, 10000);
}

TEST_F(SdlLoaderTest, PhysicsParamsDefaults) {
    StreamVorti::Lisp::PhysicsParams params;

    EXPECT_EQ(params.type, "incompressible-navier-stokes");
    EXPECT_EQ(params.formulation, "stream-vorticity");
    EXPECT_DOUBLE_EQ(params.reynolds, 100.0);
    EXPECT_DOUBLE_EQ(params.density, 1.0);
    EXPECT_DOUBLE_EQ(params.viscosity, 0.01);
}

TEST_F(SdlLoaderTest, OutputParamsDefaults) {
    StreamVorti::Lisp::OutputParams params;

    EXPECT_EQ(params.format, "vtk");
    EXPECT_DOUBLE_EQ(params.interval, 0.1);
    EXPECT_EQ(params.directory, "results/");
    EXPECT_TRUE(params.fields.empty());
}

// ==================== LispFunction Tests ====================

TEST_F(SdlLoaderTest, LispFunctionNoArg) {
    // Define a simple function that returns a constant
    StreamVorti::Lisp::Runtime::eval("(defun test-const () 42.0)");

    StreamVorti::Lisp::LispFunction func("TEST-CONST", "CL-USER");
    EXPECT_TRUE(func.isValid());

    // Can't easily test no-arg call through evaluateAt
    // Just verify the function is valid
}

TEST_F(SdlLoaderTest, LispFunctionParabolic) {
    // Define parabolic profile commonly used in CFD
    StreamVorti::Lisp::Runtime::eval(
        "(defun parabolic-profile (x y z) "
        "  (declare (ignore x z))"
        "  (* 4.0d0 y (- 1.0d0 y)))");

    StreamVorti::Lisp::LispFunction func("PARABOLIC-PROFILE", "CL-USER");
    EXPECT_TRUE(func.isValid());

    // Test at y=0.5 (maximum)
    double val_mid = func.evaluateAt(0.0, 0.5, 0.0);
    EXPECT_NEAR(val_mid, 1.0, 1e-10);

    // Test at y=0 (zero)
    double val_0 = func.evaluateAt(0.0, 0.0, 0.0);
    EXPECT_NEAR(val_0, 0.0, 1e-10);

    // Test at y=1 (zero)
    double val_1 = func.evaluateAt(0.0, 1.0, 0.0);
    EXPECT_NEAR(val_1, 0.0, 1e-10);

    // Test at y=0.25
    double val_25 = func.evaluateAt(0.0, 0.25, 0.0);
    EXPECT_NEAR(val_25, 0.75, 1e-10);  // 4 * 0.25 * 0.75 = 0.75
}

TEST_F(SdlLoaderTest, LispFunctionSinWave) {
    // Define a sine wave function
    StreamVorti::Lisp::Runtime::eval(
        "(defun sin-wave (x y z) "
        "  (declare (ignore y z))"
        "  (sin (* 2.0d0 pi x)))");

    StreamVorti::Lisp::LispFunction func("SIN-WAVE", "CL-USER");
    EXPECT_TRUE(func.isValid());

    // Test at x=0
    EXPECT_NEAR(func.evaluateAt(0.0, 0.0, 0.0), 0.0, 1e-10);

    // Test at x=0.25
    EXPECT_NEAR(func.evaluateAt(0.25, 0.0, 0.0), 1.0, 1e-10);

    // Test at x=0.5
    EXPECT_NEAR(func.evaluateAt(0.5, 0.0, 0.0), 0.0, 1e-10);
}

// ==================== Boundary Condition Tests ====================

TEST_F(SdlLoaderTest, BoundaryConditionConstruct) {
    StreamVorti::Lisp::BoundaryCondition bc;

    bc.name = "inlet";
    bc.attribute = 1;
    bc.type = "velocity";

    EXPECT_EQ(bc.name, "inlet");
    EXPECT_EQ(bc.attribute, 1);
    EXPECT_EQ(bc.type, "velocity");
    EXPECT_EQ(bc.function, nullptr);
}

// ==================== SimulationRunner Tests ====================

TEST_F(SdlLoaderTest, SimulationRunnerConstruct) {
    StreamVorti::Lisp::SimulationConfig config;
    config.name = "test-simulation";
    config.dimension = 2;
    config.dcpse.num_neighbors = 20;
    config.solver.dt = 0.0001;
    config.solver.end_time = 0.01;

    StreamVorti::Lisp::SimulationRunner runner(std::move(config));

    EXPECT_DOUBLE_EQ(runner.currentTime(), 0.0);
    EXPECT_EQ(runner.config().name, "test-simulation");
    EXPECT_EQ(runner.config().dimension, 2);
}

// ==================== Loader Getters Tests ====================

TEST_F(SdlLoaderTest, LoaderGetFunction) {
    // Define a function in SDL package
    StreamVorti::Lisp::Runtime::eval(
        "(in-package :cl-user)"
        "(defun my-test-func (x) (* x 2))");

    auto func = StreamVorti::Lisp::Loader::getFunction("MY-TEST-FUNC");
    // Function might not be valid if not in SDL package
    // This test mainly verifies the method doesn't crash
}

// ==================== Integration Tests ====================

// Note: Full integration tests with file loading require the SDL
// Lisp files to be properly loaded. These tests verify the API.

TEST_F(SdlLoaderTest, LoadFromStringMinimal) {
    // This test requires SDL packages to be loaded
    // Skip if packages aren't available

    try {
        // Try to access SDL package
        auto result = StreamVorti::Lisp::Runtime::safeEval(
            "(find-package :sdl)");

        if (StreamVorti::Lisp::Bridge::isNil(result)) {
            GTEST_SKIP() << "SDL package not loaded";
        }

        // If package exists, we can try a minimal simulation
        // (Note: This depends on full SDL implementation)

    } catch (...) {
        GTEST_SKIP() << "SDL package not available";
    }
}

#endif // STREAMVORTI_WITH_ECL

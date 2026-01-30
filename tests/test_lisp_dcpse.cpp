/**
 * @file test_lisp_dcpse.cpp
 * @brief Unit tests for DCPSE bindings to Lisp
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

#ifdef STREAMVORTI_WITH_ECL

// Include MFEM first before ECL to avoid macro conflicts
#include "mfem.hpp"

#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_mesh.hpp>
#include <StreamVorti/lisp/lisp_dcpse.hpp>

// Include DCPSE headers for complete type definitions
#include <StreamVorti/approximants/dcpse_2d.hpp>
#include <StreamVorti/approximants/dcpse_3d.hpp>

class LispDcpseTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!StreamVorti::Lisp::Runtime::isInitialized()) {
            StreamVorti::Lisp::Runtime::init();
        }
    }
};

// ==================== DCPSE Creation Tests ====================

TEST_F(LispDcpseTest, CreateDcpse2D) {
    // Create mesh
    void* mesh = sv_make_cartesian_mesh_2d(5, 5, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    // Create GridFunction
    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    // Create DCPSE
    void* dcpse = sv_make_dcpse_2d(gf, 25);
    ASSERT_NE(dcpse, nullptr);

    EXPECT_EQ(sv_dcpse_dimension(dcpse), 2);
    EXPECT_EQ(sv_dcpse_num_neighbors(dcpse), 25);
    EXPECT_EQ(sv_dcpse_num_nodes(dcpse), 36);  // (5+1)^2 = 36

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

TEST_F(LispDcpseTest, CreateDcpse3D) {
    // Create mesh
    void* mesh = sv_make_cartesian_mesh_3d(3, 3, 3, 5, 1.0, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    // Create GridFunction
    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    // Create DCPSE
    void* dcpse = sv_make_dcpse_3d(gf, 30);
    ASSERT_NE(dcpse, nullptr);

    EXPECT_EQ(sv_dcpse_dimension(dcpse), 3);
    EXPECT_EQ(sv_dcpse_num_neighbors(dcpse), 30);

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== DCPSE Update Tests ====================

TEST_F(LispDcpseTest, UpdateDcpse2D) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 20);
    ASSERT_NE(dcpse, nullptr);

    int result = sv_dcpse_update(dcpse);
    EXPECT_EQ(result, 0);  // Success

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== Derivative Application Tests ====================

TEST_F(LispDcpseTest, ApplyDx2D) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 20);
    ASSERT_NE(dcpse, nullptr);

    sv_dcpse_update(dcpse);

    int n = sv_dcpse_num_nodes(dcpse);
    std::vector<double> input(n, 1.0);  // Constant function
    std::vector<double> output(n);

    sv_dcpse_apply_dx(dcpse, input.data(), output.data(), n);

    // Derivative of constant should be approximately zero
    double max_abs = 0.0;
    for (int i = 0; i < n; ++i) {
        max_abs = std::max(max_abs, std::abs(output[i]));
    }
    EXPECT_LT(max_abs, 1e-6);

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

TEST_F(LispDcpseTest, ApplyDy2D) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 20);
    ASSERT_NE(dcpse, nullptr);

    sv_dcpse_update(dcpse);

    int n = sv_dcpse_num_nodes(dcpse);
    std::vector<double> input(n, 1.0);
    std::vector<double> output(n);

    sv_dcpse_apply_dy(dcpse, input.data(), output.data(), n);

    // Derivative of constant should be approximately zero
    double max_abs = 0.0;
    for (int i = 0; i < n; ++i) {
        max_abs = std::max(max_abs, std::abs(output[i]));
    }
    EXPECT_LT(max_abs, 1e-6);

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

TEST_F(LispDcpseTest, ApplyLaplacian2D) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 20);
    ASSERT_NE(dcpse, nullptr);

    sv_dcpse_update(dcpse);

    int n = sv_dcpse_num_nodes(dcpse);
    std::vector<double> input(n, 1.0);  // Constant function
    std::vector<double> output(n);

    sv_dcpse_apply_laplacian(dcpse, input.data(), output.data(), n);

    // Laplacian of constant should be approximately zero
    double max_abs = 0.0;
    for (int i = 0; i < n; ++i) {
        max_abs = std::max(max_abs, std::abs(output[i]));
    }
    EXPECT_LT(max_abs, 1e-6);

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== Matrix Info Tests ====================

TEST_F(LispDcpseTest, MatrixInfo2D) {
    void* mesh = sv_make_cartesian_mesh_2d(3, 3, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 15);
    ASSERT_NE(dcpse, nullptr);

    sv_dcpse_update(dcpse);

    int num_rows, num_cols, nnz;
    sv_dcpse_matrix_info(dcpse, "dx", &num_rows, &num_cols, &nnz);

    EXPECT_EQ(num_rows, 16);  // (3+1)^2 = 16
    EXPECT_EQ(num_cols, 16);
    EXPECT_GT(nnz, 0);

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== Neighbor Tests ====================

TEST_F(LispDcpseTest, GetNodeNeighbors) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 3, 1.0, 1.0);
    ASSERT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    ASSERT_NE(gf, nullptr);

    void* dcpse = sv_make_dcpse_2d(gf, 15);
    ASSERT_NE(dcpse, nullptr);

    sv_dcpse_update(dcpse);

    // Get neighbors for node 0
    int neighbors[30];
    int n;
    sv_dcpse_get_node_neighbors(dcpse, 0, neighbors, &n);

    EXPECT_GT(n, 0);
    EXPECT_LE(n, 15);  // Should not exceed num_neighbors

    // All neighbor indices should be valid
    int num_nodes = sv_dcpse_num_nodes(dcpse);
    for (int i = 0; i < n; ++i) {
        EXPECT_GE(neighbors[i], 0);
        EXPECT_LT(neighbors[i], num_nodes);
    }

    sv_free_dcpse(dcpse);
    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== DcpseWrapper C++ API Tests ====================

TEST_F(LispDcpseTest, DcpseWrapperMake2D) {
    auto mesh = std::make_unique<mfem::Mesh>(4, 4, mfem::Element::QUADRILATERAL,
                                              false, 1.0, 1.0, false);
    mfem::H1_FECollection fec(1, 2);
    mfem::FiniteElementSpace fes(mesh.get(), &fec, 1);
    mfem::GridFunction gf(&fes);

    auto dcpse = StreamVorti::Lisp::DcpseWrapper::make2D(gf, 20);
    EXPECT_NE(dcpse, nullptr);

    dcpse->Update();

    // Test derivative matrix access
    auto dx_mat = StreamVorti::Lisp::DcpseWrapper::getDerivativeMatrix(
        dcpse.get(), "dx");
    EXPECT_NE(dx_mat, nullptr);
    EXPECT_GT(dx_mat->Height(), 0);
}

#endif // STREAMVORTI_WITH_ECL

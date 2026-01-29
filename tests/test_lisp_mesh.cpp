/**
 * @file test_lisp_mesh.cpp
 * @brief Unit tests for MFEM mesh bindings to Lisp
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 */

#include <gtest/gtest.h>

#ifdef STREAMVORTI_WITH_ECL

#include <StreamVorti/lisp/ecl_runtime.hpp>
#include <StreamVorti/lisp/lisp_mesh.hpp>

#include "mfem.hpp"

class LispMeshTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!StreamVorti::Lisp::Runtime::isInitialized()) {
            StreamVorti::Lisp::Runtime::init();
        }
    }
};

// ==================== Mesh Creation Tests ====================

TEST_F(LispMeshTest, CreateCartesian2DQuad) {
    void* mesh = sv_make_cartesian_mesh_2d(5, 5, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    EXPECT_EQ(sv_mesh_dimension(mesh), 2);
    EXPECT_EQ(sv_mesh_space_dimension(mesh), 2);
    EXPECT_EQ(sv_mesh_num_vertices(mesh), 36);  // (5+1)*(5+1) = 36
    EXPECT_EQ(sv_mesh_num_elements(mesh), 25);  // 5*5 = 25

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, CreateCartesian2DTri) {
    void* mesh = sv_make_cartesian_mesh_2d(4, 4, 2, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    EXPECT_EQ(sv_mesh_dimension(mesh), 2);
    EXPECT_EQ(sv_mesh_num_vertices(mesh), 25);  // (4+1)*(4+1) = 25
    EXPECT_EQ(sv_mesh_num_elements(mesh), 32);  // 4*4*2 = 32 triangles

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, CreateCartesian3DHex) {
    void* mesh = sv_make_cartesian_mesh_3d(3, 3, 3, 5, 1.0, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    EXPECT_EQ(sv_mesh_dimension(mesh), 3);
    EXPECT_EQ(sv_mesh_space_dimension(mesh), 3);
    EXPECT_EQ(sv_mesh_num_vertices(mesh), 64);  // (3+1)^3 = 64
    EXPECT_EQ(sv_mesh_num_elements(mesh), 27);  // 3^3 = 27

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, CreateCartesian3DTet) {
    void* mesh = sv_make_cartesian_mesh_3d(2, 2, 2, 4, 1.0, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    EXPECT_EQ(sv_mesh_dimension(mesh), 3);
    EXPECT_EQ(sv_mesh_num_vertices(mesh), 27);  // (2+1)^3 = 27

    sv_free_mesh(mesh);
}

// ==================== Mesh Query Tests ====================

TEST_F(LispMeshTest, GetVertexCoordinates) {
    void* mesh = sv_make_cartesian_mesh_2d(2, 2, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    // Get all vertices
    int nv = sv_mesh_num_vertices(mesh);
    int dim = sv_mesh_space_dimension(mesh);
    EXPECT_EQ(nv, 9);
    EXPECT_EQ(dim, 2);

    std::vector<double> coords(nv * dim);
    int n;
    sv_mesh_get_vertices(mesh, coords.data(), &n);
    EXPECT_EQ(n, nv * dim);

    // First vertex should be at origin (0,0)
    // Note: MFEM ordering may vary
    double vertex0[2];
    sv_mesh_get_vertex(mesh, 0, vertex0);
    EXPECT_GE(vertex0[0], 0.0);
    EXPECT_LE(vertex0[0], 1.0);
    EXPECT_GE(vertex0[1], 0.0);
    EXPECT_LE(vertex0[1], 1.0);

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, GetElementVertices) {
    void* mesh = sv_make_cartesian_mesh_2d(2, 2, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    int vertices[4];
    int n;
    sv_mesh_get_element_vertices(mesh, 0, vertices, &n);

    EXPECT_EQ(n, 4);  // Quad has 4 vertices
    for (int i = 0; i < n; ++i) {
        EXPECT_GE(vertices[i], 0);
        EXPECT_LT(vertices[i], sv_mesh_num_vertices(mesh));
    }

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, GetBoundaryElements) {
    void* mesh = sv_make_cartesian_mesh_2d(3, 3, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    int nbe = sv_mesh_num_boundary_elements(mesh);
    EXPECT_EQ(nbe, 12);  // 3*4 = 12 boundary edges for 3x3 mesh

    sv_free_mesh(mesh);
}

// ==================== Boundary Attribute Tests ====================

TEST_F(LispMeshTest, SetBoundaryAttributeX0) {
    void* mesh = sv_make_cartesian_mesh_2d(3, 3, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    // Set attribute for x=0 boundary
    int count = sv_mesh_set_boundary_attribute(mesh, 0, 0.0, 1e-10, 10);
    EXPECT_EQ(count, 3);  // 3 boundary elements at x=0

    // Verify attribute was set
    bool found = false;
    for (int i = 0; i < sv_mesh_num_boundary_elements(mesh); ++i) {
        if (sv_mesh_get_boundary_attribute(mesh, i) == 10) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);

    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, SetBoundaryAttributeY1) {
    void* mesh = sv_make_cartesian_mesh_2d(3, 3, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    // Set attribute for y=1 boundary (top)
    int count = sv_mesh_set_boundary_attribute(mesh, 1, 1.0, 1e-10, 20);
    EXPECT_EQ(count, 3);  // 3 boundary elements at y=1

    sv_free_mesh(mesh);
}

// ==================== Mesh Refinement Tests ====================

TEST_F(LispMeshTest, UniformRefinement) {
    void* mesh = sv_make_cartesian_mesh_2d(2, 2, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    int initial_vertices = sv_mesh_num_vertices(mesh);
    int initial_elements = sv_mesh_num_elements(mesh);

    sv_mesh_refine(mesh, 1);

    // After one refinement, element count quadruples for quads
    EXPECT_GT(sv_mesh_num_vertices(mesh), initial_vertices);
    EXPECT_EQ(sv_mesh_num_elements(mesh), initial_elements * 4);

    sv_free_mesh(mesh);
}

// ==================== GridFunction Tests ====================

TEST_F(LispMeshTest, CreateGridFunction) {
    void* mesh = sv_make_cartesian_mesh_2d(3, 3, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    EXPECT_NE(gf, nullptr);

    int size = sv_grid_function_size(gf);
    EXPECT_EQ(size, sv_mesh_num_vertices(mesh));

    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

TEST_F(LispMeshTest, GridFunctionSetGetValues) {
    void* mesh = sv_make_cartesian_mesh_2d(2, 2, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    void* gf = sv_make_grid_function(mesh, 1);
    EXPECT_NE(gf, nullptr);

    int size = sv_grid_function_size(gf);
    std::vector<double> values(size);

    // Set values to node indices
    for (int i = 0; i < size; ++i) {
        values[i] = static_cast<double>(i);
    }
    sv_grid_function_set_values(gf, values.data(), size);

    // Get values back
    std::vector<double> retrieved(size);
    int n;
    sv_grid_function_get_values(gf, retrieved.data(), &n);
    EXPECT_EQ(n, size);

    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(retrieved[i], values[i]);
    }

    sv_free_grid_function(gf);
    sv_free_mesh(mesh);
}

// ==================== MeshWrapper C++ API Tests ====================

TEST_F(LispMeshTest, MeshWrapperCartesian2D) {
    auto mesh = StreamVorti::Lisp::MeshWrapper::makeCartesian2D(4, 4, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 2);
    EXPECT_EQ(mesh->GetNV(), 25);
}

TEST_F(LispMeshTest, MeshWrapperCartesian3D) {
    auto mesh = StreamVorti::Lisp::MeshWrapper::makeCartesian3D(2, 2, 2, 5, 1.0, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
}

TEST_F(LispMeshTest, MeshWrapperBoundaryPredicate) {
    auto mesh = StreamVorti::Lisp::MeshWrapper::makeCartesian2D(3, 3, 3, 1.0, 1.0);
    EXPECT_NE(mesh, nullptr);

    // Find boundary elements at x=0
    auto pred = StreamVorti::Lisp::MeshWrapper::xEquals(0.0);
    auto elements = StreamVorti::Lisp::MeshWrapper::findBoundaryElements(mesh.get(), pred);

    EXPECT_EQ(elements.size(), 3);  // 3 elements at x=0 for 3x3 mesh
}

#endif // STREAMVORTI_WITH_ECL

/**
 * @file test_finite_differences.cpp
 * @brief Unit tests for finite difference derivative matrices
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
 */

#include <gtest/gtest.h>
#include "mfem.hpp"
#include "StreamVorti/finite_differences/fd_2d.hpp"
#include "StreamVorti/finite_differences/fd_3d.hpp"

#include <cmath>

// ---------------------------------------------------------------------------
// 2D Tests
// ---------------------------------------------------------------------------

class FiniteDiff2dTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a 10x10 element Cartesian mesh on [0,1]^2
        // This gives 11x11 = 121 nodes
        mesh = new mfem::Mesh(
            mfem::Mesh::MakeCartesian2D(10, 10,
                mfem::Element::QUADRILATERAL, false, 1.0, 1.0, false));
        fec = new mfem::H1_FECollection(1, 2);
        fes = new mfem::FiniteElementSpace(mesh, fec, 1);
        gf = new mfem::GridFunction(fes);
        *gf = 0.0;
    }

    void TearDown() override {
        delete gf;
        delete fes;
        delete fec;
        delete mesh;
    }

    mfem::Mesh *mesh;
    mfem::H1_FECollection *fec;
    mfem::FiniteElementSpace *fes;
    mfem::GridFunction *gf;
};

// Test that FD matrices have correct dimensions
TEST_F(FiniteDiff2dTest, MatrixDimensions) {
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    int n = fes->GetNDofs();
    EXPECT_EQ(fd.ShapeFunctionDx().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDx().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDy().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDy().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDxx().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDxx().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDyy().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDyy().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDxy().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDxy().Width(), n);
}

// Test D(i) accessor matches named accessors
TEST_F(FiniteDiff2dTest, DAccessor) {
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    // D(0) should be the same pointer as ShapeFunctionDx()
    EXPECT_EQ(&fd.D(0), &fd.ShapeFunctionDx());
    EXPECT_EQ(&fd.D(1), &fd.ShapeFunctionDy());
}

// Test dx of a linear function f(x,y) = x
// df/dx = 1 everywhere
TEST_F(FiniteDiff2dTest, DxOfLinearX) {
    // Set gf = x coordinate
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(dudx(i), 1.0, 1e-10) << "at node " << i;
    }
}

// Test dy of a linear function f(x,y) = y
// df/dy = 1 everywhere
TEST_F(FiniteDiff2dTest, DyOfLinearY) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[1];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudy(fes->GetNDofs());
    fd.ShapeFunctionDy().Mult(*gf, dudy);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(dudy(i), 1.0, 1e-10) << "at node " << i;
    }
}

// Test dx of f(x,y) = x² should give 2x
TEST_F(FiniteDiff2dTest, DxOfQuadraticX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[0];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);

    // Central differences are exact for quadratics on interior nodes.
    // Boundary nodes use 2nd-order one-sided stencils, also exact for quadratics.
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        double expected = 2.0 * vtx[0];
        EXPECT_NEAR(dudx(i), expected, 1e-10) << "at node " << i;
    }
}

// Test dxx of f(x,y) = x² should give 2 everywhere
TEST_F(FiniteDiff2dTest, DxxOfQuadraticX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[0];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udx2(fes->GetNDofs());
    fd.ShapeFunctionDxx().Mult(*gf, d2udx2);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(d2udx2(i), 2.0, 1e-10) << "at node " << i;
    }
}

// Test dyy of f(x,y) = y² should give 2 everywhere
TEST_F(FiniteDiff2dTest, DyyOfQuadraticY) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[1] * vtx[1];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udy2(fes->GetNDofs());
    fd.ShapeFunctionDyy().Mult(*gf, d2udy2);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(d2udy2(i), 2.0, 1e-10) << "at node " << i;
    }
}

// Test dxy of f(x,y) = x*y should give 1 at interior nodes
TEST_F(FiniteDiff2dTest, DxyOfBilinear) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[1];
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udxdy(fes->GetNDofs());
    fd.ShapeFunctionDxy().Mult(*gf, d2udxdy);

    // Check interior nodes (not on boundary)
    double h = 0.1;
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        bool on_boundary = (vtx[0] < 1e-12 || vtx[0] > 1.0 - 1e-12 ||
                           vtx[1] < 1e-12 || vtx[1] > 1.0 - 1e-12);
        if (!on_boundary) {
            EXPECT_NEAR(d2udxdy(i), 1.0, 1e-10) << "at interior node " << i;
        }
    }
}

// Test that constant function has zero derivatives
TEST_F(FiniteDiff2dTest, ConstantFunctionZeroDerivatives) {
    *gf = 5.0;

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector result(fes->GetNDofs());

    fd.ShapeFunctionDx().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDy().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDxx().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDyy().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);
}

// Test Laplacian: dxx + dyy of f(x,y) = sin(pi*x)*sin(pi*y)
// Laplacian = -2*pi²*sin(pi*x)*sin(pi*y)
TEST_F(FiniteDiff2dTest, LaplacianOfSinusoidal) {
    const double pi = M_PI;
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = std::sin(pi * vtx[0]) * std::sin(pi * vtx[1]);
    }

    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dxx(fes->GetNDofs()), dyy(fes->GetNDofs());
    fd.ShapeFunctionDxx().Mult(*gf, dxx);
    fd.ShapeFunctionDyy().Mult(*gf, dyy);

    // Check interior nodes only (boundary stencils have different truncation error)
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        bool interior = (vtx[0] > 0.15 && vtx[0] < 0.85 &&
                        vtx[1] > 0.15 && vtx[1] < 0.85);
        if (interior) {
            double expected = -2.0 * pi * pi * std::sin(pi * vtx[0]) * std::sin(pi * vtx[1]);
            double computed = dxx(i) + dyy(i);
            // 2nd-order FD has O(h²) truncation error, h=0.1 so ~0.01*pi² ≈ 0.1
            EXPECT_NEAR(computed, expected, 0.15)
                << "at (" << vtx[0] << ", " << vtx[1] << ")";
        }
    }
}

// ---------------------------------------------------------------------------
// 3D Tests
// ---------------------------------------------------------------------------

class FiniteDiff3dTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a 5x5x5 element Cartesian mesh on [0,1]^3
        // This gives 6x6x6 = 216 nodes
        mesh = new mfem::Mesh(
            mfem::Mesh::MakeCartesian3D(5, 5, 5,
                mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
        fec = new mfem::H1_FECollection(1, 3);
        fes = new mfem::FiniteElementSpace(mesh, fec, 1);
        gf = new mfem::GridFunction(fes);
        *gf = 0.0;
    }

    void TearDown() override {
        delete gf;
        delete fes;
        delete fec;
        delete mesh;
    }

    mfem::Mesh *mesh;
    mfem::H1_FECollection *fec;
    mfem::FiniteElementSpace *fes;
    mfem::GridFunction *gf;
};

// Test 3D matrix dimensions
TEST_F(FiniteDiff3dTest, MatrixDimensions) {
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    int n = fes->GetNDofs();
    EXPECT_EQ(fd.ShapeFunctionDx().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDx().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDz().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDzz().Height(), n);
}

// Test dx of f(x,y,z) = x in 3D
TEST_F(FiniteDiff3dTest, DxOfLinearX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0];
    }

    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(dudx(i), 1.0, 1e-10) << "at node " << i;
    }
}

// Test dz of f(x,y,z) = z in 3D
TEST_F(FiniteDiff3dTest, DzOfLinearZ) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[2];
    }

    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector dudz(fes->GetNDofs());
    fd.ShapeFunctionDz().Mult(*gf, dudz);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(dudz(i), 1.0, 1e-10) << "at node " << i;
    }
}

// Test dzz of f(x,y,z) = z² should give 2
TEST_F(FiniteDiff3dTest, DzzOfQuadraticZ) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[2] * vtx[2];
    }

    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector d2udz2(fes->GetNDofs());
    fd.ShapeFunctionDzz().Mult(*gf, d2udz2);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        EXPECT_NEAR(d2udz2(i), 2.0, 1e-10) << "at node " << i;
    }
}

// Test 3D constant function has zero derivatives
TEST_F(FiniteDiff3dTest, ConstantFunctionZeroDerivatives) {
    *gf = 3.0;

    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector result(fes->GetNDofs());

    fd.ShapeFunctionDx().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDz().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDzz().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);
}

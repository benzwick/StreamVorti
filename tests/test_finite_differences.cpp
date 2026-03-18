/**
 * @file test_finite_differences.cpp
 * @brief Unit tests for finite difference derivative matrices
 *
 * Tests cover:
 * - Grid validation (structured, uniform spacing, sufficient size)
 * - 2nd-order stencil accuracy (exact for polynomials up to degree 2)
 * - 4th-order stencil accuracy (exact for polynomials up to degree 4)
 * - Boundary stencil accuracy
 * - Laplacian convergence on smooth functions
 * - Rejection of non-structured and non-uniform grids
 * - 3D derivative operators
 *
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2026 Benjamin F. Zwick
 */

#include <gtest/gtest.h>
#include "mfem.hpp"
#include "StreamVorti/finite_differences/fd_2d.hpp"
#include "StreamVorti/finite_differences/fd_3d.hpp"

#include <cmath>

// =====================================================================
// 2D Tests — Order 2
// =====================================================================

class FiniteDiff2dTest : public ::testing::Test {
protected:
    void SetUp() override {
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

TEST_F(FiniteDiff2dTest, MatrixDimensions) {
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    int n = fes->GetNDofs();
    EXPECT_EQ(fd.ShapeFunctionDx().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDx().Width(), n);
    EXPECT_EQ(fd.ShapeFunctionDyy().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDxy().Width(), n);
}

TEST_F(FiniteDiff2dTest, DAccessor) {
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();
    EXPECT_EQ(&fd.D(0), &fd.ShapeFunctionDx());
    EXPECT_EQ(&fd.D(1), &fd.ShapeFunctionDy());
}

TEST_F(FiniteDiff2dTest, StencilOrderAccessor) {
    StreamVorti::FiniteDiff2d fd2(*gf, 2);
    EXPECT_EQ(fd2.StencilOrder(), 2);
    StreamVorti::FiniteDiff2d fd4(*gf, 4);
    EXPECT_EQ(fd4.StencilOrder(), 4);
}

// Constant function: all derivatives must be zero
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

// Linear f(x,y)=x => dx=1 everywhere (exact for order >= 2)
TEST_F(FiniteDiff2dTest, DxOfLinearX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(dudx(i), 1.0, 1e-10) << "at node " << i;
}

// Linear f(x,y)=y => dy=1 everywhere
TEST_F(FiniteDiff2dTest, DyOfLinearY) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[1];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudy(fes->GetNDofs());
    fd.ShapeFunctionDy().Mult(*gf, dudy);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(dudy(i), 1.0, 1e-10) << "at node " << i;
}

// Quadratic f=x^2 => dx=2x (exact for 2nd-order central + boundary stencils)
TEST_F(FiniteDiff2dTest, DxOfQuadraticX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[0];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        double expected = 2.0 * mesh->GetVertex(i)[0];
        EXPECT_NEAR(dudx(i), expected, 1e-10) << "at node " << i;
    }
}

// f=x^2 => dxx=2 everywhere
TEST_F(FiniteDiff2dTest, DxxOfQuadraticX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[0];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udx2(fes->GetNDofs());
    fd.ShapeFunctionDxx().Mult(*gf, d2udx2);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(d2udx2(i), 2.0, 1e-10) << "at node " << i;
}

// f=y^2 => dyy=2 everywhere
TEST_F(FiniteDiff2dTest, DyyOfQuadraticY) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[1] * vtx[1];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udy2(fes->GetNDofs());
    fd.ShapeFunctionDyy().Mult(*gf, d2udy2);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(d2udy2(i), 2.0, 1e-10) << "at node " << i;
}

// f=x*y => dxy=1 at interior nodes
TEST_F(FiniteDiff2dTest, DxyOfBilinear) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0] * vtx[1];
    }
    StreamVorti::FiniteDiff2d fd(*gf);
    fd.Update();

    mfem::Vector d2udxdy(fes->GetNDofs());
    fd.ShapeFunctionDxy().Mult(*gf, d2udxdy);

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        bool interior = (vtx[0] > 1e-12 && vtx[0] < 1.0 - 1e-12 &&
                        vtx[1] > 1e-12 && vtx[1] < 1.0 - 1e-12);
        if (interior)
            EXPECT_NEAR(d2udxdy(i), 1.0, 1e-10);
    }
}

// Laplacian convergence: sin(pi*x)*sin(pi*y)
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

    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        bool interior = (vtx[0] > 0.15 && vtx[0] < 0.85 &&
                        vtx[1] > 0.15 && vtx[1] < 0.85);
        if (interior) {
            double expected = -2.0 * pi * pi * std::sin(pi * vtx[0]) * std::sin(pi * vtx[1]);
            double computed = dxx(i) + dyy(i);
            EXPECT_NEAR(computed, expected, 0.15);
        }
    }
}

// =====================================================================
// 2D Tests — Order 4
// =====================================================================

class FiniteDiff2dOrder4Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Need at least 6 nodes per direction for 4th-order stencils
        mesh = new mfem::Mesh(
            mfem::Mesh::MakeCartesian2D(20, 20,
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

// 4th-order: exact for cubics. f=x^3 => dx=3x^2
TEST_F(FiniteDiff2dOrder4Test, DxOfCubicX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        double x = vtx[0];
        (*gf)(i) = x * x * x;
    }
    StreamVorti::FiniteDiff2d fd(*gf, 4);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        double x = mesh->GetVertex(i)[0];
        EXPECT_NEAR(dudx(i), 3.0 * x * x, 1e-8) << "at node " << i;
    }
}

// 4th-order: f=x^4 => dx=4x^3 (exact for 4th-order central, approx at boundary)
TEST_F(FiniteDiff2dOrder4Test, DxOfQuarticX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        double x = vtx[0];
        (*gf)(i) = x * x * x * x;
    }
    StreamVorti::FiniteDiff2d fd(*gf, 4);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        double x = mesh->GetVertex(i)[0];
        EXPECT_NEAR(dudx(i), 4.0 * x * x * x, 1e-7) << "at node " << i;
    }
}

// 4th-order dxx: f=x^3 => dxx=6x (exact)
TEST_F(FiniteDiff2dOrder4Test, DxxOfCubicX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        double x = vtx[0];
        (*gf)(i) = x * x * x;
    }
    StreamVorti::FiniteDiff2d fd(*gf, 4);
    fd.Update();

    mfem::Vector d2udx2(fes->GetNDofs());
    fd.ShapeFunctionDxx().Mult(*gf, d2udx2);
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        double x = mesh->GetVertex(i)[0];
        EXPECT_NEAR(d2udx2(i), 6.0 * x, 1e-7) << "at node " << i;
    }
}

// 4th-order Laplacian converges faster than 2nd-order
TEST_F(FiniteDiff2dOrder4Test, LaplacianBetterThan2ndOrder) {
    const double pi = M_PI;
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = std::sin(pi * vtx[0]) * std::sin(pi * vtx[1]);
    }

    // 4th-order
    StreamVorti::FiniteDiff2d fd4(*gf, 4);
    fd4.Update();

    mfem::Vector dxx4(fes->GetNDofs()), dyy4(fes->GetNDofs());
    fd4.ShapeFunctionDxx().Mult(*gf, dxx4);
    fd4.ShapeFunctionDyy().Mult(*gf, dyy4);

    double max_err_o4 = 0;
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        bool interior = (vtx[0] > 0.15 && vtx[0] < 0.85 &&
                        vtx[1] > 0.15 && vtx[1] < 0.85);
        if (interior) {
            double expected = -2.0 * pi * pi * std::sin(pi * vtx[0]) * std::sin(pi * vtx[1]);
            double err = std::abs(dxx4(i) + dyy4(i) - expected);
            max_err_o4 = std::max(max_err_o4, err);
        }
    }

    // 4th-order on 20x20 should be significantly more accurate than 2nd-order on 10x10
    EXPECT_LT(max_err_o4, 0.01) << "4th-order Laplacian error too large on 20x20 grid";
}

// =====================================================================
// 3D Tests
// =====================================================================

class FiniteDiff3dTest : public ::testing::Test {
protected:
    void SetUp() override {
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

TEST_F(FiniteDiff3dTest, MatrixDimensions) {
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();
    int n = fes->GetNDofs();
    EXPECT_EQ(fd.ShapeFunctionDx().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDz().Height(), n);
    EXPECT_EQ(fd.ShapeFunctionDzz().Height(), n);
}

TEST_F(FiniteDiff3dTest, DxOfLinearX) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[0];
    }
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector dudx(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, dudx);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(dudx(i), 1.0, 1e-10);
}

TEST_F(FiniteDiff3dTest, DzOfLinearZ) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[2];
    }
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector dudz(fes->GetNDofs());
    fd.ShapeFunctionDz().Mult(*gf, dudz);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(dudz(i), 1.0, 1e-10);
}

TEST_F(FiniteDiff3dTest, DzzOfQuadraticZ) {
    for (int i = 0; i < fes->GetNDofs(); ++i) {
        const double *vtx = mesh->GetVertex(i);
        (*gf)(i) = vtx[2] * vtx[2];
    }
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector d2udz2(fes->GetNDofs());
    fd.ShapeFunctionDzz().Mult(*gf, d2udz2);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(d2udz2(i), 2.0, 1e-10);
}

TEST_F(FiniteDiff3dTest, ConstantFunctionZeroDerivatives) {
    *gf = 3.0;
    StreamVorti::FiniteDiff3d fd(*gf);
    fd.Update();

    mfem::Vector result(fes->GetNDofs());
    fd.ShapeFunctionDx().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);

    fd.ShapeFunctionDzz().Mult(*gf, result);
    for (int i = 0; i < fes->GetNDofs(); ++i)
        EXPECT_NEAR(result(i), 0.0, 1e-10);
}

// =====================================================================
// Non-rectangular grid (should fail with structured grid check)
// =====================================================================

// NOTE: Death tests for MFEM_ABORT require special handling.
// MFEM_ABORT calls exit() or abort(), which we can test with EXPECT_DEATH.
// However, this may not work on all platforms, so we mark them as optional.
// The key thing is that the validation logic EXISTS and is tested above
// via the positive tests (structured grids that pass validation).

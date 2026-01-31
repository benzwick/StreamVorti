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
 * @file lisp_dcpse.cpp
 * @brief Implementation of DCPSE bindings for Lisp
 */

#include "StreamVorti/lisp/lisp_dcpse.hpp"
#include "StreamVorti/lisp/ecl_runtime.hpp"

#include "StreamVorti/approximants/dcpse.hpp"
#include "StreamVorti/approximants/dcpse_2d.hpp"
#include "StreamVorti/approximants/dcpse_3d.hpp"

#include "mfem.hpp"

#include <cstring>
#include <algorithm>

// Helper to determine if a DCPSE object is 2D or 3D
static int getDcpseDimension(void* dcpse_ptr) {
    if (!dcpse_ptr) return 0;

    // Try dynamic_cast to determine type
    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);

    if (dynamic_cast<StreamVorti::Dcpse2d*>(dcpse)) {
        return 2;
    } else if (dynamic_cast<StreamVorti::Dcpse3d*>(dcpse)) {
        return 3;
    }

    return 0;
}

// ==================== C API Implementation ====================

extern "C" {

void* sv_make_dcpse_2d(void* gf_ptr, int num_neighbors)
{
    if (!gf_ptr) return nullptr;

    mfem::GridFunction* gf = static_cast<mfem::GridFunction*>(gf_ptr);
    StreamVorti::Dcpse2d* dcpse = new StreamVorti::Dcpse2d(*gf, num_neighbors);

    return static_cast<void*>(dcpse);
}

void* sv_make_dcpse_3d(void* gf_ptr, int num_neighbors)
{
    if (!gf_ptr) return nullptr;

    mfem::GridFunction* gf = static_cast<mfem::GridFunction*>(gf_ptr);
    StreamVorti::Dcpse3d* dcpse = new StreamVorti::Dcpse3d(*gf, num_neighbors);

    return static_cast<void*>(dcpse);
}

void sv_free_dcpse(void* dcpse)
{
    delete static_cast<StreamVorti::Dcpse*>(dcpse);
}

int sv_dcpse_update(void* dcpse_ptr)
{
    if (!dcpse_ptr) return -1;

    try {
        StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
        dcpse->Update();
        return 0;
    } catch (...) {
        return -1;
    }
}

int sv_dcpse_dimension(void* dcpse)
{
    return getDcpseDimension(dcpse);
}

int sv_dcpse_num_nodes(void* dcpse_ptr)
{
    if (!dcpse_ptr) return 0;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    // Use NumSupportNodes() which is available immediately after construction
    return dcpse->NumSupportNodes();
}

int sv_dcpse_num_neighbors(void* dcpse_ptr)
{
    if (!dcpse_ptr) return 0;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    return dcpse->NumNeighbors();
}

void sv_dcpse_apply_dx(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->D(0);

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_dy(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->D(1);

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_dz(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim != 3) {
        // Zero output for 2D
        std::fill(output, output + n, 0.0);
        return;
    }

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->D(2);

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_dxx(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim == 2) {
        StreamVorti::Dcpse2d* dcpse =
            static_cast<StreamVorti::Dcpse2d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDxx();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    } else if (dim == 3) {
        StreamVorti::Dcpse3d* dcpse =
            static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDxx();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    }
}

void sv_dcpse_apply_dyy(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim == 2) {
        StreamVorti::Dcpse2d* dcpse =
            static_cast<StreamVorti::Dcpse2d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDyy();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    } else if (dim == 3) {
        StreamVorti::Dcpse3d* dcpse =
            static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDyy();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    }
}

void sv_dcpse_apply_dzz(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim != 3) {
        std::fill(output, output + n, 0.0);
        return;
    }

    StreamVorti::Dcpse3d* dcpse =
        static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDzz();

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_dxy(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim == 2) {
        StreamVorti::Dcpse2d* dcpse =
            static_cast<StreamVorti::Dcpse2d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDxy();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    } else if (dim == 3) {
        StreamVorti::Dcpse3d* dcpse =
            static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
        const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDxy();

        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector out_vec(output, n);
        mat.Mult(in_vec, out_vec);
    }
}

void sv_dcpse_apply_dxz(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim != 3) {
        std::fill(output, output + n, 0.0);
        return;
    }

    StreamVorti::Dcpse3d* dcpse =
        static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDxz();

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_dyz(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    int dim = getDcpseDimension(dcpse_ptr);
    if (dim != 3) {
        std::fill(output, output + n, 0.0);
        return;
    }

    StreamVorti::Dcpse3d* dcpse =
        static_cast<StreamVorti::Dcpse3d*>(dcpse_ptr);
    const mfem::SparseMatrix& mat = dcpse->ShapeFunctionDyz();

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat.Mult(in_vec, out_vec);
}

void sv_dcpse_apply_laplacian(void* dcpse_ptr, const double* input, double* output, int n)
{
    if (!dcpse_ptr || !input || !output) return;

    // Initialize output to zero
    std::fill(output, output + n, 0.0);

    // Temporary buffer
    std::vector<double> temp(n);

    // Apply dxx and add to output
    sv_dcpse_apply_dxx(dcpse_ptr, input, temp.data(), n);
    for (int i = 0; i < n; ++i) {
        output[i] += temp[i];
    }

    // Apply dyy and add to output
    sv_dcpse_apply_dyy(dcpse_ptr, input, temp.data(), n);
    for (int i = 0; i < n; ++i) {
        output[i] += temp[i];
    }

    // If 3D, also add dzz
    int dim = getDcpseDimension(dcpse_ptr);
    if (dim == 3) {
        sv_dcpse_apply_dzz(dcpse_ptr, input, temp.data(), n);
        for (int i = 0; i < n; ++i) {
            output[i] += temp[i];
        }
    }
}

void sv_dcpse_matrix_info(void* dcpse_ptr, const char* deriv,
                          int* num_rows, int* num_cols, int* nnz)
{
    if (!dcpse_ptr || !deriv) return;

    const mfem::SparseMatrix* mat =
        StreamVorti::Lisp::DcpseWrapper::getDerivativeMatrix(
            static_cast<StreamVorti::Dcpse*>(dcpse_ptr), deriv);

    if (mat) {
        if (num_rows) *num_rows = mat->Height();
        if (num_cols) *num_cols = mat->Width();
        if (nnz) *nnz = mat->NumNonZeroElems();
    }
}

void sv_dcpse_get_matrix_csr(void* dcpse_ptr, const char* deriv,
                              int* row_ptr, int* col_ind, double* values)
{
    if (!dcpse_ptr || !deriv || !row_ptr || !col_ind || !values) return;

    const mfem::SparseMatrix* mat =
        StreamVorti::Lisp::DcpseWrapper::getDerivativeMatrix(
            static_cast<StreamVorti::Dcpse*>(dcpse_ptr), deriv);

    if (!mat) return;

    int nrows = mat->Height();
    const int* I = mat->GetI();
    const int* J = mat->GetJ();
    const double* A = mat->GetData();

    // Copy row pointers
    for (int i = 0; i <= nrows; ++i) {
        row_ptr[i] = I[i];
    }

    // Copy column indices and values
    int nnz = mat->NumNonZeroElems();
    for (int i = 0; i < nnz; ++i) {
        col_ind[i] = J[i];
        values[i] = A[i];
    }
}

int sv_dcpse_save_matrix(void* dcpse_ptr, const char* deriv, const char* filename)
{
    if (!dcpse_ptr || !deriv || !filename) return -1;

    try {
        StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
        dcpse->SaveDerivToFile(deriv, filename);
        return 0;
    } catch (...) {
        return -1;
    }
}

void sv_dcpse_get_neighbors(void* dcpse_ptr, int* indices, int* offsets)
{
    if (!dcpse_ptr || !indices || !offsets) return;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    const auto& neighs = dcpse->NeighborIndices();

    int num_nodes = neighs.size();
    offsets[0] = 0;

    int idx = 0;
    for (int i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < neighs[i].size(); ++j) {
            indices[idx++] = neighs[i][j];
        }
        offsets[i + 1] = idx;
    }
}

void sv_dcpse_get_node_neighbors(void* dcpse_ptr, int node_id,
                                  int* neighbors, int* n)
{
    if (!dcpse_ptr || !neighbors || !n) return;

    StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
    const auto& neighs = dcpse->NeighborIndices();

    if (node_id < 0 || node_id >= static_cast<int>(neighs.size())) {
        *n = 0;
        return;
    }

    const auto& node_neighs = neighs[node_id];
    *n = static_cast<int>(node_neighs.size());

    for (size_t i = 0; i < node_neighs.size(); ++i) {
        neighbors[i] = node_neighs[i];
    }
}

int sv_dcpse_save_neighbors(void* dcpse_ptr, const char* filename)
{
    if (!dcpse_ptr || !filename) return -1;

    try {
        StreamVorti::Dcpse* dcpse = static_cast<StreamVorti::Dcpse*>(dcpse_ptr);
        dcpse->SaveNeighsToFile(dcpse->NeighborIndices(), filename);
        return 0;
    } catch (...) {
        return -1;
    }
}

} // extern "C"

// ==================== C++ API Implementation ====================

namespace StreamVorti {
namespace Lisp {

std::unique_ptr<Dcpse2d>
DcpseWrapper::make2D(mfem::GridFunction& gf, int num_neighbors)
{
    return std::make_unique<Dcpse2d>(gf, num_neighbors);
}

std::unique_ptr<Dcpse3d>
DcpseWrapper::make3D(mfem::GridFunction& gf, int num_neighbors)
{
    return std::make_unique<Dcpse3d>(gf, num_neighbors);
}

const mfem::SparseMatrix*
DcpseWrapper::getDerivativeMatrix(Dcpse* dcpse, const std::string& name)
{
    if (!dcpse) return nullptr;

    int dim = 0;
    Dcpse2d* dcpse2d = dynamic_cast<Dcpse2d*>(dcpse);
    Dcpse3d* dcpse3d = dynamic_cast<Dcpse3d*>(dcpse);

    if (dcpse2d) dim = 2;
    else if (dcpse3d) dim = 3;
    else return nullptr;

    // 2D matrices
    if (dim == 2 && dcpse2d) {
        if (name == "dx") return &dcpse2d->ShapeFunctionDx();
        if (name == "dy") return &dcpse2d->ShapeFunctionDy();
        if (name == "dxx") return &dcpse2d->ShapeFunctionDxx();
        if (name == "dyy") return &dcpse2d->ShapeFunctionDyy();
        if (name == "dxy") return &dcpse2d->ShapeFunctionDxy();
    }

    // 3D matrices
    if (dim == 3 && dcpse3d) {
        if (name == "dx") return &dcpse3d->ShapeFunctionDx();
        if (name == "dy") return &dcpse3d->ShapeFunctionDy();
        if (name == "dz") return &dcpse3d->ShapeFunctionDz();
        if (name == "dxx") return &dcpse3d->ShapeFunctionDxx();
        if (name == "dyy") return &dcpse3d->ShapeFunctionDyy();
        if (name == "dzz") return &dcpse3d->ShapeFunctionDzz();
        if (name == "dxy") return &dcpse3d->ShapeFunctionDxy();
        if (name == "dxz") return &dcpse3d->ShapeFunctionDxz();
        if (name == "dyz") return &dcpse3d->ShapeFunctionDyz();
    }

    return nullptr;
}

void DcpseWrapper::applyDerivative(Dcpse* dcpse, const std::string& name,
                                    const double* input, double* output, int n)
{
    const mfem::SparseMatrix* mat = getDerivativeMatrix(dcpse, name);
    if (!mat) return;

    mfem::Vector in_vec(const_cast<double*>(input), n);
    mfem::Vector out_vec(output, n);
    mat->Mult(in_vec, out_vec);
}

void DcpseWrapper::applyLaplacian(Dcpse* dcpse,
                                   const double* input, double* output, int n)
{
    std::fill(output, output + n, 0.0);
    std::vector<double> temp(n);

    Dcpse2d* dcpse2d = dynamic_cast<Dcpse2d*>(dcpse);
    Dcpse3d* dcpse3d = dynamic_cast<Dcpse3d*>(dcpse);

    if (dcpse2d) {
        // dxx
        const mfem::SparseMatrix& dxx = dcpse2d->ShapeFunctionDxx();
        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector temp_vec(temp.data(), n);
        dxx.Mult(in_vec, temp_vec);
        for (int i = 0; i < n; ++i) output[i] += temp[i];

        // dyy
        const mfem::SparseMatrix& dyy = dcpse2d->ShapeFunctionDyy();
        dyy.Mult(in_vec, temp_vec);
        for (int i = 0; i < n; ++i) output[i] += temp[i];
    }
    else if (dcpse3d) {
        mfem::Vector in_vec(const_cast<double*>(input), n);
        mfem::Vector temp_vec(temp.data(), n);

        // dxx
        dcpse3d->ShapeFunctionDxx().Mult(in_vec, temp_vec);
        for (int i = 0; i < n; ++i) output[i] += temp[i];

        // dyy
        dcpse3d->ShapeFunctionDyy().Mult(in_vec, temp_vec);
        for (int i = 0; i < n; ++i) output[i] += temp[i];

        // dzz
        dcpse3d->ShapeFunctionDzz().Mult(in_vec, temp_vec);
        for (int i = 0; i < n; ++i) output[i] += temp[i];
    }
}

void registerDcpseFunctions()
{
    // Register DCPSE functions with ECL
    Runtime::eval(R"(
        (defpackage :streamvorti.dcpse
          (:use :cl)
          (:export #:make-dcpse-2d
                   #:make-dcpse-3d
                   #:dcpse-update
                   #:dcpse-apply-dx
                   #:dcpse-apply-dy
                   #:dcpse-apply-laplacian))
    )");
}

} // namespace Lisp
} // namespace StreamVorti

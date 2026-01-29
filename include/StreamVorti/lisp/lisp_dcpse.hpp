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
 * @file lisp_dcpse.hpp
 * @brief DCPSE operations exposed to Lisp via CFFI
 *
 * Provides C-style functions for DCPSE derivative computation
 * that can be called from Common Lisp through CFFI.
 */

#ifndef STREAMVORTI_LISP_LISP_DCPSE_HPP_
#define STREAMVORTI_LISP_LISP_DCPSE_HPP_

#include <cstdint>

// C-style API for CFFI
#ifdef __cplusplus
extern "C" {
#endif

// ==================== DCPSE Creation ====================

/**
 * @brief Create 2D DCPSE derivative operator
 *
 * @param gf Pointer to GridFunction
 * @param num_neighbors Number of neighbors for support domain
 * @return Opaque pointer to Dcpse2d object
 */
void* sv_make_dcpse_2d(void* gf, int num_neighbors);

/**
 * @brief Create 3D DCPSE derivative operator
 *
 * @param gf Pointer to GridFunction
 * @param num_neighbors Number of neighbors for support domain
 * @return Opaque pointer to Dcpse3d object
 */
void* sv_make_dcpse_3d(void* gf, int num_neighbors);

/**
 * @brief Free DCPSE object memory
 */
void sv_free_dcpse(void* dcpse);

// ==================== DCPSE Operations ====================

/**
 * @brief Compute/update derivative matrices
 *
 * @param dcpse Pointer to DCPSE object
 * @return 0 on success, -1 on error
 */
int sv_dcpse_update(void* dcpse);

/**
 * @brief Get dimension of DCPSE (2 or 3)
 */
int sv_dcpse_dimension(void* dcpse);

/**
 * @brief Get number of nodes
 */
int sv_dcpse_num_nodes(void* dcpse);

/**
 * @brief Get number of neighbors per node
 */
int sv_dcpse_num_neighbors(void* dcpse);

// ==================== Derivative Application ====================

/**
 * @brief Apply dx operator to vector
 *
 * @param dcpse Pointer to DCPSE object
 * @param input Input values array
 * @param output Output values array (must be pre-allocated)
 * @param n Length of arrays
 */
void sv_dcpse_apply_dx(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dy operator to vector
 */
void sv_dcpse_apply_dy(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dz operator to vector (3D only)
 */
void sv_dcpse_apply_dz(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dxx operator to vector
 */
void sv_dcpse_apply_dxx(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dyy operator to vector
 */
void sv_dcpse_apply_dyy(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dzz operator to vector (3D only)
 */
void sv_dcpse_apply_dzz(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dxy operator to vector
 */
void sv_dcpse_apply_dxy(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dxz operator to vector (3D only)
 */
void sv_dcpse_apply_dxz(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply dyz operator to vector (3D only)
 */
void sv_dcpse_apply_dyz(void* dcpse, const double* input, double* output, int n);

/**
 * @brief Apply Laplacian operator (dxx + dyy [+ dzz]) to vector
 */
void sv_dcpse_apply_laplacian(void* dcpse, const double* input, double* output, int n);

// ==================== Matrix Access ====================

/**
 * @brief Get sparse matrix info for derivative operator
 *
 * @param dcpse Pointer to DCPSE object
 * @param deriv Derivative name: "dx", "dy", "dz", "dxx", etc.
 * @param num_rows Output: number of rows
 * @param num_cols Output: number of columns
 * @param nnz Output: number of non-zeros
 */
void sv_dcpse_matrix_info(void* dcpse, const char* deriv,
                          int* num_rows, int* num_cols, int* nnz);

/**
 * @brief Get sparse matrix in CSR format
 *
 * @param dcpse Pointer to DCPSE object
 * @param deriv Derivative name: "dx", "dy", etc.
 * @param row_ptr Output: row pointers (size = num_rows + 1)
 * @param col_ind Output: column indices (size = nnz)
 * @param values Output: matrix values (size = nnz)
 */
void sv_dcpse_get_matrix_csr(void* dcpse, const char* deriv,
                              int* row_ptr, int* col_ind, double* values);

/**
 * @brief Save derivative matrix to file
 *
 * @param dcpse Pointer to DCPSE object
 * @param deriv Derivative name
 * @param filename Output file path
 * @return 0 on success, -1 on error
 */
int sv_dcpse_save_matrix(void* dcpse, const char* deriv, const char* filename);

// ==================== Neighbor Information ====================

/**
 * @brief Get neighbor indices for all nodes
 *
 * @param dcpse Pointer to DCPSE object
 * @param indices Output: flat array of neighbor indices
 * @param offsets Output: offsets for each node (size = num_nodes + 1)
 */
void sv_dcpse_get_neighbors(void* dcpse, int* indices, int* offsets);

/**
 * @brief Get neighbor indices for a single node
 *
 * @param dcpse Pointer to DCPSE object
 * @param node_id Node index
 * @param neighbors Output: array of neighbor indices
 * @param n Output: number of neighbors
 */
void sv_dcpse_get_node_neighbors(void* dcpse, int node_id,
                                  int* neighbors, int* n);

/**
 * @brief Save neighbor information to file
 *
 * @param dcpse Pointer to DCPSE object
 * @param filename Output file path
 * @return 0 on success, -1 on error
 */
int sv_dcpse_save_neighbors(void* dcpse, const char* filename);

#ifdef __cplusplus
}
#endif

// C++ API
#ifdef __cplusplus

#include <memory>
#include <string>
#include <vector>

namespace mfem {
    class GridFunction;
    class SparseMatrix;
}

namespace StreamVorti {
    class Dcpse;
    class Dcpse2d;
    class Dcpse3d;

namespace Lisp {

/**
 * @class DcpseWrapper
 * @brief C++ wrapper for DCPSE with Lisp integration
 */
class DcpseWrapper {
public:
    /**
     * @brief Create 2D DCPSE
     */
    static std::unique_ptr<Dcpse2d>
    make2D(mfem::GridFunction& gf, int num_neighbors);

    /**
     * @brief Create 3D DCPSE
     */
    static std::unique_ptr<Dcpse3d>
    make3D(mfem::GridFunction& gf, int num_neighbors);

    /**
     * @brief Get derivative matrix by name
     */
    static const mfem::SparseMatrix*
    getDerivativeMatrix(Dcpse* dcpse, const std::string& name);

    /**
     * @brief Apply derivative operator
     */
    static void applyDerivative(Dcpse* dcpse, const std::string& name,
                                const double* input, double* output, int n);

    /**
     * @brief Apply Laplacian (sum of second derivatives)
     */
    static void applyLaplacian(Dcpse* dcpse,
                               const double* input, double* output, int n);
};

/**
 * @brief Register DCPSE functions with ECL
 */
void registerDcpseFunctions();

} // namespace Lisp
} // namespace StreamVorti

#endif // __cplusplus

#endif // STREAMVORTI_LISP_LISP_DCPSE_HPP_

/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017 Konstantinos A. Mountris
 * Copyright (C) 2020-2025 Benjamin F. Zwick
 * Copyright (C) 2026 Weizheng (Will) Li
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
 *
 * Contributors (alphabetically):
 *      George C. BOURANTAS
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

#ifndef STREAMVORTI_APPROXIMANTS_PAR_DCPSE_2D_HPP_
#define STREAMVORTI_APPROXIMANTS_PAR_DCPSE_2D_HPP_

#include "StreamVorti/approximants/par_dcpse.hpp"

#ifdef MFEM_USE_MPI

#include "mfem.hpp"

namespace StreamVorti {

/*!
 * \class ParDcpse2d
 * \brief Parallel DC-PSE derivatives in 2D.
 *
 * This class computes DCPSE derivative operators in parallel using ghost
 * node communication. The output is HypreParMatrix objects suitable for
 * parallel linear algebra operations.
 *
 * Matrix structure:
 * - Diagonal block: Local rows × Local columns (owned DOFs)
 * - Off-diagonal block: Local rows × Ghost columns (non-owned DOFs)
 *
 * Algorithm:
 * 1. Call Update() to exchange ghosts and compute stencils
 * 2. For each local node, compute DCPSE coefficients using extended k-NN
 * 3. Populate diagonal/off-diagonal blocks based on neighbor ownership
 * 4. Construct HypreParMatrix from blocks
 */
class ParDcpse2d : public ParDcpse
{
public:
    /*!
     * \brief ParDcpse2d constructor.
     * \param gf Parallel grid function containing nodal coordinates
     * \param num_neighbors Number of neighbors for k-NN search
     */
    ParDcpse2d(mfem::ParGridFunction &gf, int num_neighbors);

    /*!
     * \brief ParDcpse2d destructor.
     */
    virtual ~ParDcpse2d();

    /*!
     * \brief Update derivative matrices.
     *
     * This performs:
     * 1. Ghost node exchange via ParSupportDomain::Update()
     * 2. DCPSE stencil computation for all local nodes
     * 3. HypreParMatrix assembly from diagonal/off-diagonal blocks
     */
    void Update() override;

    /*!
     * \brief Save derivative matrix to file.
     * \param deriv Derivative name ("dx", "dy", "dxx", "dyy", "dxy")
     * \param filename Output file path
     */
    void SaveDerivToFile(const std::string &deriv, const std::string &filename) const override;

    /*!
     * \brief Get derivative matrix by index.
     * \param i Derivative index (0=dx, 1=dy)
     * \return Reference to HypreParMatrix
     */
    const mfem::HypreParMatrix & D(int i) const override {
        switch (i)
        {
        case 0: return *this->sh_func_dx_;
        case 1: return *this->sh_func_dy_;
        }
        MFEM_ABORT( "ParDcpse2d::D: Index " << i << " out of bounds (must be 0 or 1 for 2D)." );
    }

    /*!
     * \brief Get first derivative in x direction.
     */
    inline const mfem::HypreParMatrix & ShapeFunctionDx() const {
        return *this->sh_func_dx_;
    }

    /*!
     * \brief Get first derivative in y direction.
     */
    inline const mfem::HypreParMatrix & ShapeFunctionDy() const {
        return *this->sh_func_dy_;
    }

    /*!
     * \brief Get second derivative in xx direction.
     */
    inline const mfem::HypreParMatrix & ShapeFunctionDxx() const {
        return *this->sh_func_dxx_;
    }

    /*!
     * \brief Get second derivative in yy direction.
     */
    inline const mfem::HypreParMatrix & ShapeFunctionDyy() const {
        return *this->sh_func_dyy_;
    }

    /*!
     * \brief Get mixed second derivative in xy direction.
     */
    inline const mfem::HypreParMatrix & ShapeFunctionDxy() const {
        return *this->sh_func_dxy_;
    }

private:
    /*!
     * \brief Build derivative matrices using extended k-NN with ghosts.
     *
     * This is the core DCPSE algorithm, adapted from serial Dcpse2d::Update().
     * The main difference is assembly into diagonal/off-diagonal blocks.
     */
    void BuildDerivativeMatrices();

    /*!
     * \brief Get column map for HypreParMatrix off-diagonal block.
     *
     * This creates a mapping from ghost global DOF indices to column indices
     * in the off-diagonal block.
     *
     * \return Array of global DOF indices for each ghost node
     */
    HYPRE_BigInt* GetOffDiagonalColumnMap();

    // Derivative matrices (owned by this class)
    mfem::HypreParMatrix* sh_func_dx_;   /*!< 1st derivative wrt x */
    mfem::HypreParMatrix* sh_func_dy_;   /*!< 1st derivative wrt y */
    mfem::HypreParMatrix* sh_func_dxx_;  /*!< 2nd derivative wrt xx */
    mfem::HypreParMatrix* sh_func_dyy_;  /*!< 2nd derivative wrt yy */
    mfem::HypreParMatrix* sh_func_dxy_;  /*!< Mixed 2nd derivative wrt xy */

    // Storage for diagonal and off-diagonal blocks (must outlive HypreParMatrix)
    mfem::SparseMatrix *diag_dx_, *offdiag_dx_;
    mfem::SparseMatrix *diag_dy_, *offdiag_dy_;
    mfem::SparseMatrix *diag_dxx_, *offdiag_dxx_;
    mfem::SparseMatrix *diag_dyy_, *offdiag_dyy_;
    mfem::SparseMatrix *diag_dxy_, *offdiag_dxy_;

    // Column map for off-diagonal blocks
    HYPRE_BigInt* col_map_offd_;

    // Global row/column partitioning
    HYPRE_BigInt* row_starts_;
    HYPRE_BigInt* col_starts_;

    // Input ParFiniteElement Space (for correct matrix partitioning)
    // NOTE: par_fespace_ (inherited from ParSupportDomain) is for COORDINATES (vector field)
    // but derivative matrices must match the INPUT field's partitioning (scalar/vector)
    mfem::ParFiniteElementSpace* input_pfes_;
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_APPROXIMANTS_PAR_DCPSE_2D_HPP_

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
 *
 * Contributors (alphabetically):
 *      George C. BOURANTAS
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

#include "StreamVorti/approximants/par_dcpse_2d.hpp"

#ifdef MFEM_USE_MPI

#include <filesystem>
#include <cfloat>
#include <climits>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace StreamVorti {

ParDcpse2d::ParDcpse2d(mfem::ParGridFunction &gf, int num_neighbors)
    : ParDcpse(gf, num_neighbors),
      sh_func_dx_(nullptr),
      sh_func_dy_(nullptr),
      sh_func_dxx_(nullptr),
      sh_func_dyy_(nullptr),
      sh_func_dxy_(nullptr),
      diag_dx_(nullptr), offdiag_dx_(nullptr),
      diag_dy_(nullptr), offdiag_dy_(nullptr),
      diag_dxx_(nullptr), offdiag_dxx_(nullptr),
      diag_dyy_(nullptr), offdiag_dyy_(nullptr),
      diag_dxy_(nullptr), offdiag_dxy_(nullptr),
      col_map_offd_(nullptr),
      row_starts_(nullptr),
      col_starts_(nullptr),
      input_pfes_(nullptr)
{
    // Get row/column partitioning from INPUT ParGridFunction's FEspace
    // NOT from coordinate FEspace (par_fespace_)!
    //
    // Reason: par_fespace_ is created by ParSupportDomain for coordinate storage
    // and may have different dimension (2D vector vs 1D scalar) than input field.
    // Derivative matrices must match input field's partitioning.
    input_pfes_ = gf.ParFESpace();
    row_starts_ = input_pfes_->GetTrueDofOffsets();
    col_starts_ = input_pfes_->GetTrueDofOffsets();
}

ParDcpse2d::~ParDcpse2d()
{
    // Clean up HypreParMatrix objects
    delete sh_func_dx_;
    delete sh_func_dy_;
    delete sh_func_dxx_;
    delete sh_func_dyy_;
    delete sh_func_dxy_;

    // Clean up diagonal/off-diagonal blocks
    delete diag_dx_; delete offdiag_dx_;
    delete diag_dy_; delete offdiag_dy_;
    delete diag_dxx_; delete offdiag_dxx_;
    delete diag_dyy_; delete offdiag_dyy_;
    delete diag_dxy_; delete offdiag_dxy_;

    // Clean up column map
    delete[] col_map_offd_;
}

void ParDcpse2d::Update()
{
    mfem::StopWatch timer;
    timer.Start();

    if (rank_ == 0) {
        std::cout << "ParDcpse2d::Update: Starting parallel DCPSE computation" << std::endl;
    }

    // Step 1: Exchange ghost nodes and build extended k-NN tree
    // This calls ParSupportDomain::Update() which does:
    // - ComputeSupportRadiuses() using local nodes
    // - ExchangeGhostCoordinates() via MPI
    // - BuildExtendedKnnTree() with local + ghost nodes
    ParSupportDomain::Update();

    // Step 2: Build derivative matrices using extended k-NN
    BuildDerivativeMatrices();

    if (rank_ == 0) {
        std::cout << "ParDcpse2d::Update: Completed in " << timer.RealTime() << " s" << std::endl;
    }
}

void ParDcpse2d::BuildDerivativeMatrices()
{
    mfem::StopWatch timer;
    timer.Start();

    // Get extended neighbor indices (includes ghosts)
    std::vector<std::vector<int> > support_nodes_ids = this->NeighborIndices();

    std::vector<double> support_radiuses = this->SupportRadiuses();

    int nnodes_local = this->NumLocalNodes();
    int num_ghosts = this->NumGhostNodes();

    // Initialize diagonal and off-diagonal blocks
    // Diagonal: nnodes_local × nnodes_local
    // Off-diagonal: nnodes_local × num_ghosts (nullptr if no ghosts)
    diag_dx_ = new mfem::SparseMatrix(nnodes_local, nnodes_local);
    offdiag_dx_ = (num_ghosts > 0) ? new mfem::SparseMatrix(nnodes_local, num_ghosts) : nullptr;

    diag_dy_ = new mfem::SparseMatrix(nnodes_local, nnodes_local);
    offdiag_dy_ = (num_ghosts > 0) ? new mfem::SparseMatrix(nnodes_local, num_ghosts) : nullptr;

    diag_dxx_ = new mfem::SparseMatrix(nnodes_local, nnodes_local);
    offdiag_dxx_ = (num_ghosts > 0) ? new mfem::SparseMatrix(nnodes_local, num_ghosts) : nullptr;

    diag_dyy_ = new mfem::SparseMatrix(nnodes_local, nnodes_local);
    offdiag_dyy_ = (num_ghosts > 0) ? new mfem::SparseMatrix(nnodes_local, num_ghosts) : nullptr;

    diag_dxy_ = new mfem::SparseMatrix(nnodes_local, nnodes_local);
    offdiag_dxy_ = (num_ghosts > 0) ? new mfem::SparseMatrix(nnodes_local, num_ghosts) : nullptr;

    // Track condition number statistics
    int min_supp_nodes = INT_MAX;
    int max_supp_nodes = INT_MIN;
    double min_cond_A1 = DBL_MAX;
    double max_cond_A1 = DBL_MIN;
    bool local_abort = false;

    // Helper lambda to get node coordinates (local or ghost)
    // For local nodes, use tdof_to_ldof_ to convert true DOF index → local DOF (vertex) index,
    // because local_coords_ data is indexed by local vertex, not true DOF.
    auto GetNodeCoord = [&](int idx) -> std::pair<double, double> {
        if (idx < nnodes_local) {
            // Local node: idx is true DOF index, convert to local DOF for coordinate lookup
            mfem::ParFiniteElementSpace* pfes = this->par_fespace_;
            int ldof = this->tdof_to_ldof_[idx];
            double x = (*this->local_coords_)(pfes->DofToVDof(ldof, 0));
            double y = (*this->local_coords_)(pfes->DofToVDof(ldof, 1));
            return {x, y};
        } else {
            // Ghost node
            const double* coords = this->GhostNodeCoord(idx);
            return {coords[0], coords[1]};
        }
    };

    // Iterate over all LOCAL nodes only (we only compute stencils for owned DOFs)
    for (int node_id = 0; node_id < nnodes_local; ++node_id)
    {
        int nsupp = support_nodes_ids[node_id].size();
        min_supp_nodes = std::min(min_supp_nodes, nsupp);
        max_supp_nodes = std::max(max_supp_nodes, nsupp);

        if (nsupp < 1)
        {
            MFEM_WARNING("ParDcpse2d: Node " << node_id << " on rank " << rank_ << " has no neighbors. Abort!");
            local_abort = true;
            goto finish;
        }

        auto [node_X, node_Y] = GetNodeCoord(node_id);

        // Monomial basis - Vandermonde matrices
        Eigen::MatrixXd V1(nsupp, 6);
        Eigen::MatrixXd V2(nsupp, 5);
        Eigen::VectorXd expWd(nsupp);
        std::vector<Eigen::Triplet<double> > expW;
        expW.reserve(nsupp);

        double epsilon = 0.3 * support_radiuses[node_id];

        // Iterate over node's neighbors (includes ghosts)
        for (size_t it = 0; it < nsupp; ++it)
        {
            int neigh_id = support_nodes_ids[node_id][it];
            auto [neigh_X, neigh_Y] = GetNodeCoord(neigh_id);

            double x = (node_X - neigh_X) / epsilon;
            double y = (node_Y - neigh_Y) / epsilon;

            // Set V1 vandermonde matrix for 1st derivatives
            V1.coeffRef(it, 0) = 1.;
            V1.coeffRef(it, 1) = x;
            V1.coeffRef(it, 2) = y;
            V1.coeffRef(it, 3) = x*x;
            V1.coeffRef(it, 4) = x*y;
            V1.coeffRef(it, 5) = y*y;

            double neigh_dist_sq = (neigh_X - node_X) * (neigh_X - node_X) +
                                   (neigh_Y - node_Y) * (neigh_Y - node_Y);

            // Set expWd vector
            double Wd = -(neigh_dist_sq / (epsilon*epsilon));
            expWd(it) = std::exp(Wd);

            // Set expW triplet
            double val = -(neigh_dist_sq / (2.*epsilon*epsilon));
            double temp = std::exp(val);
            expW.emplace_back(Eigen::Triplet<double>(it, it, temp));

            // Set V2 vandermonde matrix for 2nd derivatives
            V2(it, 0) = x;
            V2(it, 1) = y;
            V2(it, 2) = x*x;
            V2(it, 3) = x*y;
            V2(it, 4) = y*y;
        }

        Eigen::SparseMatrix<double> E(nsupp, nsupp);
        E.setFromTriplets(expW.begin(), expW.end());

        // 1st order
        Eigen::Matrix<double, 6, 6> A1;
        Eigen::MatrixXd B1(nsupp, 6);
        B1 = E*V1;
        A1 = B1.transpose()*B1;

        // Condition number check
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A1);
            double condA1 = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            min_cond_A1 = std::min(min_cond_A1, condA1);
            max_cond_A1 = std::max(max_cond_A1, condA1);

            if (condA1 > this->cond_A_limit_abort || condA1 <= 0.0)
            {
                MFEM_WARNING("ParDcpse2d: cond(A1) = " << condA1
                             << " exceeds limit at node " << node_id
                             << " on rank " << rank_ << ". Abort!");
                local_abort = true;
                goto finish;
            }
        }

        // x & y derivatives vectors
        Eigen::Matrix<double, 6, 1> bx, by;
        bx.setZero(); bx(1) = -1.;
        by.setZero(); by(2) = -1.;

        // Solve linear systems
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 6, 6> > sol1(A1);
        Eigen::Matrix<double, 6, 1> aTx = sol1.solve(bx);
        Eigen::Matrix<double, 6, 1> aTy = sol1.solve(by);

        Eigen::VectorXd coeffsx = (V1*aTx).cwiseProduct(expWd);
        Eigen::VectorXd coeffsy = (V1*aTy).cwiseProduct(expWd);

        // 2nd order
        Eigen::Matrix<double, 5, 5> A2;
        Eigen::MatrixXd B2(nsupp, 5);
        B2 = E*V2;
        A2 = B2.transpose()*B2;

        // xx, yy, xy derivatives vectors
        Eigen::Matrix<double, 5, 1> bxx, byy, bxy;
        bxx.setZero(); bxx(2) = 2.;
        byy.setZero(); byy(4) = 2.;
        bxy.setZero(); bxy(3) = 1.;

        // Solve linear systems
        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 5, 5> > sol2(A2);
        Eigen::VectorXd aTxx = sol2.solve(bxx);
        Eigen::VectorXd aTyy = sol2.solve(byy);
        Eigen::VectorXd aTxy = sol2.solve(bxy);

        Eigen::VectorXd coeffsxx = (V2*aTxx).cwiseProduct(expWd);
        Eigen::VectorXd coeffsyy = (V2*aTyy).cwiseProduct(expWd);
        Eigen::VectorXd coeffsxy = (V2*aTxy).cwiseProduct(expWd);

        // Shape functions derivatives components
        Eigen::VectorXd valX(nsupp), valY(nsupp);
        Eigen::VectorXd valXX(nsupp), valYY(nsupp), valXY(nsupp);
        valX.setZero(); valY.setZero();
        valXX.setZero(); valYY.setZero(); valXY.setZero();

        double sumCoeffsx = 0., sumCoeffsy = 0.;
        double sumCoeffsxx = 0., sumCoeffsyy = 0., sumCoeffsxy = 0.;

        for (size_t it = 0; it < nsupp; ++it)
        {
            valX(it) = coeffsx(it) / epsilon;
            sumCoeffsx += coeffsx(it);

            valY(it) = coeffsy(it) / epsilon;
            sumCoeffsy += coeffsy(it);

            valXX(it) = coeffsxx(it) / (epsilon*epsilon);
            sumCoeffsxx += coeffsxx(it);

            valYY(it) = coeffsyy(it) / (epsilon*epsilon);
            sumCoeffsyy += coeffsyy(it);

            valXY(it) = coeffsxy(it) / (epsilon*epsilon);
            sumCoeffsxy += coeffsxy(it);
        }

        // Add sum corrections to diagonal
        for (size_t it = 0; it < nsupp; ++it)
        {
            int neigh_id = support_nodes_ids[node_id][it];
            if (neigh_id == node_id) {
                valX(it) = valX(it) + (sumCoeffsx/epsilon);
                valY(it) = valY(it) + (sumCoeffsy/epsilon);
                valXX(it) = valXX(it) - (sumCoeffsxx/(epsilon*epsilon));
                valYY(it) = valYY(it) - (sumCoeffsyy/(epsilon*epsilon));
                valXY(it) = valXY(it) - (sumCoeffsxy/(epsilon*epsilon));
            }
        }

        // Store in sparse matrices (diagonal or off-diagonal)
        for (size_t it = 0; it < nsupp; ++it)
        {
            int neigh_id = support_nodes_ids[node_id][it];

            if (neigh_id < nnodes_local) {
                // Local neighbor -> diagonal block
                diag_dx_->Add(node_id, neigh_id, valX(it));
                diag_dy_->Add(node_id, neigh_id, valY(it));
                diag_dxx_->Add(node_id, neigh_id, valXX(it));
                diag_dyy_->Add(node_id, neigh_id, valYY(it));
                diag_dxy_->Add(node_id, neigh_id, valXY(it));
            } else if (offdiag_dx_) {
                // Ghost neighbor -> off-diagonal block (only if ghosts exist)
                int ghost_idx = neigh_id - nnodes_local;
                offdiag_dx_->Add(node_id, ghost_idx, valX(it));
                offdiag_dy_->Add(node_id, ghost_idx, valY(it));
                offdiag_dxx_->Add(node_id, ghost_idx, valXX(it));
                offdiag_dyy_->Add(node_id, ghost_idx, valYY(it));
                offdiag_dxy_->Add(node_id, ghost_idx, valXY(it));
            }
        }
    } // End iterate over all local nodes

finish:
    // Check if any rank needs to abort
    bool global_abort = false;
    MPI_Allreduce(&local_abort, &global_abort, 1, MPI_C_BOOL, MPI_LOR, comm_);
    if (global_abort) {
        MFEM_ABORT("ParDcpse2d: Ill-conditioned DCPSE stencils detected across ranks");
    }

    // Finalize sparse matrices
    diag_dx_->Finalize(); if (offdiag_dx_) offdiag_dx_->Finalize();
    diag_dy_->Finalize(); if (offdiag_dy_) offdiag_dy_->Finalize();
    diag_dxx_->Finalize(); if (offdiag_dxx_) offdiag_dxx_->Finalize();
    diag_dyy_->Finalize(); if (offdiag_dyy_) offdiag_dyy_->Finalize();
    diag_dxy_->Finalize(); if (offdiag_dxy_) offdiag_dxy_->Finalize();

    // Get column map for off-diagonal blocks
    col_map_offd_ = GetOffDiagonalColumnMap();

    // Build HypreParMatrix objects
    // Global matrix size: use GlobalTrueVSize(), NOT row_starts_[nranks_].
    // GetTrueDofOffsets() returns a 2-element array [local_start, local_end];
    // indexing beyond index 1 is an out-of-bounds read.
    HYPRE_BigInt global_size = (HYPRE_BigInt)input_pfes_->GlobalTrueVSize();

    // Use different constructor depending on whether ghosts exist
    if (num_ghosts > 0) {
        // Multi-rank case: use constructor with off-diagonal blocks
        sh_func_dx_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                               row_starts_, col_starts_,
                                               diag_dx_, offdiag_dx_, col_map_offd_);
        sh_func_dy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                               row_starts_, col_starts_,
                                               diag_dy_, offdiag_dy_, col_map_offd_);
        sh_func_dxx_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_,
                                                diag_dxx_, offdiag_dxx_, col_map_offd_);
        sh_func_dyy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_,
                                                diag_dyy_, offdiag_dyy_, col_map_offd_);
        sh_func_dxy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_,
                                                diag_dxy_, offdiag_dxy_, col_map_offd_);
    } else {
        // Single-rank or no-ghost case: use the 6-argument block-diagonal constructor.
        // The 8-argument constructor requires non-null offd and cmap; use the
        // simpler version when there are no off-diagonal (ghost) contributions.
        sh_func_dx_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                               row_starts_, col_starts_, diag_dx_);
        sh_func_dy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                               row_starts_, col_starts_, diag_dy_);
        sh_func_dxx_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_, diag_dxx_);
        sh_func_dyy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_, diag_dyy_);
        sh_func_dxy_ = new mfem::HypreParMatrix(comm_, global_size, global_size,
                                                row_starts_, col_starts_, diag_dxy_);
    }

    // Report statistics
    int global_min_supp, global_max_supp;
    double global_min_cond, global_max_cond;
    MPI_Reduce(&min_supp_nodes, &global_min_supp, 1, MPI_INT, MPI_MIN, 0, comm_);
    MPI_Reduce(&max_supp_nodes, &global_max_supp, 1, MPI_INT, MPI_MAX, 0, comm_);
    MPI_Reduce(&min_cond_A1, &global_min_cond, 1, MPI_DOUBLE, MPI_MIN, 0, comm_);
    MPI_Reduce(&max_cond_A1, &global_max_cond, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);

    if (rank_ == 0) {
        std::cout << "  Support nodes per node: [" << global_min_supp << ", "
                  << global_max_supp << "]" << std::endl;
        std::cout << "  Condition number range: [" << global_min_cond << ", "
                  << global_max_cond << "]" << std::endl;
        std::cout << "  Matrix assembly completed in " << timer.RealTime() << " s" << std::endl;
    }
}

HYPRE_BigInt* ParDcpse2d::GetOffDiagonalColumnMap()
{
    int num_ghosts = this->NumGhostNodes();
    if (num_ghosts == 0) {
        return nullptr;
    }

    // ParSupportDomain stores GLOBAL geometry-based node indices.
    // For H1 order-1 elements:
    //   - Scalar field (vdim=1): global_dof = global_node_idx
    //   - Vector field (vdim=d, byNODES): global_dof_comp[i] = global_node_idx + i * total_nodes
    //
    // Since DCPSE operates on scalar fields (one component at a time), we use vdim=1.
    HYPRE_BigInt* col_map = new HYPRE_BigInt[num_ghosts];

    for (int g = 0; g < num_ghosts; ++g) {
        // Get the GLOBAL node index for this ghost
        HYPRE_BigInt global_node_idx = this->GetGhostNodeIndex(this->NumLocalNodes() + g);

        // For scalar H1 order-1 fields, global DOF = global node index
        col_map[g] = global_node_idx;
    }

    return col_map;
}

void ParDcpse2d::SaveDerivToFile(const std::string &deriv, const std::string &filename) const
{
    // TODO: Implement parallel matrix output
    if (rank_ == 0) {
        std::cout << "ParDcpse2d::SaveDerivToFile: Not yet implemented for parallel matrices" << std::endl;
    }
}

} // namespace StreamVorti

#endif // MFEM_USE_MPI

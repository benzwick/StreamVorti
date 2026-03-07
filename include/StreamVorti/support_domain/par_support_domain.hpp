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
 *      Weizheng (Will) LI
 *      Konstantinos A. MOUNTRIS
 *      Benjamin F. ZWICK
 */

#ifndef STREAMVORTI_SUPPORT_DOMAIN_PAR_SUPPORT_DOMAIN_HPP_
#define STREAMVORTI_SUPPORT_DOMAIN_PAR_SUPPORT_DOMAIN_HPP_

#include "StreamVorti/support_domain/support_domain.hpp"

#ifdef MFEM_USE_MPI

#include "mfem.hpp"
#include <mpi.h>
#include <vector>
#include <map>

namespace StreamVorti {

/*!
 * \class ParSupportDomain
 * \brief Parallel support domain with ghost node communication for DCPSE.
 *
 * This class extends SupportDomain to handle ghost nodes from neighboring
 * MPI ranks. Ghost nodes are needed when a local node's support domain
 * extends beyond the local partition boundary.
 *
 * Ghost exchange workflow:
 * 1. Identify boundary nodes that may need ghost data
 * 2. Exchange coordinates via MPI between neighboring ranks
 * 3. Build extended k-NN tree including local + ghost nodes
 * 4. Use extended tree for DCPSE stencil computation
 */
class ParSupportDomain : public SupportDomain
{
public:
    /*!
     * \brief ParSupportDomain constructor for parallel DCPSE.
     * \param gf Parallel grid function containing nodal coordinates
     * \param num_neighbors Number of neighbors for k-NN search
     *
     * The constructor:
     * - Stores local coordinates from ParGridFunction
     * - Initializes MPI communicator from gf's ParMesh
     * - Does NOT perform ghost exchange (call Update() separately)
     */
    ParSupportDomain(const mfem::ParGridFunction &gf, int num_neighbors);

    /*!
     * \brief ParSupportDomain destructor.
     */
    virtual ~ParSupportDomain();

    /*!
     * \brief Update ghost nodes and rebuild extended k-NN tree.
     *
     * This method should be called:
     * - After construction to perform initial ghost exchange
     * - After mesh adaptation or particle movement
     *
     * Steps:
     * 1. Compute support radiuses using local nodes only
     * 2. Identify boundary nodes that may need ghosts
     * 3. Exchange ghost coordinates via MPI
     * 4. Build extended k-NN tree with local + ghost nodes
     */
    void Update();

    /*!
     * \brief Compute support radiuses using local coordinates.
     * \param neighs_num Number of neighbors for k-NN search
     *
     * Override base class method to use parallel coordinate storage.
     */
    void ComputeSupportRadiuses(const std::size_t &neighs_num) override;

    /*!
     * \brief Get neighbor indices using extended k-NN tree (local + ghost).
     *
     * This overrides the serial NeighborIndices() to use the extended
     * point cloud that includes ghost nodes.
     *
     * \return Vector of neighbor indices for each local node.
     *         Indices < num_local are local nodes.
     *         Indices >= num_local are ghost nodes (use GhostNodeCoord).
     */
    const std::vector<std::vector<int> > NeighborIndices();

    /*!
     * \brief Get global DOF index for a LOCAL node (coordinate space).
     * \param local_idx Local index (must be < num_local_nodes_)
     * \return Global DOF index in coordinate space
     * \note For ghost nodes, use GetGhostNodeIndex() instead
     */
    HYPRE_BigInt LocalToGlobal(int local_idx) const;

    /*!
     * \brief Get global geometry-based node index for a ghost node.
     * \param ghost_idx Ghost extended index (>= num_local_nodes_)
     * \return Global node index (FEspace-independent, geometry-based)
     */
    HYPRE_BigInt GetGhostNodeIndex(int ghost_idx) const;

    /*!
     * \brief Get coordinates of a ghost node.
     * \param ghost_idx Ghost index (>= num_local_nodes_)
     * \return Pointer to [x, y, z] coordinates (always 3D, unused dims are 0)
     */
    const double* GhostNodeCoord(int ghost_idx) const;

    /*!
     * \brief Get total number of nodes (local + ghost).
     */
    inline int NumExtendedNodes() const {
        return num_local_nodes_ + ghost_coords_.size();
    }

    /*!
     * \brief Get number of local nodes (excluding ghosts).
     */
    inline int NumLocalNodes() const {
        return num_local_nodes_;
    }

    /*!
     * \brief Get number of ghost nodes.
     */
    inline int NumGhostNodes() const {
        return ghost_node_indices_.size();
    }

    /*!
     * \brief Check if index refers to a ghost node.
     */
    inline bool IsGhost(int idx) const {
        return idx >= num_local_nodes_;
    }

    /*!
     * \brief Get MPI rank.
     */
    inline int Rank() const { return rank_; }

    /*!
     * \brief Get number of MPI ranks.
     */
    inline int NumRanks() const { return nranks_; }

protected:
    // Allow derived classes (ParDcpse) to access coordinate data
    mfem::ParFiniteElementSpace* par_fespace_;/*!< Parallel FE space for coordinates */
    mfem::ParGridFunction* local_coords_;     /*!< Local node coordinates [x, y, z] */
    std::vector<int> tdof_to_ldof_;           /*!< Maps scalar true DOF → scalar local DOF (= local vertex).
                                                   local_coords_ is indexed by local DOF, not true DOF. */

private:
    /*!
     * \brief Exchange ghost coordinates with neighboring ranks.
     *
     * Algorithm:
     * 1. Identify boundary nodes (nodes with support radius extending beyond local partition)
     * 2. For each boundary node, determine which ranks might own needed ghosts
     * 3. Use MPI point-to-point communication to exchange coordinate data
     * 4. Store received ghost nodes in ghost_coords_ and ghost_global_ids_
     */
    void ExchangeGhostCoordinates();

    /*!
     * \brief Identify which local nodes are on partition boundaries.
     *
     * A node is a boundary candidate if its support radius might extend
     * into a neighboring partition.
     *
     * \return Vector of local indices for boundary nodes
     */
    std::vector<int> IdentifyBoundaryNodes();

    /*!
     * \brief Build extended k-NN search tree with local + ghost nodes.
     *
     * This creates a CGAL k-d tree containing both local and ghost nodes,
     * used for computing extended support domains.
     */
    void BuildExtendedKnnTree();

    MPI_Comm comm_;                           /*!< MPI communicator */
    int rank_;                                /*!< MPI rank */
    int nranks_;                              /*!< Number of MPI ranks */

    int dim_;                                 /*!< Spatial dimension (2D or 3D) */
    int num_local_nodes_;                     /*!< Number of local nodes (excluding ghosts) */

    std::vector<HYPRE_BigInt> local_to_global_; /*!< Local to global DOF mapping (coordinate space) */

    std::vector<double> ghost_coords_;        /*!< Ghost node coordinates [x1,y1,z1, x2,y2,z2, ...] */
    std::vector<HYPRE_BigInt> ghost_node_indices_; /*!< GLOBAL node indices for ghost nodes (geometry-based, FEspace-independent) */

    bool ghost_exchange_done_;                /*!< Flag indicating ghost exchange completed */
};

} // namespace StreamVorti

#endif // MFEM_USE_MPI

#endif // STREAMVORTI_SUPPORT_DOMAIN_PAR_SUPPORT_DOMAIN_HPP_

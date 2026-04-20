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

#include "StreamVorti/support_domain/par_support_domain.hpp"

#ifdef MFEM_USE_MPI

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>

namespace StreamVorti {

ParSupportDomain::ParSupportDomain(const mfem::ParGridFunction &gf, int num_neighbors)
    : SupportDomain(gf, *gf.ParFESpace()->GetParMesh()),  // Call base constructor (currently stub)
      comm_(gf.ParFESpace()->GetComm()),
      rank_(0),
      nranks_(1),
      dim_(gf.ParFESpace()->GetParMesh()->Dimension()),
      num_local_nodes_(gf.ParFESpace()->GetTrueVSize()),
      par_fespace_(nullptr),
      local_coords_(nullptr),
      ghost_exchange_done_(false)
{
    // Get MPI rank and size
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &nranks_);

    // Validate dimension
    if (dim_ < 1 || dim_ > 3) {
        MFEM_ABORT("ParSupportDomain: Mesh dimension " << dim_ << " not supported (must be 1D, 2D, or 3D)");
    }

    // Validate H1 finite element space
    const mfem::FiniteElementCollection *fec = gf.ParFESpace()->FEColl();
    if (dynamic_cast<const mfem::H1_FECollection*>(fec) == nullptr) {
        MFEM_ABORT("ParSupportDomain: Grid function FE space is not H1.");
    }

    // Store number of neighbors
    this->num_neighbors_ = num_neighbors;

    // Create parallel finite element space for storing coordinates
    mfem::ParMesh* pmesh = gf.ParFESpace()->GetParMesh();
    par_fespace_ = new mfem::ParFiniteElementSpace(pmesh, fec, dim_, mfem::Ordering::byVDIM);

    // Extract local node coordinates from mesh
    local_coords_ = new mfem::ParGridFunction(par_fespace_);
    pmesh->GetNodes(*local_coords_);

    // Build local to global NODE mapping using INPUT gf's partitioning
    // For H1 order-1, global node index = global DOF index (for scalar fields)
    mfem::ParFiniteElementSpace* input_pfes = gf.ParFESpace();
    HYPRE_BigInt* input_tdof_offsets = input_pfes->GetTrueDofOffsets();
    // NOTE: GetTrueDofOffsets() returns rank-local array starting at current rank's offset
    HYPRE_BigInt local_node_offset = input_tdof_offsets[0];  // First element is this rank's offset

    local_to_global_.resize(num_local_nodes_);
    for (int i = 0; i < num_local_nodes_; ++i) {
        local_to_global_[i] = local_node_offset + i;  // Global node index = global DOF for scalar H1
    }

    // Build inverse mapping: scalar true DOF → scalar local DOF (= local vertex index).
    //
    // local_coords_ is a ParGridFunction whose data is stored in LOCAL DOF order
    // (i.e., indexed by local mesh vertex index, which includes shared vertices from
    // other ranks). In contrast, iteration over owned nodes uses TRUE DOF indices
    // 0..num_local_nodes_-1. On multi-rank runs, shared vertices may appear BEFORE
    // owned vertices in the local mesh ordering, so true DOF i ≠ local vertex i.
    // This mapping corrects the index when accessing local_coords_.
    {
        int nldofs = input_pfes->GetNDofs();  // total local scalar DOFs (owned + shared)
        tdof_to_ldof_.assign(num_local_nodes_, -1);
        for (int ldof = 0; ldof < nldofs; ++ldof) {
            int tdof = input_pfes->GetLocalTDofNumber(ldof);
            if (tdof >= 0) {
                tdof_to_ldof_[tdof] = ldof;
            }
        }
    }

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: Initialized on " << nranks_ << " ranks" << std::endl;
        std::cout << "  Dimension: " << dim_ << "D" << std::endl;
        std::cout << "  Local nodes on rank 0: " << num_local_nodes_ << std::endl;
        std::cout << "  Neighbors per node: " << num_neighbors_ << std::endl;
    }
}


ParSupportDomain::~ParSupportDomain()
{
    delete local_coords_;
    delete par_fespace_;
}


void ParSupportDomain::ComputeSupportRadiuses(const std::size_t &neighs_num)
{
    // Override base class implementation to use local_coords_ instead of support_nodes_
    mfem::StopWatch timer;
    timer.Start();

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: computing support radiuses for " << num_local_nodes_ << " local nodes" << std::endl;
    }

    if (neighs_num < 1) {
        MFEM_ABORT("Number of neighbors = " << neighs_num << " < 1.");
    }

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3d;
    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    // Access the base class protected member
    this->support_radiuses_.clear();
    this->support_radiuses_.reserve(num_local_nodes_);

    // Build point cloud from local coordinates
    std::vector<Point_3d> support_coords;
    support_coords.reserve(num_local_nodes_);

    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];  // local DOF (= local vertex) for true DOF i
        double x = 0.0, y = 0.0, z = 0.0;
        if (dim_ > 0) x = (*local_coords_)(par_fespace_->DofToVDof(ldof, 0));
        if (dim_ > 1) y = (*local_coords_)(par_fespace_->DofToVDof(ldof, 1));
        if (dim_ > 2) z = (*local_coords_)(par_fespace_->DofToVDof(ldof, 2));
        support_coords.emplace_back(Point_3d(x, y, z));
    }

    // Build k-NN tree
    Tree tree(support_coords.begin(), support_coords.end());

    // Compute support radius for each local node
    for (int i = 0; i < num_local_nodes_; ++i) {
        Point_3d query = support_coords[i];

        // Find k nearest neighbors
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;
        distances.reserve(neighs_num);

        // Collect distances to all neighbors
        for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
            distances.emplace_back(std::sqrt(it->second));
        }

        // Sort distances and take the farthest neighbor distance as support radius
        std::sort(distances.begin(), distances.end());
        this->support_radiuses_.emplace_back(distances.back());
    }

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: support radius computation completed in "
                  << timer.RealTime() << " s" << std::endl;
    }
}


void ParSupportDomain::RecomputeSupportRadiusesExtended(const std::size_t &neighs_num)
{
    // Recompute support radii using ALL nodes (local + ghost) so that
    // nodes near partition boundaries get the same k-NN distance as in
    // serial. This ensures partition-independent DCPSE operators.
    mfem::StopWatch timer;
    timer.Start();

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3d;
    typedef CGAL::Search_traits_3<K> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    // Build point cloud from local + ghost coordinates
    std::vector<Point_3d> all_coords;
    int num_ghosts = ghost_node_indices_.size();
    all_coords.reserve(num_local_nodes_ + num_ghosts);

    // Local nodes
    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];
        double x = 0.0, y = 0.0, z = 0.0;
        if (dim_ > 0) x = (*local_coords_)(par_fespace_->DofToVDof(ldof, 0));
        if (dim_ > 1) y = (*local_coords_)(par_fespace_->DofToVDof(ldof, 1));
        if (dim_ > 2) z = (*local_coords_)(par_fespace_->DofToVDof(ldof, 2));
        all_coords.emplace_back(Point_3d(x, y, z));
    }

    // Ghost nodes
    for (int g = 0; g < num_ghosts; ++g) {
        double x = ghost_coords_[g * 3 + 0];
        double y = ghost_coords_[g * 3 + 1];
        double z = ghost_coords_[g * 3 + 2];
        all_coords.emplace_back(Point_3d(x, y, z));
    }

    // Build k-NN tree over all nodes
    Tree tree(all_coords.begin(), all_coords.end());

    // Recompute support radius for each local node
    this->support_radiuses_.clear();
    this->support_radiuses_.reserve(num_local_nodes_);

    for (int i = 0; i < num_local_nodes_; ++i) {
        Point_3d query = all_coords[i];
        Neighbor_search search(tree, query, neighs_num);

        std::vector<double> distances;
        distances.reserve(neighs_num);
        for (auto it = search.begin(); it != search.end(); ++it) {
            distances.emplace_back(std::sqrt(it->second));
        }
        std::sort(distances.begin(), distances.end());
        this->support_radiuses_.emplace_back(distances.back());
    }

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: recomputed support radii with extended nodes in "
                  << timer.RealTime() << " s" << std::endl;
    }
}


void ParSupportDomain::Update()
{
    mfem::StopWatch timer;
    timer.Start();

    if (rank_ == 0) {
        std::cout << "ParSupportDomain::Update: Starting ghost exchange" << std::endl;
    }

    // Step 1: Compute support radiuses using LOCAL nodes only.
    // This gives an upper bound (local-only k-NN distance >= global k-NN
    // distance) used to determine which ghost nodes to fetch.
    ComputeSupportRadiuses(num_neighbors_);

    // Step 2: Exchange ghost node coordinates with neighboring ranks
    ExchangeGhostCoordinates();

    // Step 3: Recompute support radiuses using LOCAL + GHOST nodes.
    // The initial local-only radii overestimate the true support radius
    // for nodes near partition boundaries (their actual nearest neighbors
    // are on other ranks). Recomputing with the full extended node set
    // gives radii identical to serial, ensuring partition-independent results.
    RecomputeSupportRadiusesExtended(num_neighbors_);

    // Step 4: Build extended k-NN tree (local + ghost nodes)
    BuildExtendedKnnTree();

    ghost_exchange_done_ = true;

    if (rank_ == 0) {
        std::cout << "ParSupportDomain::Update: Completed in " << timer.RealTime() << " s" << std::endl;
    }
}


void ParSupportDomain::ExchangeGhostCoordinates()
{
    mfem::StopWatch timer;
    timer.Start();

    ghost_coords_.clear();
    ghost_node_indices_.clear();

    int nprocs;
    MPI_Comm_size(comm_, &nprocs);

    if (nprocs == 1) {
        if (rank_ == 0) {
            std::cout << "  Ghost exchange: single rank, no exchange needed ("
                      << timer.RealTime() << " s)" << std::endl;
        }
        ghost_exchange_done_ = true;
        return;
    }

    // Determine MPI type for HYPRE_BigInt
    MPI_Datatype mpi_hypre_big_int;
    if (sizeof(HYPRE_BigInt) == sizeof(int))            mpi_hypre_big_int = MPI_INT;
    else if (sizeof(HYPRE_BigInt) == sizeof(long))      mpi_hypre_big_int = MPI_LONG;
    else if (sizeof(HYPRE_BigInt) == sizeof(long long))  mpi_hypre_big_int = MPI_LONG_LONG;
    else MFEM_ABORT("ParSupportDomain: Unsupported HYPRE_BigInt size");

    // Prepare local data: coordinates (3D-padded) and global IDs
    std::vector<double> send_coords(num_local_nodes_ * 3);
    std::vector<HYPRE_BigInt> send_global_ids(num_local_nodes_);
    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];
        for (int d = 0; d < dim_; ++d)
            send_coords[i * 3 + d] = (*local_coords_)(par_fespace_->DofToVDof(ldof, d));
        for (int d = dim_; d < 3; ++d)
            send_coords[i * 3 + d] = 0.0;
        send_global_ids[i] = local_to_global_[i];
    }

    // Allgather: collect all node coordinates and global IDs from all ranks.
    // Guarantees every rank has access to every node, producing partition-
    // independent DCPSE operators identical to serial.
    // For large-scale runs (>O(10^5) nodes), replace with a neighbor-only
    // exchange using geometric bounding-box overlap (see issue #45).
    std::vector<int> all_counts(nprocs);
    MPI_Allgather(&num_local_nodes_, 1, MPI_INT,
                  all_counts.data(), 1, MPI_INT, comm_);

    int total_nodes = 0;
    std::vector<int> coord_counts(nprocs), coord_displs(nprocs);
    std::vector<int> id_displs(nprocs);
    for (int r = 0; r < nprocs; ++r) {
        coord_counts[r] = all_counts[r] * 3;
        coord_displs[r] = total_nodes * 3;
        id_displs[r]    = total_nodes;
        total_nodes     += all_counts[r];
    }

    std::vector<double>       all_coords(total_nodes * 3);
    std::vector<HYPRE_BigInt> all_global_ids(total_nodes);

    MPI_Allgatherv(send_coords.data(), num_local_nodes_ * 3, MPI_DOUBLE,
                   all_coords.data(), coord_counts.data(), coord_displs.data(),
                   MPI_DOUBLE, comm_);
    MPI_Allgatherv(send_global_ids.data(), num_local_nodes_, mpi_hypre_big_int,
                   all_global_ids.data(), all_counts.data(), id_displs.data(),
                   mpi_hypre_big_int, comm_);

    // Keep ALL remote nodes as ghosts (no distance filtering).
    // This ensures the extended k-d tree contains the full global point set,
    // so the fuzzy sphere neighbor search returns identical results to serial
    // regardless of partition topology.
    // For large-scale runs, add distance-based filtering with a generous
    // margin beyond the support radius to avoid borderline exclusions.
    std::set<HYPRE_BigInt> local_global_set(local_to_global_.begin(),
                                             local_to_global_.end());
    for (int i = 0; i < total_nodes; ++i) {
        HYPRE_BigInt gid = all_global_ids[i];
        if (local_global_set.count(gid) > 0) continue;  // skip local nodes
        ghost_coords_.push_back(all_coords[i * 3 + 0]);
        ghost_coords_.push_back(all_coords[i * 3 + 1]);
        ghost_coords_.push_back(all_coords[i * 3 + 2]);
        ghost_node_indices_.push_back(gid);
    }

    // Sort ghost arrays by global node index (HYPRE requirement)
    {
        int ng = static_cast<int>(ghost_node_indices_.size());
        std::vector<int> perm(ng);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&](int a, int b) {
            return ghost_node_indices_[a] < ghost_node_indices_[b];
        });
        std::vector<HYPRE_BigInt> si(ng);
        std::vector<double> sc(ng * 3);
        for (int i = 0; i < ng; ++i) {
            si[i]          = ghost_node_indices_[perm[i]];
            sc[i * 3 + 0]  = ghost_coords_[perm[i] * 3 + 0];
            sc[i * 3 + 1]  = ghost_coords_[perm[i] * 3 + 1];
            sc[i * 3 + 2]  = ghost_coords_[perm[i] * 3 + 2];
        }
        ghost_node_indices_ = std::move(si);
        ghost_coords_ = std::move(sc);
    }

    int num_ghosts = static_cast<int>(ghost_node_indices_.size());
    int total_ghosts = 0;
    MPI_Reduce(&num_ghosts, &total_ghosts, 1, MPI_INT, MPI_SUM, 0, comm_);

    if (rank_ == 0) {
        std::cout << "  Ghost exchange (allgather) completed in "
                  << timer.RealTime() << " s" << std::endl;
        std::cout << "  Total ghost nodes across all ranks: " << total_ghosts << std::endl;
        std::cout << "  Ghost nodes on rank 0: " << num_ghosts << std::endl;
    }
}


std::vector<int> ParSupportDomain::IdentifyBoundaryNodes()
{
    // Identify nodes whose support radius extends beyond the local partition bounding box.
    // A node is a boundary candidate if it lies within support_radius of any bounding-box face.
    // This is conservative (may include some interior nodes) but cheap to compute.

    std::vector<int> boundary_nodes;

    // Compute local bounding box
    double local_min[3] = {std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::max()};
    double local_max[3] = {std::numeric_limits<double>::lowest(),
                           std::numeric_limits<double>::lowest(),
                           std::numeric_limits<double>::lowest()};

    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];
        for (int d = 0; d < dim_; ++d) {
            double coord = (*local_coords_)(par_fespace_->DofToVDof(ldof, d));
            local_min[d] = std::min(local_min[d], coord);
            local_max[d] = std::max(local_max[d], coord);
        }
    }

    // For each node, check if support radius extends beyond local bbox
    const std::vector<double>& radiuses = this->SupportRadiuses();
    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];
        double radius = radiuses[i];

        bool near_boundary = false;
        for (int d = 0; d < dim_; ++d) {
            double coord = (*local_coords_)(par_fespace_->DofToVDof(ldof, d));
            if (coord - radius < local_min[d] + 1e-10 ||
                coord + radius > local_max[d] - 1e-10) {
                near_boundary = true;
                break;
            }
        }

        if (near_boundary) {
            boundary_nodes.push_back(i);
        }
    }

    return boundary_nodes;
}


void ParSupportDomain::BuildExtendedKnnTree()
{
    // This method is a placeholder for future optimization.
    // Currently, k-NN search is performed on-demand in NeighborIndices()
    // using both local and ghost coordinates.
    //
    // Future optimization: Pre-build a persistent CGAL k-d tree and store it
    // as a member variable for faster repeated queries.
}


const std::vector<std::vector<int> > ParSupportDomain::NeighborIndices()
{
    mfem::StopWatch timer;
    timer.Start();

    if (!ghost_exchange_done_) {
        MFEM_ABORT("ParSupportDomain::NeighborIndices: Must call Update() before querying neighbors");
    }

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: Computing neighbor indices with extended tree" << std::endl;
    }

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3d;
    typedef std::tuple<Point_3d, int> Point_and_int;
    typedef CGAL::Search_traits_3<K> Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
                                        CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
                                        Traits_base> Traits;
    typedef CGAL::Kd_tree<Traits> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

    // Build list of all points (local + ghost) with their extended indices
    std::list<Point_and_int> points_and_indices;

    // Add local nodes (indices 0 to num_local_nodes_ - 1)
    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];  // local DOF (= local vertex) for true DOF i
        double x = 0.0, y = 0.0, z = 0.0;
        if (dim_ > 0) x = (*local_coords_)(par_fespace_->DofToVDof(ldof, 0));
        if (dim_ > 1) y = (*local_coords_)(par_fespace_->DofToVDof(ldof, 1));
        if (dim_ > 2) z = (*local_coords_)(par_fespace_->DofToVDof(ldof, 2));
        points_and_indices.emplace_back(std::make_tuple(Point_3d(x, y, z), i));
    }

    // Add ghost nodes (indices num_local_nodes_ to num_local_nodes_ + num_ghosts - 1)
    int num_ghosts = ghost_node_indices_.size();
    for (int g = 0; g < num_ghosts; ++g) {
        double x = ghost_coords_[g * 3 + 0];
        double y = ghost_coords_[g * 3 + 1];
        double z = ghost_coords_[g * 3 + 2];
        int extended_idx = num_local_nodes_ + g;
        points_and_indices.emplace_back(std::make_tuple(Point_3d(x, y, z), extended_idx));
    }

    // Build CGAL k-d tree with all points
    Tree tree(points_and_indices.begin(), points_and_indices.end());

    // For each LOCAL node, find neighbors within support radius
    std::vector<std::vector<int> > neighbor_indices(num_local_nodes_);
    std::vector<Point_and_int> domain_nodes;

    for (int i = 0; i < num_local_nodes_; ++i) {
        int ldof = tdof_to_ldof_[i];  // local DOF (= local vertex) for true DOF i
        double x = 0.0, y = 0.0, z = 0.0;
        if (dim_ > 0) x = (*local_coords_)(par_fespace_->DofToVDof(ldof, 0));
        if (dim_ > 1) y = (*local_coords_)(par_fespace_->DofToVDof(ldof, 1));
        if (dim_ > 2) z = (*local_coords_)(par_fespace_->DofToVDof(ldof, 2));
        Point_3d center(x, y, z);

        double radius = this->SupportRadiuses()[i];
        Fuzzy_sphere fs(center, radius);

        // Search for neighbors
        domain_nodes.clear();
        tree.search(std::back_inserter(domain_nodes), fs);

        // Extract indices
        std::vector<int> neigh_indices;
        neigh_indices.reserve(domain_nodes.size());
        for (const auto& d_node : domain_nodes) {
            neigh_indices.push_back(std::get<1>(d_node));
        }

        neighbor_indices[i] = neigh_indices;
    }

    if (rank_ == 0) {
        std::cout << "ParSupportDomain: Neighbor computation completed in "
                  << timer.RealTime() << " s" << std::endl;
    }

    return neighbor_indices;
}


HYPRE_BigInt ParSupportDomain::LocalToGlobal(int local_idx) const
{
    if (local_idx < 0 || local_idx >= num_local_nodes_) {
        MFEM_ABORT("ParSupportDomain::LocalToGlobal: Invalid local index " << local_idx
                   << " (must be 0 to " << num_local_nodes_-1 << ")");
    }

    // This returns coordinate-space DOF ID for LOCAL nodes only
    // For ghost nodes, use GetGhostNodeIndex() to get the geometry-based node index
    return local_to_global_[local_idx];
}


HYPRE_BigInt ParSupportDomain::GetGhostNodeIndex(int ghost_idx) const
{
    if (!IsGhost(ghost_idx)) {
        MFEM_ABORT("ParSupportDomain::GetGhostNodeIndex: Index " << ghost_idx << " is not a ghost node");
    }

    int g = ghost_idx - num_local_nodes_;
    if (g < 0 || g >= (int)ghost_node_indices_.size()) {
        MFEM_ABORT("ParSupportDomain::GetGhostNodeIndex: Ghost index " << g << " out of range");
    }

    return ghost_node_indices_[g];
}


const double* ParSupportDomain::GhostNodeCoord(int ghost_idx) const
{
    if (!IsGhost(ghost_idx)) {
        MFEM_ABORT("ParSupportDomain::GhostNodeCoord: Index " << ghost_idx << " is not a ghost node");
    }

    int g = ghost_idx - num_local_nodes_;
    if (g < 0 || g >= (int)ghost_node_indices_.size()) {
        MFEM_ABORT("ParSupportDomain::GhostNodeCoord: Ghost index " << g << " out of range");
    }

    return &ghost_coords_[g * 3];
}

} // namespace StreamVorti

#endif // MFEM_USE_MPI

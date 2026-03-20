# Parallel DCPSE Implementation Design

StreamVorti's parallel solver (`StreamVorti_par`) distributes the
stream-function/vorticity formulation across MPI ranks using MFEM's parallel
infrastructure and HYPRE's distributed linear algebra. This document describes
the full design: class hierarchy, data flow, key algorithms, and the MFEM
parallel concepts that make it correct.

---

## Table of Contents

1. [Class Hierarchy](#1-class-hierarchy)
2. [MFEM Parallel Concepts](#2-mfem-parallel-concepts)
3. [ParSupportDomain — Ghost Exchange](#3-parsupportdomain--ghost-exchange)
4. [ParDcpse2d — Parallel Derivative Matrices](#4-pardcpse2d--parallel-derivative-matrices)
5. [streamvorti_par.cpp — Time-Stepping Solver](#5-streamvorti_parcpp--time-stepping-solver)
6. [Data Flow: One Timestep](#6-data-flow-one-timestep)
7. [Critical Invariant: True DOF vs Local DOF](#7-critical-invariant-true-dof-vs-local-dof)
8. [Solver Choice](#8-solver-choice)
9. [Post-Processing: Parallel Centerline Output](#9-post-processing-parallel-centerline-output)
10. [Known Limitations and Future Work](#10-known-limitations-and-future-work)

---

## 1. Class Hierarchy

```
SupportDomain          (serial, src/support_domain/)
└── ParSupportDomain   (parallel, src/support_domain/par_support_domain.cpp)
    └── ParDcpse       (parallel base, src/approximants/par_dcpse.cpp)
        └── ParDcpse2d (parallel 2D, src/approximants/par_dcpse_2d.cpp)
```

### Inheritance responsibilities

| Class | Owns |
| --- | --- |
| `SupportDomain` | Serial k-NN tree, support radius computation (stub for parallel) |
| `ParSupportDomain` | MPI communicator, local coordinates, ghost exchange, `tdof_to_ldof_` mapping |
| `ParDcpse` | MPI rank/size mirrored from `ParSupportDomain`, condition-number limits |
| `ParDcpse2d` | Five `HypreParMatrix` derivative operators (dx, dy, dxx, dyy, dxy), diagonal/off-diagonal `SparseMatrix` blocks, `col_map_offd_` |

All parallel classes are compiled only when `MFEM_USE_MPI` is defined.

---

## 2. MFEM Parallel Concepts

Understanding these three index spaces is essential to all parallel code.

### 2.1 Local DOF vs True DOF

For a `ParFiniteElementSpace` (H1, order 1) on a parallel mesh:

| Term | Meaning | Size |
| --- | --- | --- |
| **Local DOF** (`ldof`) | Index into the local mesh vertex array. Includes both owned vertices and shared ("ghost") vertices from neighbouring ranks. Shared vertices appear at **lower** local indices (lower global vertex number first). | `GetNDofs()` |
| **True DOF** (`tdof`) | Index of an **owned** vertex. These are the entries in solution vectors (`vorticity`, `streamfunction`, …). | `GetTrueVSize()` |

The mapping between the two is **not** the identity for `np > 1`:

```cpp
// true DOF → local DOF (built in ParSupportDomain constructor)
int nldofs = pfes->GetNDofs();
tdof_to_ldof_.assign(num_local_nodes_, -1);
for (int ldof = 0; ldof < nldofs; ++ldof) {
    int tdof = pfes->GetLocalTDofNumber(ldof);
    if (tdof >= 0) tdof_to_ldof_[tdof] = ldof;
}
```

**Rule:** whenever you need the physical coordinates of true DOF `t`, use
`mesh->GetVertex(tdof_to_ldof_[t])` or
`(*local_coords_)(pfes->DofToVDof(tdof_to_ldof_[t], d))`. Never use
`mesh->GetVertex(t)` directly.

### 2.2 `ParGridFunction` coordinate storage

`local_coords_` is a `ParGridFunction` on a vector FE space of dimension `dim`.
Its data array is laid out by `byVDIM` ordering: `[x₀, y₀, x₁, y₁, …]` where
the index is the **local DOF** (local vertex), not the true DOF.

`DofToVDof(ldof, d)` converts a local scalar DOF and component index to the
flat data index in the `ParGridFunction` data array.

### 2.3 Global DOF offsets and `HypreParMatrix` partitioning

`GetTrueDofOffsets()` returns a 2-element array `[start, end)` for this rank.
`GlobalTrueVSize()` is the total number of true DOFs across all ranks.

`HypreParMatrix` stores the matrix in HYPRE's distributed CSR format:

```
[ diag block | off-diag block ]
```

- **Diagonal block:** rows and columns owned by this rank (local × local).
- **Off-diagonal block:** rows owned by this rank, columns owned by other
  ranks (local × num_ghosts). Column indices are stored as **local ghost
  indices** (0-based); the mapping to global column indices is in
  `col_map_offd_`, which **must be strictly ascending**.

---

## 3. ParSupportDomain — Ghost Exchange

`ParSupportDomain::Update()` is the entry point. It runs three sequential
steps before any derivative computation:

```
Update()
├── ComputeSupportRadiuses()   — local k-NN, finds max support radius per node
├── ExchangeGhostCoordinates() — neighbour-only MPI, collects remote nodes
└── BuildExtendedKnnTree()     — placeholder; tree is built on-demand in NeighborIndices()
```

### 3.1 ComputeSupportRadiuses

Builds a CGAL k-d tree from local owned nodes (using `tdof_to_ldof_` for
correct coordinate access) and computes the k-th nearest neighbour distance as
the support radius for each node. This radius determines how far the DCPSE
stencil of each node extends and therefore which ghost nodes are needed.

### 3.2 ExchangeGhostCoordinates

**Algorithm:**

1. Call `pmesh->ExchangeFaceNbrData()` (idempotent) to populate
   `face_nbr_group`. Use `pmesh->GetNFaceNeighbors()` and
   `pmesh->GetFaceNbrRank(fn)` to identify the set of face-neighbour ranks.

2. Exchange local node counts with each face neighbour via `MPI_Isend` /
   `MPI_Irecv` (phase A).

3. Exchange all local coordinates and global node indices with each face
   neighbour (phase B). Only face-neighbour ranks communicate; diagonal-only
   neighbours are excluded (acceptable because DCPSE support radii are much
   smaller than partition width in practice).

4. Each rank filters the received nodes: a remote node is kept as a ghost if
   it falls within the support radius of any *boundary node* (a local node
   whose support radius extends to or beyond the local bounding box).

5. Ghost arrays (`ghost_coords_`, `ghost_node_indices_`) are sorted ascending
   by global index so that `col_map_offd_` satisfies HYPRE's requirement.

**Ghost node indexing:** ghost nodes are addressed with extended index
`num_local_nodes_ + g` (0-based `g`). The mapping
`ghost_node_indices_[g]` gives the global true DOF of ghost `g`, which becomes
`col_map_offd_[g]`.

**Communication pattern:** neighbour-only — O(N_face_neighbours × N_local)
instead of O(N_total). Communication is restricted to face-neighbour ranks
using MFEM's `ParMesh::GetNFaceNeighbors()` / `GetFaceNbrRank(fn)` API.

### 3.3 NeighborIndices

Builds a CGAL k-d tree over all extended nodes (local + ghost) and performs a
fuzzy-sphere (radius) query for each local node. Returns
`vector<vector<int>> neighbor_indices` where indices `< num_local` are local
nodes and indices `>= num_local` are ghosts.

---

## 4. ParDcpse2d — Parallel Derivative Matrices

### 4.1 Matrix structure

For each of the five derivative operators (∂/∂x, ∂/∂y, ∂²/∂x², ∂²/∂y²,
∂²/∂x∂y), `ParDcpse2d` assembles a global `HypreParMatrix` from two local
`SparseMatrix` blocks:

```
HypreParMatrix = [ diag (nnodes_local × nnodes_local) |
                   offdiag (nnodes_local × num_ghosts) ]
```

Columns of `diag` are true DOF indices on this rank (0-based local).
Columns of `offdiag` are ghost indices (0-based local ghost), with the actual
global column indices stored in `col_map_offd_`.

All five matrices share the same `col_map_offd_` (ghost columns are identical
for all derivative directions).

### 4.2 DCPSE stencil computation (BuildDerivativeMatrices)

For each local node `i`:

1. Retrieve extended neighbours from `NeighborIndices()`.
2. Look up physical coordinates via `GetNodeCoord(idx)`:
   - If `idx < nnodes_local`: use `tdof_to_ldof_[idx]` for the local DOF,
     then read from `local_coords_`.
   - If `idx >= nnodes_local`: read directly from `ghost_coords_`.
3. Build Vandermonde matrices V1 (6×nsupp, for 1st-order) and V2 (5×nsupp,
   for 2nd-order) in normalised coordinates `(Δx/ε, Δy/ε)`.
4. Build Gaussian weight matrix E = diag(exp(-|xᵢ - xⱼ|²/(2ε²))).
5. Solve the DCPSE normal equations: A = (EV)ᵀ(EV), b = eₖ (unit vector for
   derivative direction k).
6. DCPSE coefficients: `c = A⁻¹ Vᵀ E · exp(-|xᵢ - xⱼ|²/ε²)`.
7. For each coefficient `c[j]`:
   - If neighbour `j` is local → insert into `diag` block at `(i, j)`.
   - If neighbour `j` is a ghost → insert into `offdiag` block at
     `(i, ghost_local_idx)`, where `ghost_local_idx` is the 0-based position
     of this ghost in `ghost_node_indices_`.

### 4.3 HypreParMatrix construction

After all stencils are assembled:

```cpp
col_map_offd_ = GetOffDiagonalColumnMap();  // copies ghost_node_indices_ to HYPRE_BigInt*

if (num_ghosts > 0) {
    sh_func_dx_ = new HypreParMatrix(comm_, global_size, global_size,
                                     row_starts_, col_starts_,
                                     diag_dx_, offdiag_dx_, col_map_offd_);
} else {
    sh_func_dx_ = new HypreParMatrix(comm_, global_size, global_size,
                                     row_starts_, col_starts_, diag_dx_);
}
```

The single-rank constructor is used when there are no ghosts (`np=1`), because
HYPRE requires `offdiag` and `col_map_offd` to be non-null for the 8-argument
constructor.

### 4.4 Condition number check

After building A1 for each node, the ratio `λ_max / λ_min` (from Eigen's
`JacobiSVD`) is checked against `cond_A_limit_abort = 1e30`. If any rank
has an ill-conditioned stencil, `MPI_Allreduce(MPI_LOR)` propagates the abort
flag and all ranks abort cleanly.

---

## 5. streamvorti_par.cpp — Time-Stepping Solver

### 5.1 Initialisation

```
MPI initialise
Parse CLI args
CreateOrLoadMesh → ParMesh (partitioned via METIS or default)
ParFiniteElementSpace fes (H1, order 1, scalar)
ParGridFunction coords_gf → passed to ParDcpse2d
InitialiseParDCPSE → ParDcpse2d::Update() → ghost exchange + matrix assembly
Build Laplacian: laplacian_matrix = ParAdd(&dxx_matrix, &dyy_matrix)
                 EliminateRows(neumann_rows) + SetRow via GetDiag/GetOffd
                 OperatorHandle::EliminateRowsCols(dirichlet_tdof_list)
Create HypreGMRES solver on laplacian_matrix
IdentifyBoundaryNodesPar → boundary/interior node lists (true DOF indices)
```

**Boundary identification (`IdentifyBoundaryNodesPar`):** iterates all local
mesh vertices (`mesh.GetNV()`), calls `pfes.GetLocalTDofNumber(v)` to skip
non-owned vertices, and classifies owned vertices by their physical coordinates
(bottom/top/left/right by ε-tolerance). This correctly handles the case where
shared vertices appear at low local indices and could otherwise be missed.

### 5.2 Time-stepping loop (Euler explicit)

Each timestep:

```
STEP 1 — Compute derivatives (matrix-vector products, parallel):
    ∂ψ/∂y = dy_matrix.Mult(ψ)     → dpsi_dy  (= u)
    ∂ψ/∂x = dx_matrix.Mult(ψ)     → dpsi_dx  (= −v)
    ∂ω/∂x = dx_matrix.Mult(ω)     → domega_dx
    ∂ω/∂y = dy_matrix.Mult(ω)     → domega_dy
    ∂²ω/∂x² = dxx_matrix.Mult(ω)  → d2omega_dx2
    ∂²ω/∂y² = dyy_matrix.Mult(ω)  → d2omega_dy2

STEP 2 — Update vorticity (interior nodes only, local loop):
    ω += dt × (1/Re × (∂²ω/∂x² + ∂²ω/∂y²) − (∂ψ/∂y × ∂ω/∂x − ∂ψ/∂x × ∂ω/∂y))

STEP 3 — Solve Poisson equation (parallel, GMRES):
    −∇²ψ = ω  →  laplacian_matrix × ψ = rhs
    Neumann RHS: rhs[neumann] = 0  (∂ψ/∂n = 0)
    Dirichlet RHS: EliminateBC corrects rhs for boundary values and
                   eliminated column contributions

STEP 4 — Apply velocity boundary conditions (local loop):
    u_velocity = dy_matrix.Mult(ψ)
    v_velocity = −dx_matrix.Mult(ψ)
    Top nodes: u = 1, v = 0  (moving lid)
    Other walls: u = v = 0

STEP 5 — Update boundary vorticity:
    ω_boundary = ∂v/∂x − ∂u/∂y  (from velocity with BCs applied)
```

`HypreParMatrix::Mult` handles all inter-rank communication internally via
HYPRE's `comm_pkg`: for each call the owning rank of each column sends its
solution-vector value to all ranks that have that column in their off-diagonal
block.

### 5.3 Parallel Laplacian and boundary conditions

The Poisson system uses mixed Dirichlet/Neumann BCs. Neumann rows are
replaced first (before Dirichlet column elimination):

```cpp
// Sum DCPSE second-derivative matrices
HypreParMatrix* laplacian_matrix = mfem::ParAdd(&dxx_matrix, &dyy_matrix);
*laplacian_matrix *= -1.0;

// STEP 1: Replace Neumann rows with derivative operators.
// EliminateRows zeros the rows; then SetRow on GetDiag/GetOffd views
// writes the derivative stencil through to the HYPRE internal data.
laplacian_matrix->EliminateRows(neumann_rows);
for (auto& [idx, axis] : neumann_psi_info) {
    // Diagonal block: GetDiag → SetRow writes through
    // Off-diagonal block: GetOffd → SetRow writes through
}

// STEP 2: Eliminate Dirichlet DOFs via OperatorHandle.
// Stores eliminated column entries in Ae_h for per-timestep RHS correction.
mfem::Array<int> ess_tdof_list;
for (int idx : dirichlet_psi_nodes) ess_tdof_list.Append(idx);
mfem::OperatorHandle A_h(laplacian_matrix, false);
mfem::OperatorHandle Ae_h;
Ae_h.EliminateRowsCols(A_h, ess_tdof_list);
```

`dirichlet_psi_nodes` and `neumann_psi_info` contain **true DOF** indices,
consistent with solver vector indexing (`mfem::Vector` of size `TrueVSize`).

Wall ψ constants are computed on the serial mesh before partitioning (using
MFEM boundary element connectivity), then applied to local nodes after
`ParMesh` creation. This ensures correct ψ propagation regardless of which
rank owns the wall-inlet corner vertices.

---

## 6. Data Flow: One Timestep

```
                       All ranks                        Rank 0 only
                      ─────────────                    ─────────────
ψ[TrueVSize]  ──→  HypreParMatrix::Mult  ──→  dpsi_dy, dpsi_dx
ω[TrueVSize]  ──→  HypreParMatrix::Mult  ──→  domega_dx/dy, d2omega_dx2/dy2
                   (inter-rank comms via HYPRE comm_pkg)

domega_*, dpsi_*  ──→  local Euler loop (interior nodes)  ──→  ω_new

rhs = ω_new; rhs[neumann] = 0; EliminateBC(Ae, ess, psi_bc, rhs)
HypreGMRES.Mult(rhs, ψ_new)  ──→  ψ_new
(~80–100 GMRES iterations; global all-reduce per iteration)

ψ_new  ──→  HypreParMatrix::Mult  ──→  u, v
BCs applied locally
u, v  ──→  HypreParMatrix::Mult  ──→  du_dy, dv_dx
ω_boundary = dv_dx − du_dy  (applied to boundary nodes)
```

---

## 7. Critical Invariant: True DOF vs Local DOF

This invariant must be preserved throughout the codebase:

| Operation | Correct index type |
| --- | --- |
| Access `mfem::Vector vorticity[i]` | True DOF `i` |
| Access `(*local_coords_)(DofToVDof(?, d))` | Local DOF `?` = `tdof_to_ldof_[tdof]` |
| Access `mesh->GetVertex(?)` | Local vertex `?` = `tdof_to_ldof_[tdof]` |
| `diag_dx_->Add(row, col, val)` | Both row and col are **true DOF** indices (0-based local) |
| `offdiag_dx_->Add(row, ghost_col, val)` | Row is true DOF, col is ghost index (0-based within this rank's ghost list) |
| `ess_tdof_list.Append(idx)` | True DOF `idx` |

`tdof_to_ldof_` is built once in `ParSupportDomain`'s constructor and stored
as a `protected` member so `ParDcpse2d` can access it directly.

---

## 8. Solver Choice

The DCPSE Laplacian (dxx + dyy) assembled from meshless stencils is **not
symmetric** in general. Symmetry depends on the mesh partition; some partitions
produce a nearly-symmetric matrix (CG converges), others do not (CG stalls).

| Solver | Symmetry required | Result |
| --- | --- | --- |
| `HypreBoomerAMG` | Yes (M-matrix) | Fails at setup; DCPSE matrix is not an M-matrix |
| `HyprePCG` | Yes (SPD) | Stalls for some partitions (e.g. np=2 on ≥40×40 grid) |
| `HypreGMRES` | No | Converges in ~75–100 iterations for all tested cases |

**Default:** `HypreGMRES`. Pass `-solver cg` or `-solver hypre` on the
command line only if you have verified convergence for your specific problem.

---

## 9. Post-Processing: Parallel Centerline Output

`ExtractCenterline` collects velocity values along a line (e.g. x=0.5) and
writes them to a `.dat` file for comparison with benchmark data.

**Parallel-correct implementation:**

1. Each rank iterates its **local DOFs** (`pfes->GetNDofs()`), not true DOFs.
2. For each local DOF `ldof`: check ownership via `GetLocalTDofNumber(ldof)`.
   Skip if `tdof < 0` (shared, owned by another rank).
3. Get physical coordinate from `mesh->GetVertex(ldof)` (local vertex = local
   DOF for P1 elements).
4. Get velocity value from `u_velocity[tdof]` (true DOF).
5. Gather all `(coord, u, v)` tuples to rank 0 via `MPI_Gatherv`.
6. Rank 0 sorts by coordinate and writes the file once.

This avoids the two bugs present in a naïve implementation:
- Mixing local DOF coordinates with true DOF velocity values.
- All ranks writing to the same file simultaneously (data race).

---

## 10. Known Limitations and Future Work

### Ghost exchange: static topology assumed

Ghost exchange runs once during `Update()` and the ghost set is fixed for the
entire simulation. This is correct for a static mesh (lid-driven cavity) but
must be called again after any particle movement or mesh refinement.

### 3D: only 2D derivative matrices implemented

`ParDcpse2d` computes 2D derivatives only. A `ParDcpse3d` subclass following
the same pattern would be needed for 3D problems.

### OpenMP: hybrid MPI+OpenMP available but not performance-tuned

The vorticity update loop in `streamvorti_par.cpp` is guarded with
`#ifdef _OPENMP / #else` to avoid a double-update race. OpenMP is compiled in
via `libomp` (Homebrew) and the `OpenMP::OpenMP_CXX` CMake target. The current
`#pragma omp parallel for` covers only the local Euler update; the CGAL k-NN
search and HYPRE solves are not threaded. A tuning pass (thread-count sweep,
affinity pinning) has not been done.

### Adaptive timestepping with Gershgorin

The `-ats` flag and Gershgorin criterion for adaptive `dt` are implemented but
use `HypreParMatrix::GetDiag()` to access the diagonal of the DCPSE Laplacian.
This is correct but performs an extra HYPRE operation per check; the check
frequency is controlled by `-gsf`.

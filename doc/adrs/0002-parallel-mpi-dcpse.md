# ADR-0002: Parallel MPI DCPSE via MFEM/HYPRE

**Status:** Accepted
**Date:** 2026-03-03
**Authors:** Will Li, Benjamin F. Zwick

---

## Context

StreamVorti solves the 2D stream-function/vorticity formulation of the
Navier-Stokes equations using explicit timestepping and DCPSE (Direct
Collocation Particle Strength Exchange) derivative operators. The serial solver
(`StreamVorti`) builds each derivative as a `SparseMatrix` and solves the
Poisson equation for the stream function using a sequential GMRES solver.

For production-scale problems (fine meshes, high Reynolds numbers, long
integration times), the serial solver is too slow. We need to distribute both:

1. **The DCPSE matrices** — currently single-rank `SparseMatrix` objects.
2. **The Poisson solve** — currently a sequential `GMRESSolver`.

MFEM already provides the parallel infrastructure we need:

- `ParMesh` / `ParFiniteElementSpace` — distributed mesh + DOF management
- `HypreParMatrix` — HYPRE's distributed CSR format, used by all MFEM parallel
  solvers
- `HypreGMRES` / `HyprePCG` — parallel Krylov solvers backed by HYPRE
- MPI face-neighbour communication primitives (`GetFaceNbrRank`, etc.)

The challenge is that DCPSE stencils are meshless: the neighbours of node `i`
are determined by a k-NN radius search, not by mesh connectivity. When the
partition boundary cuts through a node's support domain, the stencil requires
coordinate data from neighbouring MPI ranks — data that MFEM does not
automatically provide for meshless operations.

---

## Decision

### 1. Class hierarchy: inherit from MFEM parallel types

Introduce a parallel class hierarchy that mirrors the serial one:

```
SupportDomain          (serial)
└── ParSupportDomain   (parallel: ghost exchange, extended k-NN)
    └── ParDcpse       (parallel base: HypreParMatrix interface)
        └── ParDcpse2d (parallel 2D: assemble 5 derivative HypreParMatrix)
```

All parallel classes are gated on `MFEM_USE_MPI` and compiled only when MFEM
was built with MPI support.

### 2. Ghost exchange: neighbour-only MPI point-to-point

When a local node's DCPSE stencil extends across a partition boundary, it needs
the coordinates of nodes owned by a neighbouring rank. We call these
**ghost nodes**.

**Approach:** use `ParMesh::ExchangeFaceNbrData()` to identify which MPI ranks
share a mesh face with the local partition, then exchange coordinate data only
with those ranks via `MPI_Isend` / `MPI_Irecv`.

**Why not all-to-all (`MPI_Allgatherv`)?**
An all-to-all broadcast sends O(N_total) data to every rank regardless of
proximity. For meshless methods with local support domains, only face-adjacent
partitions can possibly contribute ghost nodes — diagonal neighbours are
excluded because DCPSE support radii are small relative to partition width.
Neighbour-only communication is O(N_face_neighbours × N_local), which scales
correctly with the number of ranks.

**Two-phase protocol:**

- **Phase A** (size exchange): each rank sends its node count to each face
  neighbour so receive buffers can be sized without a separate metadata message.
- **Phase B** (data exchange): send the full local coordinate array and global
  node ID array to each face neighbour; receive the same from each face
  neighbour.

Received nodes are filtered by distance: a remote node is kept as a ghost only
if it falls within the support radius of at least one boundary-candidate local
node. Ghost arrays are sorted ascending by global index because HYPRE requires
`col_map_offd` to be strictly monotone.

### 3. True DOF vs Local DOF: explicit mapping

MFEM's parallel DOF spaces use two index spaces that are **not** interchangeable
on multi-rank runs:

| Index space | Meaning | Size |
|---|---|---|
| **True DOF** (`tdof`) | Owned vertex on this rank | `GetTrueVSize()` |
| **Local DOF** (`ldof`) | All vertices in the local mesh, including shared vertices from neighbouring ranks | `GetNDofs()` |

Shared vertices appear at lower local indices (sorted by global vertex number).
On `np > 1`, true DOF `i` ≠ local DOF `i`.

We build an explicit `tdof_to_ldof_` mapping once in `ParSupportDomain`'s
constructor (via `GetLocalTDofNumber`) and use it for every coordinate lookup.
This mapping is stored as a `protected` member so `ParDcpse2d` can access it
directly without re-deriving it.

**Rule enforced throughout the codebase:** `mesh->GetVertex(t)` is never
called with a true DOF index. Only `mesh->GetVertex(tdof_to_ldof_[t])` or
`(*local_coords_)(DofToVDof(tdof_to_ldof_[t], d))` is used.

### 4. HypreParMatrix assembly: diagonal + off-diagonal blocks

Each DCPSE derivative operator is assembled into a `HypreParMatrix` via its
two-block structure:

```
HypreParMatrix = [ diag (n_local × n_local) | offdiag (n_local × n_ghost) ]
```

- **Diagonal block:** stencil coefficients where the neighbour is a locally
  owned node. Stored as `SparseMatrix` with true-DOF row and column indices.
- **Off-diagonal block:** stencil coefficients where the neighbour is a ghost
  node. Stored as `SparseMatrix` with true-DOF row and ghost-local column
  indices. The `col_map_offd_` array maps ghost-local → global DOF.

All five derivative operators (∂x, ∂y, ∂²x, ∂²y, ∂²xy) share the same
`col_map_offd_` because the ghost set is determined by geometry, not by the
derivative direction.

### 5. Parallel Poisson solve: HypreGMRES

The DCPSE Laplacian (`dxx + dyy`) is assembled via `mfem::ParAdd` and boundary
conditions are applied in two steps:

1. **Neumann rows** (pressure/outflow boundaries) are replaced with the normal
   derivative operator via `GetDiag`/`GetOffd` CSR block access.
2. **Dirichlet rows/columns** are eliminated via `OperatorHandle::EliminateRowsCols`,
   which stores the eliminated column entries for per-timestep RHS correction.

The resulting distributed matrix is solved with `HypreGMRES`.

**Why GMRES, not PCG or AMG?**

The DCPSE Laplacian is not symmetric in general. Symmetry depends on the
specific mesh partition; some partitions produce a nearly-symmetric matrix
(PCG converges in ~100 iterations), but others do not (PCG stalls at residual
~1e-4 with 500-iteration cap). BoomerAMG fails at setup because the matrix is
not an M-matrix. GMRES converges in 75–100 iterations for all tested rank
counts and grid sizes.

### 6. OpenMP for local loops

The vorticity update loop (interior nodes only, no inter-rank communication)
is parallelised with `#pragma omp parallel for` inside a `#ifdef _OPENMP`
guard. The guard is necessary to avoid a latent double-update bug: without it,
the loop body would execute twice when `_OPENMP` is defined but the pragma is
absent.

OpenMP support on macOS requires `libomp` from Homebrew. The CMakeLists.txt
explicitly searches the Homebrew opt paths (`/usr/local/opt/libomp` for x86,
`/opt/homebrew/opt/libomp` for ARM) before the default system paths.

### 7. Parallel centerline output

`ExtractCenterline` collects the u-velocity profile along a vertical or
horizontal line for Ghia et al. benchmark comparison. The parallel-correct
implementation:

1. Each rank collects its local owned-DOF contributions (using
   `GetLocalTDofNumber` to skip shared DOFs).
2. Data is gathered to rank 0 via `MPI_Gatherv`.
3. Rank 0 sorts by coordinate and writes the file once.

This avoids the two bugs in a naïve implementation: mixing local-DOF
coordinates with true-DOF velocity values, and all ranks writing to the same
file simultaneously.

---

## Consequences

### Accepted

- Serial solver (`StreamVorti`) is completely unchanged. The parallel solver
  (`StreamVorti_par`) is a separate executable that shares the library but
  uses distinct parallel classes.
- `HypreParMatrix::Mult` handles all inter-rank communication automatically
  via HYPRE's `comm_pkg`; application code does not need explicit ghost syncs
  during time-stepping.
- Validated against the Ghia et al. (1982) Re=100 lid-driven cavity benchmark:
  L2 error 0.30% (np=1) and 0.39% (np=4), both within the 5% threshold.

### Known limitations

- **Diagonal neighbours:** the face-neighbour exchange excludes ranks that
  share only a corner (diagonal). For standard DCPSE stencil sizes (30
  neighbours on a 40×40 grid, support radius ~5 grid spacings), partition
  widths are much larger than the support radius and diagonal exchange is not
  needed. Grids coarser than ~10×10 nodes per partition with very large stencil
  counts may require a 2-hop exchange.
- **Static ghost set:** ghost exchange runs once during `Update()`. The ghost
  set is fixed for the entire simulation. Particle movement or dynamic mesh
  refinement would require calling `Update()` again.
- **2D only:** `ParDcpse2d` implements the five 2D derivative operators.
  A `ParDcpse3d` subclass would follow the same pattern.
- **HypreGMRES is unpreconditioned:** no preconditioner is applied. For larger
  problems or higher Reynolds numbers, an ILU or AMG preconditioner matched
  to the DCPSE stencil structure may be needed.

### Performance

On a 40×40 grid (Re=100, dt=0.001, 5000 steps):

| Ranks | ms/iter | Speedup |
|-------|---------|---------|
| np=1  | 5.2     | 1.00×   |
| np=2  | 3.5     | 1.49×   |
| np=4  | 3.3     | 1.58×   |

Scaling is communication-bound at small grid sizes. The 80×80 grid achieves
2.98× speedup at np=4.

---

## Alternatives considered

### All-to-all ghost exchange (`MPI_Allgatherv`)

**Implemented first as Phase 1.** Simple to implement: every rank broadcasts
its full coordinate array. Correct, but O(N_total) data per rank; does not
scale with the number of ranks. Replaced by neighbour-only exchange (this ADR).

### MFEM `GroupCommunicator` for shared-vertex data

MFEM provides a `GroupCommunicator` that communicates only across shared mesh
entities (vertices, edges, faces). It operates on the `ldof` space and transfers
only the shared-vertex values. This is more efficient than our approach for
the data that it covers, but it only communicates shared vertex coordinates —
not the coordinates of all local nodes that might fall within a remote node's
support domain. Interior nodes of one partition can be within the support radius
of a boundary node on the adjacent partition, so shared-vertex communication
alone is not sufficient.

### MFEM `ParGridFunction::ExchangeFaceNbrData`

This is MFEM's built-in mechanism for exchanging grid-function values across
face-neighbour boundaries (used in DG methods). It is the correct approach for
**field values** but operates on the `ldof` space and cannot easily be adapted
to carry our custom ghost coordinate arrays and global node ID metadata.

---

## References

- Bourantas et al. (2016), "Using DC PSE operator discretization in Eulerian
  meshless collocation methods improves their robustness in complex geometries"
- Ghia, Ghia, Shin (1982), "High-Re solutions for incompressible flow using the
  Navier-Stokes equations and a multigrid method"
- MFEM parallel examples: `examples/ex1p.cpp`, `miniapps/meshing/`
- HYPRE Reference Manual (v2.30): `IJ_mv/IJMatrix.h`, `parcsr_mv/`
- `docs/PARALLEL_DESIGN.md` — detailed implementation notes for this feature
- `docs/PARALLEL_RESULTS.md` — full validation and scaling data

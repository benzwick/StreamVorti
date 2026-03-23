# Boundary Conditions for Stream Function-Vorticity Formulation

This document describes how boundary conditions are implemented in
StreamVorti's stream function-vorticity (ψ-ω) solver.

## Governing Equations

The 2D incompressible Navier-Stokes equations in stream
function-vorticity form:

1. **Vorticity transport** (explicit time stepping):

   ∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = (1/Re)·(∂²ω/∂x² + ∂²ω/∂y²)

2. **Stream function Poisson equation** (implicit, solved each timestep):

   ∇²ψ = -ω

3. **Velocity recovery** from stream function:

   u = ∂ψ/∂y
   v = -∂ψ/∂x

4. **Boundary vorticity** (derived from velocity BCs):

   ω|_boundary = ∂v/∂x - ∂u/∂y

## Time Loop Structure

Each timestep follows this sequence:

```
1. Update vorticity (explicit Euler):
     ω^(n+1) = ω^(n) + dt·[(1/Re)·∇²ω - u·∂ω/∂x - v·∂ω/∂y]

2. Apply vorticity BCs at boundaries:
     ω|_b = ∂v/∂x - ∂u/∂y  (computed from velocity field)

3. Solve stream function Poisson equation:
     -∇²ψ = ω  (negated for positive definiteness)
     with Dirichlet/Neumann BCs on ψ

4. Recover velocity from stream function:
     u = ∂ψ/∂y,  v = -∂ψ/∂x

5. Apply velocity BCs at boundaries:
     Overwrite u,v at no-slip/velocity nodes
     Leave u,v at pressure/outflow/slip nodes
```

## Boundary Condition Types

| BC Type    | ψ Condition                | Velocity Step                                |
|------------|----------------------------|----------------------------------------------|
| `no-slip`  | Dirichlet ψ = constant     | u = 0, v = 0                                 |
| `velocity` | Dirichlet ψ = ∫u dy        | u = u_prescribed, v = v_prescribed            |
| `pressure` | Neumann ∂ψ/∂n = 0          | Leave computed values (natural outflow)        |
| `outflow`  | Neumann ∂ψ/∂n = 0          | Leave computed values (natural outflow)        |
| `slip`     | Dirichlet ψ = constant     | Leave computed values (zero normal velocity)   |

### Stream Function at Velocity Boundaries (Inlet)

For a velocity BC on a vertical boundary at x = x₀ with prescribed u(y):

    u = ∂ψ/∂y  →  ψ(x₀, y) = ∫₀ʸ u(x₀, η) dη

For a velocity BC on a horizontal boundary at y = y₀ with prescribed v(x):

    v = -∂ψ/∂x  →  ψ(x, y₀) = ψ_ref - ∫₀ˣ v(ξ, y₀) dξ

The integration is performed numerically using composite Simpson's
rule with 20 sub-intervals per node. Each node computes its ψ value
independently — there are no inter-node dependencies, so the same
code works in both serial and parallel without any communication or
node sorting.

For example, for a vertical inlet at x = x₀ with node at y:

    ψ(y) ≈ Σₖ (h/6) · [u(x₀, aₖ) + 4·u(x₀, mₖ) + u(x₀, bₖ)]

where h = y/20, aₖ = k·h, bₖ = (k+1)·h, mₖ = (aₖ+bₖ)/2.

This approach generalizes to any velocity function defined in SDL,
not just analytically integrable profiles. The von Karman reference
code uses the analytical integral directly, but we support arbitrary
Lisp-defined velocity functions.

### Stream Function on No-Slip Walls

ψ = constant on any impermeable wall (from impermeability u·n = 0,
not from no-slip). Each connected wall segment has its own constant,
determined by mass conservation:

    ψ(B) - ψ(A) = volumetric flow rate between points A and B

The wall constant equals ψ at the corner vertex where the wall meets
a velocity boundary. The solver discovers this boundary topology
using MFEM's boundary element connectivity (`GetBdrElementVertices`):

1. Build a map from each vertex to the set of boundary attributes
   meeting at that vertex
2. Corner vertices have 2+ attributes
3. For each wall BC, find corners shared with velocity BCs
4. Compute ψ at the corner using the velocity function (Simpson's rule)
5. Set all wall nodes to that constant

For the lid-driven cavity (closed domain), all walls share the same
constant ψ = 0, which is the trivial case (no velocity boundaries
to propagate from).

For channel flow:
- Bottom wall at y=0 → corner with inlet at (x₀, 0) → ψ = ∫₀⁰ u dη = 0
- Top wall at y=H → corner with inlet at (x₀, H) → ψ = ∫₀ᴴ u dη = Q

In parallel, the corner detection is performed on the serial mesh
before partitioning, since corner vertices may not be on every rank.

### Neumann BCs for Pressure/Outflow

For pressure or outflow boundaries, the stream function satisfies
∂ψ/∂n = 0 (zero-gradient condition). This is implemented by replacing
the Laplacian matrix row at Neumann nodes with the normal derivative
operator:

- At x = x_outlet: row ← ∂/∂x operator (from DCPSE dx matrix)
- At y = y_outlet: row ← ∂/∂y operator (from DCPSE dy matrix)

The normal direction is determined from `predicate_axis` on the
`BoundaryCondition` struct.

The RHS for Neumann rows is set to 0 at each timestep.

## Implementation Details

### Laplacian Matrix Assembly

The Poisson system is:

    -∇²ψ = ω  (negated for positive definiteness with CG solver)

The matrix is assembled once and factorized before the time loop:

```cpp
// Laplacian = -(dxx + dyy)
laplacian_matrix = new SparseMatrix(dxx_matrix);
laplacian_matrix->Add(1.0, dyy_matrix);
*laplacian_matrix *= -1.0;

// STEP 1: Neumann rows first — replace with normal derivative operator
for (auto& [idx, axis] : neumann_psi_info) {
    laplacian_matrix->EliminateRow(idx);
    // Copy derivative operator row
    deriv.GetRow(idx, cols, vals);
    for (int j = 0; j < cols.Size(); ++j)
        laplacian_matrix->Set(idx, cols[j], vals(j));
}

// STEP 2: Dirichlet rows — replace with identity
for (int idx : dirichlet_psi_nodes) {
    laplacian_matrix->EliminateRow(idx);
    laplacian_matrix->Set(idx, idx, 1.0);
}
```

Neumann rows are processed first so that Dirichlet wins at shared corner
nodes (where two boundary attributes share a vertex). The RHS is set
in the same order: Neumann nodes to 0 first, then Dirichlet nodes to
ψ_prescribed.

The serial solver uses row-only elimination (no column elimination)
because DCPSE matrices have non-symmetric sparsity and
`SparseMatrix::EliminateRowsCols` assumes symmetric sparsity. The
parallel solver uses `OperatorHandle::EliminateRowsCols` with
`HypreParMatrix`, which handles non-symmetric matrices correctly
(see [Parallel Implementation Notes](#parallel-implementation-notes) below).

### Corner Node Handling

MFEM assigns boundary attributes to boundary *elements* (edges in 2D),
not vertices. At corners, two edges with different attributes share a
vertex, so corner nodes appear in multiple attribute lists.

Both the matrix and the RHS are processed Neumann-first, Dirichlet-second,
so Dirichlet always wins at shared corners. This ensures wall-outlet corner
nodes get the wall's ψ value rather than the Neumann condition.

For the cavity, all corners are no-slip regardless of order. For
channels, corner nodes between outlet and walls get the wall's ψ value
(propagated from the inlet integration endpoint).

### Parallel Implementation Notes

1. **Neumann row replacement** — `HypreParMatrix` does not have a direct
   `SetRow` method, but its internal CSR blocks are accessible via
   `GetDiag(SparseMatrix&)` and `GetOffd(SparseMatrix&, HYPRE_BigInt*&)`.
   These return non-owning views that write through to the underlying
   HYPRE data. Neumann rows are first zeroed with `EliminateRows()`,
   then replaced with derivative operator entries via `SetRow` on the
   diagonal and off-diagonal block views. All DCPSE matrices share
   identical sparsity, so `SetRow` always matches `GetRow` column count.

2. **Dirichlet elimination** — after Neumann row replacement, Dirichlet
   DOFs are eliminated via `OperatorHandle::EliminateRowsCols`, which
   calls `HypreParMatrix::EliminateRowsCols` internally. This works
   correctly for non-symmetric DCPSE matrices (unlike the serial
   `SparseMatrix` version, which assumes symmetric sparsity). The
   eliminated column entries are stored in a separate matrix and used
   each timestep by `EliminateBC` to correct the RHS.

3. **Wall ψ constants** — computed on the serial mesh before partitioning
   using MFEM boundary element connectivity. This avoids the problem
   of corner vertices being on different ranks.

4. **Per-node ψ integration** — each rank computes ψ at its own inlet
   nodes independently using Simpson's rule. No MPI communication needed.

### Simply Connected Domains Only

Multiply connected domains (e.g., channel with internal cylinder
obstacle) require additional pressure single-valuedness equations per
internal boundary. The ψ value on an internal boundary is an unknown
determined by the constraint that pressure must be single-valued
around the obstacle. This is not yet implemented. See issue #25.

## SDL Examples

### Lid-Driven Cavity (Closed Domain)

```lisp
(bc lid    :velocity :u 1 :v 0)
(bc bottom :no-slip)
(bc left   :no-slip)
(bc right  :no-slip)
```

All boundaries: Dirichlet ψ = 0.

### Poiseuille Channel (Open Domain)

```lisp
(bc inlet  :velocity :u (fn (x y) (* 4 y (- 1 y))) :v 0)
(bc outlet :pressure 0)
(bc top    :no-slip)
(bc bottom :no-slip)
```

- Inlet: Dirichlet ψ(y) = 2y² - 4y³/3
- Outlet: Neumann ∂ψ/∂x = 0
- Bottom: Dirichlet ψ = 0
- Top: Dirichlet ψ = 2/3 (flow rate Q)

## References

- Von Karman cylinder reference code:
  `von_karman/von_karman/cylinder/src/ns_solver.cpp`
- DCPSE 2D reference (MATLAB):
  `2020-12-dcpse/dcpse2d-1-julia-bz/Lid_driven_streamf_vorticity_DCPSE.m`
- Ghia, Ghia & Shin (1982) — lid-driven cavity benchmark

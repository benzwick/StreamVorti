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

No-slip walls have ψ = constant along the wall. The constant is
determined by continuity with adjacent velocity boundaries at
shared corner nodes:

- Bottom wall → ψ = 0 (reference, from inlet node at y=0)
- Top wall → ψ = Q (total flow rate, from inlet node at y=H)

where Q = ∫₀ᴴ u_inlet(y) dy is the total inlet flow rate.

For the lid-driven cavity (closed domain), all walls share the same
constant ψ = 0, which is the trivial case.

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

// Dirichlet rows: replace with identity
for (int idx : dirichlet_psi_nodes) {
    laplacian_matrix->EliminateRow(idx);
    laplacian_matrix->Set(idx, idx, 1.0);
}

// Neumann rows: replace with normal derivative operator
for (auto& [idx, axis] : neumann_psi_info) {
    laplacian_matrix->EliminateRow(idx);
    // Copy derivative operator row
    deriv.GetRow(idx, cols, vals);
    for (int j = 0; j < cols.Size(); ++j)
        laplacian_matrix->Set(idx, cols[j], vals(j));
}
```

### Corner Node Handling

MFEM assigns boundary attributes to boundary *elements* (edges in 2D),
not vertices. At corners, two edges with different attributes share a
vertex, so corner nodes appear in multiple attribute lists. The last
BC processed wins at these shared nodes — this is standard FEM practice.

For the cavity, all corners are no-slip regardless of order. For
channels, corner nodes between inlet and walls get the wall's ψ value
(propagated from the inlet integration endpoint).

### Parallel Limitations

The parallel solver (`streamvorti_par.cpp`) uses `HypreParMatrix` for
the Laplacian, which does not support arbitrary row replacement.
Neumann BCs are not yet implemented in parallel — pressure/outflow
outlets fall back to Dirichlet ψ = 0. Non-zero Dirichlet ψ at
inlet/walls is fully supported in parallel.

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

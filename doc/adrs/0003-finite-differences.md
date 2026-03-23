# ADR-0003: Finite Differences as Spatial Discretization

**Status:** Proposed
**Date:** 2026-03-18
**Authors:** Benjamin Zwick

## Context

StreamVorti currently uses DC PSE (Direct Collocation Particle Strength Exchange)
as its spatial discretization method. DCPSE is a meshless method that works on
scattered node distributions and computes derivative operators via weighted least
squares on local monomial bases. While powerful and flexible, DCPSE has overhead
from k-NN neighbor searches (via CGAL) and dense local system solves (via Eigen)
that is unnecessary when the problem domain is a structured Cartesian grid.

For structured grids, classical finite differences (FD) provide:
- Simple, well-understood stencils with known truncation error
- No need for neighbor search (grid topology is implicit)
- No need for dense linear algebra (stencil weights are analytical)
- Fast matrix assembly
- Easy verification against analytical solutions

Many benchmark problems in computational fluid dynamics (e.g., lid-driven cavity)
use structured grids, making FD a natural and efficient choice for these cases.

## Decision

We add a finite difference module (`FiniteDiff`, `FiniteDiff2d`, `FiniteDiff3d`)
that provides the same sparse derivative matrix interface as DCPSE but constructs
matrices using standard FD stencils on structured Cartesian MFEM meshes.

### Design

**Class hierarchy:**
```
FiniteDiff (base class)
├── FiniteDiff2d
└── FiniteDiff3d
```

The classes do **not** inherit from `SupportDomain` or `Dcpse` because FD does
not require CGAL neighbor search or Eigen dense solves. Instead, `FiniteDiff`
provides the same user-facing interface:

- `Update()` — builds the derivative matrices
- `D(int i)` — returns the i-th first derivative matrix
- `ShapeFunctionDx()`, `ShapeFunctionDy()`, etc. — named accessors
- `SaveDerivToFile()` — exports matrices in MATLAB format

**Stencils used (2nd-order, default):**

| Derivative | Interior (central) | Boundary (one-sided) | Order |
|-----------|-------------------|---------------------|-------|
| d/dx | (-1, 0, 1) / 2h | (-3, 4, -1) / 2h | O(h²) |
| d²/dx² | (1, -2, 1) / h² | (2, -5, 4, -1) / h² | O(h²) |
| d²/dxdy | (f₊₊ - f₊₋ - f₋₊ + f₋₋) / 4hxhy | one-sided fallback | O(h²) |

**Stencils used (4th-order, optional):**

| Derivative | Interior (central) | Boundary (one-sided) | Order |
|-----------|-------------------|---------------------|-------|
| d/dx | (1, -8, 0, 8, -1) / 12h | (-25, 48, -36, 16, -3) / 12h | O(h⁴) |
| d²/dx² | (-1, 16, -30, 16, -1) / 12h² | (45, -154, 214, -156, 61, -10) / 12h² | O(h⁴) |

Near-boundary (1 node from edge) stencils are also 4th-order accurate.

**Grid validation:**

The `Update()` method performs rigorous validation before assembling stencils:

1. **Node coordinate extraction** from GridFunction (not mesh vertices),
   supporting higher-order FE spaces for visualization.
2. **Structured grid detection** via unique coordinate analysis with
   mesh-aware tolerances scaled to the coordinate range.
3. **Uniform spacing verification** checking that max(|h_i - h_mean|)/h_mean
   < 1e-6 for all intervals along each axis.
4. **Minimum grid size check** ensuring sufficient nodes for the stencil
   order (4 for O(h²), 6 for O(h⁴)).
5. **Complete coverage check** verifying every grid position has exactly
   one DOF (no missing or duplicate nodes).

**Parallel support:**

`ParFiniteDiff2d` constructs `HypreParMatrix` derivative operators for use
with MFEM parallel solvers. Each rank detects the global grid structure via
MPI collective operations, then assembles its local rows. This follows the
same pattern as `ParDcpse2d` but without ghost exchange (FD stencils
reference global DOF indices directly).

### Usage

```cpp
// Create structured mesh and FE space
mfem::Mesh mesh(mfem::Mesh::MakeCartesian2D(40, 40,
    mfem::Element::QUADRILATERAL, false, 1.0, 1.0, false));
mfem::H1_FECollection fec(1, 2);
mfem::FiniteElementSpace fes(&mesh, &fec, 1);
mfem::GridFunction gf(&fes);

// Build FD derivative matrices
StreamVorti::FiniteDiff2d fd(gf);
fd.Update();

// Use exactly like DCPSE
const mfem::SparseMatrix &dx = fd.ShapeFunctionDx();
const mfem::SparseMatrix &dxx = fd.ShapeFunctionDxx();
mfem::Vector dudx(fes.GetNDofs());
dx.Mult(gf, dudx);
```

### Integration with SDL

In the simulation definition language (ADR-0001), FD is selected via the
`:fdm` method keyword:

```lisp
(method :fdm)
```

This is distinct from `:dcpse` (meshless) and `:fem` (finite elements).

## Consequences

### Advantages

1. **Performance**: No CGAL or Eigen dependencies for the FD module. Matrix
   assembly is O(N) with small constant.

2. **Simplicity**: FD stencils are analytically known, making the code easy to
   verify and understand.

3. **Compatibility**: Same `SparseMatrix` interface as DCPSE, so existing
   time-stepping and solver code works without modification.

4. **Verification**: Exact for polynomials up to the stencil order, making
   convergence testing straightforward.

5. **No new dependencies**: Uses only MFEM (already required).

### Disadvantages

1. **Structured grids only**: Unlike DCPSE, FD requires a Cartesian grid.
   Unstructured meshes will fail with an error.

2. **Uniform spacing**: Current implementation assumes uniform grid spacing
   in each direction. Non-uniform grids would need modified stencils.

3. **No shared base class with DCPSE**: `FiniteDiff` and `Dcpse` have similar
   interfaces but no common abstract base. This is intentional — FD does not
   need `SupportDomain` and forcing a common base would add artificial coupling.
   A future refactoring could introduce a pure interface class if needed.

### Future Work

- Non-uniform grid spacing support (stretching functions)
- Compact (Padé) schemes for higher accuracy at same stencil width
- 4th-order stencils for 3D (currently only 2nd-order in 3D)
- `ParFiniteDiff3d` for parallel 3D simulations

## References

- LeVeque, R.J. (2007). *Finite Difference Methods for Ordinary and Partial
  Differential Equations*. SIAM.
- Strikwerda, J.C. (2004). *Finite Difference Schemes and Partial Differential
  Equations*. 2nd ed. SIAM.

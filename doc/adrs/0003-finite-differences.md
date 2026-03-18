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

**Stencils used:**

| Derivative | Interior (central) | Boundary (one-sided) | Order |
|-----------|-------------------|---------------------|-------|
| d/dx | (-1, 0, 1) / 2h | (-3, 4, -1) / 2h | O(h²) |
| d²/dx² | (1, -2, 1) / h² | (2, -5, 4, -1) / h² | O(h²) |
| d²/dxdy | (f₊₊ - f₊₋ - f₋₊ + f₋₋) / 4hxhy | one-sided fallback | O(h²) |

All stencils are second-order accurate.

**Grid detection:**

The `Update()` method automatically detects the grid structure from the MFEM
mesh by extracting unique x, y (and z) coordinates, verifying the node count
matches a structured grid, and building an index map from grid indices to DOF
indices. This works regardless of MFEM's internal DOF ordering.

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

- Non-uniform grid spacing support
- Higher-order stencils (4th-order compact schemes)
- Parallel (MPI) FD via `ParFiniteDiff2d` following the `ParDcpse` pattern
- Integration into the SDL method selection

## References

- LeVeque, R.J. (2007). *Finite Difference Methods for Ordinary and Partial
  Differential Equations*. SIAM.
- Strikwerda, J.C. (2004). *Finite Difference Schemes and Partial Differential
  Equations*. 2nd ed. SIAM.

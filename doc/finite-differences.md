# Finite Differences Module

## Overview

The finite differences (FD) module provides spatial derivative operators for
structured Cartesian grids. It builds sparse matrices that approximate first
and second partial derivatives using standard finite difference stencils.

The FD module is an alternative to DCPSE for problems on structured grids,
offering simpler implementation and faster matrix assembly.

## Classes

### `FiniteDiff` (base class)

Abstract base class in `include/StreamVorti/finite_differences/fd.hpp`.

**Constructor:**
```cpp
FiniteDiff(mfem::GridFunction &gf);
```

**Pure virtual methods:**
- `void Update()` — Compute derivative matrices
- `const mfem::SparseMatrix & D(int i) const` — Get i-th first derivative
- `void SaveDerivToFile(const std::string &deriv, const std::string &filename) const`

### `FiniteDiff2d`

2D finite differences in `include/StreamVorti/finite_differences/fd_2d.hpp`.

**Derivative matrices:**
- `ShapeFunctionDx()` — d/dx
- `ShapeFunctionDy()` — d/dy
- `ShapeFunctionDxx()` — d²/dx²
- `ShapeFunctionDyy()` — d²/dy²
- `ShapeFunctionDxy()` — d²/dxdy

### `FiniteDiff3d`

3D finite differences in `include/StreamVorti/finite_differences/fd_3d.hpp`.

**Additional derivative matrices (beyond 2D):**
- `ShapeFunctionDz()` — d/dz
- `ShapeFunctionDzz()` — d²/dz²
- `ShapeFunctionDxz()` — d²/dxdz
- `ShapeFunctionDyz()` — d²/dydz

## Stencils

### First Derivatives

**Interior nodes** (central difference, 2nd-order):

```
df/dx ≈ (f(i+1,j) - f(i-1,j)) / (2h)
```

**Boundary nodes** (one-sided, 2nd-order):

```
Forward:  df/dx ≈ (-3f(i,j) + 4f(i+1,j) - f(i+2,j)) / (2h)
Backward: df/dx ≈ (3f(i,j) - 4f(i-1,j) + f(i-2,j)) / (2h)
```

### Second Derivatives

**Interior nodes** (central, 2nd-order):

```
d²f/dx² ≈ (f(i+1,j) - 2f(i,j) + f(i-1,j)) / h²
```

**Boundary nodes** (one-sided, 2nd-order):

```
Forward:  d²f/dx² ≈ (2f(i) - 5f(i+1) + 4f(i+2) - f(i+3)) / h²
Backward: d²f/dx² ≈ (2f(i) - 5f(i-1) + 4f(i-2) - f(i-3)) / h²
```

### Mixed Derivatives

**Interior nodes:**

```
d²f/dxdy ≈ (f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1)) / (4·hx·hy)
```

At boundaries, the stencil falls back to available one-sided neighbors.

## Usage Example

```cpp
#include "StreamVorti/mfem_main.hpp"

// Create a 40x40 structured mesh on [0,1]²
mfem::Mesh mesh(mfem::Mesh::MakeCartesian2D(
    40, 40, mfem::Element::QUADRILATERAL, false, 1.0, 1.0, false));

// H1 order-1 finite element space
mfem::H1_FECollection fec(1, 2);
mfem::FiniteElementSpace fes(&mesh, &fec, 1);
mfem::GridFunction u(&fes);

// Initialize u with some values...

// Build FD derivative matrices
StreamVorti::FiniteDiff2d fd(u);
fd.Update();

// Compute derivatives
mfem::Vector dudx(fes.GetNDofs());
mfem::Vector d2udx2(fes.GetNDofs());
fd.ShapeFunctionDx().Mult(u, dudx);
fd.ShapeFunctionDxx().Mult(u, d2udx2);

// Laplacian: ∇²u = d²u/dx² + d²u/dy²
mfem::Vector laplacian(fes.GetNDofs());
fd.ShapeFunctionDxx().Mult(u, laplacian);
mfem::Vector dyy_result(fes.GetNDofs());
fd.ShapeFunctionDyy().Mult(u, dyy_result);
laplacian += dyy_result;

// Save derivative matrix to file
fd.SaveDerivToFile("dx", "output_dat/fd_dx.txt");
```

## Comparison with DCPSE

| Feature | Finite Differences | DCPSE |
|---------|-------------------|-------|
| Grid requirement | Structured Cartesian | Any (scattered nodes) |
| Dependencies | MFEM only | MFEM + CGAL + Eigen |
| Assembly cost | O(N) | O(N·k²) where k = neighbors |
| Accuracy order | 2nd order | Depends on basis/neighbors |
| Boundary handling | One-sided stencils | Same as interior |
| Implementation | Analytical stencils | Weighted least squares |

## Requirements

- MFEM library
- Structured Cartesian mesh (created via `Mesh::MakeCartesian2D` or `MakeCartesian3D`)
- H1 finite element space, order 1

The mesh must have a regular grid structure. The module verifies this during
`Update()` and will abort with an error if the grid is not structured.

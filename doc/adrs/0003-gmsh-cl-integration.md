# ADR-0003: Gmsh-CL Integration for Mesh Generation

**Status:** Proposed
**Date:** 2026-03-23
**Authors:** Benjamin F. Zwick

---

## Context

StreamVorti currently supports only Cartesian meshes generated via MFEM and
meshes loaded from files. Complex geometries (e.g., flow around a cylinder)
require unstructured mesh generation. The `demo/cylinder.lisp` SDL file
referenced a `:generator :gmsh` syntax that was never implemented.

[gmsh-cl](https://github.com/benzwick/gmsh-cl) provides Common Lisp CFFI
bindings to Gmsh, with a clean API for OpenCASCADE geometry, meshing, and
visualization.

Three approaches were considered:

1. **C++ Gmsh API** — call Gmsh's C API from C++ after translating the Lisp
   geometry tree. Requires maintaining a geometry translator layer.
2. **StreamVorti geometry → gmsh-cl translator** — auto-translate StreamVorti's
   geometry classes (rectangle, circle, difference-shape) to gmsh-cl OCC calls.
   Adds an abstraction layer that duplicates gmsh-cl's already clean API.
3. **Direct gmsh-cl in SDL** — load gmsh-cl as an ASDF system in the embedded
   ECL runtime. Users write gmsh-cl code directly in SDL files.

## Decision

Option (3): load gmsh-cl in embedded ECL and expose it directly to SDL files
via a `(domain :gmsh ...)` macro form.

The simulation macro compiles `(domain :gmsh body...)` into code that:
1. Checks for the `GMSH` package (clear error if not available)
2. Calls `gmsh:initialize` / `gmsh:finalize` (with `unwind-protect`)
3. Adds a model and executes the user's gmsh-cl body forms
4. Writes the mesh to a temporary `.msh` file
5. Creates a `domain-data` with `:type :file`

The existing C++ `extractMesh` path handles the rest — `get-type` returns
`"loaded"`, `get-path` returns the file, and `MeshWrapper::loadFromFile`
loads the Gmsh mesh into MFEM (which natively reads `.msh` format).

gmsh-cl symbols are resolved at runtime via `INTERN` to avoid reader errors
when the `GMSH` package is not loaded (`STREAMVORTI_WITH_GMSH=OFF`).

## Consequences

- **No new abstraction layers.** gmsh-cl's API is the geometry API. Users write
  `(occ:rectangle ...)`, `(occ:cut ...)`, `(gmsh/mesh:generate ...)` directly.
- **Optional dependency.** Guarded by `STREAMVORTI_WITH_GMSH` CMake option.
  Requires ECL (`STREAMVORTI_WITH_ECL`), libgmsh at runtime, and CFFI for ECL.
- **gmsh-cl as git submodule** at `_reference/gmsh-cl`, consistent with MFEM
  and HYPRE.
- **CFFI + ECL risk.** CFFI lists ECL as supported but it is less tested than
  SBCL. If CFFI proves unreliable in embedded ECL, fallback would be ECL's
  native FFI.
- **No C++ changes required** for the basic integration. The existing
  `:file` → `"loaded"` → `loadFromFile` path works unchanged.

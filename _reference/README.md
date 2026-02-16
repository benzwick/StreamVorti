# Reference Sources

Git submodules and links to StreamVorti's dependencies for code reference.

## Submodules (checked out here)

| Library | Version | Directory | Description |
|---------|---------|-----------|-------------|
| [MFEM](https://mfem.org/) | 4.8 | `mfem/` | Finite element methods library. Provides mesh handling, linear algebra, grid functions. Core dependency. |
| [HYPRE](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) | 2.19.0 | `hypre/` | Scalable linear solvers and multigrid methods. Required by MFEM for parallel builds. |
| [vgplot](https://github.com/volkers/vgplot) | latest | `vgplot/` | Common Lisp gnuplot interface for plotting. Used for visualization scripts. |

To initialize submodules after cloning:
```bash
git submodule update --init
```

## Other Dependencies (not submoduled)

### Required

| Library | Version | Importance | Repository | Documentation |
|---------|---------|------------|------------|---------------|
| [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) | 4.0.3 | **Critical** - mesh partitioning for MFEM | [GitHub](https://github.com/KarypisLab/METIS) | [Manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf) |
| [CGAL](https://www.cgal.org/) | 5.0+ (6.0 recommended) | **Critical** - computational geometry algorithms, Delaunay triangulation, spatial searching | [GitHub](https://github.com/CGAL/cgal) | [Documentation](https://doc.cgal.org/latest/Manual/index.html) |
| [Eigen](https://eigen.tuxfamily.org/) | 3.4 | **Critical** - linear algebra (matrices, vectors, solvers) | [GitLab](https://gitlab.com/libeigen/eigen) | [Documentation](https://eigen.tuxfamily.org/dox/) |
| [OpenMP](https://www.openmp.org/) | (system) | **Required** - shared-memory parallelism | N/A | [Specification](https://www.openmp.org/specifications/) |

### Optional (Recommended)

| Library | Version | Importance | Repository | Documentation |
|---------|---------|------------|------------|---------------|
| [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) | (system) | **High** - UMFPACK/KLU direct solvers, default linear solver in StreamVorti | [GitHub](https://github.com/DrTimothyAldenDavis/SuiteSparse) | [UMFPACK User Guide](https://people.engr.tamu.edu/davis/suitesparse.html) |
| [ECL](https://ecl.common-lisp.dev/) | (system) | **High** - Embeddable Common Lisp, enables SDL simulation files | [GitLab](https://gitlab.com/embeddable-common-lisp/ecl) | [Manual](https://ecl.common-lisp.dev/static/manual/) |
| [SBCL](http://www.sbcl.org/) | (system) | **Medium** - Steel Bank Common Lisp, used for validation scripts and Lisp tests | [GitHub](https://github.com/sbcl/sbcl) | [Manual](http://www.sbcl.org/manual/) |
| [LAPACK](https://www.netlib.org/lapack/) | (system) | **Medium** - dense linear algebra routines, used by MFEM | [GitHub](https://github.com/Reference-LAPACK/lapack) | [Documentation](https://www.netlib.org/lapack/explore-html/) |
| [OpenMPI](https://www.open-mpi.org/) | (system) | **Medium** - MPI for distributed-memory parallelism | [GitHub](https://github.com/open-mpi/ompi) | [Documentation](https://www.open-mpi.org/doc/) |
| [zlib](https://zlib.net/) | (system) | **Low** - compression, used by MFEM for compressed mesh files | [GitHub](https://github.com/madler/zlib) | [Manual](https://zlib.net/manual.html) |

### Future Integration Candidates

| Library | Version | Purpose | Repository | Documentation |
|---------|---------|---------|------------|---------------|
| [Gmsh](https://gmsh.info/) | 4.13+ | Mesh generation with CSG support (needed for cylinder demo) | [GitLab](https://gitlab.onelab.info/gmsh/gmsh) | [Documentation](https://gmsh.info/doc/texinfo/gmsh.html), [C API](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-API) |
| [PETSc](https://petsc.org/) | 3.21+ | Scalable scientific computing (solvers, preconditioners) | [GitLab](https://gitlab.com/petsc/petsc) | [Documentation](https://petsc.org/release/manual/) |
| [Clasp](https://github.com/clasp-developers/clasp) | latest | Common Lisp on LLVM with direct C++ interop | [GitHub](https://github.com/clasp-developers/clasp) | [Wiki](https://github.com/clasp-developers/clasp/wiki) |
| [CFFI](https://cffi.common-lisp.dev/) | latest | Portable CL foreign function interface (ECL/SBCL/Clasp) | [GitHub](https://github.com/cffi/cffi) | [Manual](https://cffi.common-lisp.dev/manual/cffi-manual.html) |

## Build Versions Used in CI

From `.github/workflows/build-from-source.yml`:
- CMake 3.30+
- C++17 standard
- MFEM 4.8 with `+metis +lapack +suitesparse +zlib`
- HYPRE 2.19.0
- METIS 4.0.3

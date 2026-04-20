# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Dependencies — local build in `_reference/`

All major dependencies (MFEM, HYPRE) are built locally from the git submodules in `_reference/`. **Do not look in `/opt/mfem/` or other system locations** — the canonical build is local.

### 1. Build HYPRE (one-time)

```bash
cd _reference/hypre/src/cmbuild
cmake .. -DHYPRE_WITH_MPI=ON -DCMAKE_BUILD_TYPE=Release
make -j6
make install   # installs to _reference/hypre/src/hypre/
```

### 2. Build MFEM (one-time)

MFEM is built with MPI, MUMPS (parallel direct solver), SuperLU_DIST, and SuiteSparse (UMFPack).

The Debian system MUMPS does not link against ParMETIS, so we override
`MUMPS_REQUIRED_PACKAGES` to drop that dependency. Similarly, we provide
explicit paths for SuperLU_DIST and `MUMPS_LIBRARIES` because MFEM's default
search paths assume self-built tarballs in sibling directories. The explicit
`MUMPS_LIBRARIES` list also avoids MFEM's autodetection adding `MPI::MPI_C`
(unimported) as a transitive ScaLAPACK dependency.

```bash
cd _reference/mfem
mkdir build && cd build
cmake .. \
  -DMFEM_USE_MPI=YES \
  -DMFEM_USE_MUMPS=YES \
  -DMFEM_USE_SUPERLU=YES \
  -DMFEM_USE_SUITESPARSE=YES \
  -DHYPRE_DIR=$(realpath ../../hypre/src/hypre) \
  -DSuperLUDist_INCLUDE_DIRS=/usr/include/superlu-dist \
  -DSuperLUDist_LIBRARIES=/usr/lib/x86_64-linux-gnu/libsuperlu_dist.so \
  -DMUMPS_REQUIRED_PACKAGES="MPI_Fortran;METIS;LAPACK;BLAS" \
  -DMUMPS_LIBRARIES="/usr/lib/x86_64-linux-gnu/libdmumps.so;/usr/lib/x86_64-linux-gnu/libmumps_common.so;/usr/lib/x86_64-linux-gnu/libpord.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so" \
  -DCMAKE_BUILD_TYPE=Release
make -j6
```

### 3. Build StreamVorti

```bash
mkdir build && cd build
cmake -DMFEM_DIR=$(realpath ../_reference/mfem/build) -DCMAKE_BUILD_TYPE=Debug ..
make -j6
```

With ECL (enables SDL/Lisp simulation files):
```bash
cmake -DMFEM_DIR=$(realpath ../_reference/mfem/build) -DCMAKE_BUILD_TYPE=Debug -DSTREAMVORTI_WITH_ECL=ON ..
```

With Gmsh (enables unstructured mesh generation via gmsh-cl, requires ECL):
```bash
# 1. Initialize gmsh-cl submodule and install CL dependencies
git submodule update --init _reference/gmsh-cl
cd _reference/gmsh-cl && ocicl install

# 2. Build libgmsh from gmsh-cl's Gmsh submodule
cd _reference/gmsh-cl
git submodule update --init _reference/gmsh
cd _reference/gmsh && mkdir build && cd build
cmake .. -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=Release
make -j6

# 3. Build StreamVorti with Gmsh support
cd $PROJECT_ROOT/build
cmake -DMFEM_DIR=$(realpath ../_reference/mfem/build) -DCMAKE_BUILD_TYPE=Debug \
      -DSTREAMVORTI_WITH_ECL=ON -DSTREAMVORTI_WITH_GMSH=ON ..
make -j6
```
libgmsh is found automatically from the gmsh-cl submodule build directory.

Run tests:
```bash
cmake -DSTREAMVORTI_BUILD_TESTS=ON ..
make -j6
ctest --output-on-failure   # C++ tests (Google Test)
make test-lisp              # Lisp tests (requires SBCL or ECL)
```

The build produces:
- `StreamVorti` - Full solver (from `streamvorti.cpp`) - time-stepping, ParaView output, SDL support
- `MfemRun` - Derivative export tool (from `mfem_main.cpp`) - outputs DCPSE derivatives for MATLAB/Python
- `libStreamVorti_static.a` - Static library in `build/lib/StreamVorti/`

## System dependencies (apt packages)

**Required:**
- `libcgal-dev` — Computational Geometry Algorithms Library
- `libeigen3-dev` — linear algebra
- `libopenmpi-dev` — MPI
- `libsuitesparse-dev` — UMFPack and friends (for MFEM_USE_SUITESPARSE)
- `libsuperlu-dist-dev` — distributed parallel sparse direct solver (for MFEM_USE_SUPERLU)
- `libmumps-dev` — parallel sparse direct solver (for MFEM_USE_MUMPS)
- `libscalapack-openmpi-dev` — required by MUMPS
- `libmetis-dev` — mesh partitioning

Install everything in one go:
```bash
sudo apt install libcgal-dev libeigen3-dev libopenmpi-dev libsuitesparse-dev \
                 libsuperlu-dist-dev libmumps-dev libscalapack-openmpi-dev \
                 libmetis-dev
```

OpenMP comes with the system compiler.

**Debian Linux Note:** If you encounter `fatal error: Eigen/Dense: No such file or directory`, create symlinks:
```bash
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
```

## Architecture Overview

StreamVorti is a C++ library for solving PDEs using explicit methods, specifically implementing the Strong-Form Meshless Stream Function - Vorticity formulation.

**Core Components:**
- `src/approximants/` - DCPSE (Direct Collocation Particle Strength Exchange) implementations
  - `dcpse.hpp/cpp` - Base DCPSE functionality
  - `dcpse_2d.hpp/cpp` - 2D-specific DCPSE methods
  - `dcpse_3d.hpp/cpp` - 3D-specific DCPSE methods
- `src/support_domain/` - Support domain calculations for meshless methods
- `include/StreamVorti/mfem_main.hpp` - Main header that includes all modules

**Executables:**

1. `StreamVorti` (from `streamvorti.cpp`) - Full solver
   - CLI mode: `StreamVorti -dim 2 -nx 40 -ny 40 -nn 25 -Re 100 -dt 0.001 -ft 10.0 -pv`
   - SDL mode: `StreamVorti -f demo/cavity.lisp -lp lisp -pv` (requires ECL build)
   - `-lp lisp` tells the binary where to find the Lisp SDL files; needed when not running from repo root

2. `MfemRun` (from `mfem_main.cpp`) - Derivative export tool
   - Usage: `MfemRun -sd -sn`
   - Outputs DCPSE derivative matrices for use by MATLAB/Python

**Library Structure:**
- Header-only style with implementation in `src/` compiled into static library
- Modular design with each component as separate CMake target
- All functionality accessible through single `mfem_main.hpp` header

**Reference Sources:**
- `_reference/` contains git submodules of key dependencies (MFEM, HYPRE, vgplot) for code reference
- `_reference/README.md` has links and version info for all dependencies
- Initialize with: `git submodule update --init`

The codebase focuses on meshless methods for computational fluid dynamics, with MFEM providing mesh handling and linear algebra operations while StreamVorti implements the meshless approximation schemes.

## Validation

Validate results against Ghia et al. (1982) benchmark data:

1. Run simulation to steady state (CLI mode):
```bash
./StreamVorti -nx 40 -ny 40 -Re 100 -dt 0.001 -ft 30.0 -pv
```

2. Run validation script (requires SBCL):
```bash
# CLI mode output uses prefix "mfem_square10x10":
sbcl --load lisp/compare-ghia.lisp \
     --eval '(streamvorti.validation:run-validation :results-file "output_dat/mfem_square10x10_u_centerline_x0.5.dat" :reynolds 100)'

# SDL mode output uses the simulation name (e.g., "lid-driven-cavity"):
sbcl --load lisp/compare-ghia.lisp \
     --eval '(streamvorti.validation:run-validation :results-file "output_dat/lid-driven-cavity_u_centerline_x0.5.dat" :reynolds 100)'
```

## Git Workflow

**CRITICAL:** Never use `git add -A` or `git add .`. Always follow these steps:

1. `git diff` - Review changes first
2. `git add <specific-files>` - Stage specific files
3. `git commit -m "message"` - Commit with descriptive message

Each step must be separate commands, not chained.

**Separate commits for unrelated changes:** Always create separate commits for logically unrelated changes. Each commit should be atomic and address a single concern.

**NEVER debug via CI commits.** Do not push speculative fixes to see if they work in CI. This pollutes the git history with garbage commits. Instead:
1. Test locally first
2. Read the documentation and source code to understand how things work
3. If stuck, ask the user for help rather than guessing

## Code Quality

**No hacks or workarounds.** When using a library:
1. Read the documentation and source code first
2. Look at working examples (e.g., demo files, tests)
3. Write clean, idiomatic code that uses the library as intended
4. If something doesn't work, understand WHY before attempting fixes
5. Never add sleeps, polling loops, or other band-aids to paper over misunderstanding
6. Ask for help rather than piling on hacks
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

Build the project using CMake:
```bash
mkdir build
cd build
cmake -DMFEM_DIR=/opt/mfem/mfem-4.8 -DCMAKE_BUILD_TYPE=Debug ..
make -j6
```

The build produces two executables:
- `MfemRun` - Basic MFEM demo (from `mfem_main.cpp`) - computes and saves derivatives only
- `StreamVorti` - Full lid-driven cavity solver (from `streamvorti.cpp`) - time-stepping solver with ParaView output, validation features
- `libStreamVorti_static.a` - Static library in `build/lib/StreamVorti/`

**Use `StreamVorti` for the complete lid-driven cavity solver.**

## Dependencies

**Required:**
- MFEM library (specify path with `-DMFEM_DIR`)
- CGAL (Computational Geometry Algorithms Library)
- Eigen3 (linear algebra library)

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

1. `MfemRun` (from `mfem_main.cpp`) - Basic derivative computation demo
   - Usage: `MfemRun -sd -sn`
   - Computes and saves DCPSE derivative matrices
   - No time-stepping solver

2. `StreamVorti` (from `streamvorti.cpp`) - Full lid-driven cavity solver
   - Usage: `StreamVorti -dim 2 -nx 40 -ny 40 -nn 25 -Re 100 -dt 0.001 -ft 10.0 -pv`
   - Implements stream function-vorticity formulation
   - Time-stepping solver with steady-state detection
   - ParaView VTK output
   - Centerline extraction for Ghia validation
   - Optional SDL/Lisp support with `-DSTREAMVORTI_WITH_ECL=ON` (requires `ecl libecl-dev`)
   - SDL usage: `StreamVorti -f demo/cavity.lisp -pv`

**Configuration:**
- Command-line options for grid, Re, dt, etc.
- Optional `.ini` files for DCPSE parameters
- SDL files (`.lisp`) for complete simulation definition (requires ECL)

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

1. Run simulation to steady state:
```bash
./StreamVorti -nx 40 -ny 40 -Re 100 -dt 0.001 -ft 30.0 -pv
```

2. Run validation script (requires SBCL):
```bash
sbcl --load lisp/compare-ghia.lisp \
     --eval '(streamvorti.validation:run-validation :results-file "output_dat/mfem_square10x10_u_centerline_x0.5.dat" :reynolds 100)'
```

Reference data files in `data/`:
- `ghia_1982_u.txt` - u-velocity along x=0.5 (Re=100-10000)
- `ghia_1982_v.txt` - v-velocity along y=0.5
- `erturk_2005_u.txt` - u-velocity for high Re (1000-21000)

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
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

The build produces:
- `MfemRun` - Main MFEM simulation executable
- `libStreamVorti_static.a` - Static library in `build/lib/StreamVorti/`

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
- `include/StreamVorti/stream_vorti.hpp` - Main header that includes all modules

**Main Executable:**
- `mfem_main.cpp` - MFEM integration example/demo
- Usage: `MfemRun -sd -sn`

**Configuration:**
- Uses `.ini` files for configuration (see `demo/mfem.ini`)
- DCPSE parameters: cutoff radius, support radius, output file specifications
- Config includes neighbor computation and derivative approximation settings

**Library Structure:**
- Header-only style with implementation in `src/` compiled into static library
- Modular design with each component as separate CMake target
- All functionality accessible through single `stream_vorti.hpp` header

The codebase focuses on meshless methods for computational fluid dynamics, with MFEM providing mesh handling and linear algebra operations while StreamVorti implements the meshless approximation schemes.
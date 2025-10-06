#!/bin/bash
set -e

MFEM_DIR=${1:-"/opt/mfem/mfem-4.8"}
BUILD_TYPE=${2:-"Release"}

echo "Building StreamVorti with MFEM from ${MFEM_DIR}..."

# Create build directory
mkdir -p build
cd build

# Configure StreamVorti with MFEM
cmake \
  -DMFEM_DIR=${MFEM_DIR} \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DCMAKE_CXX_FLAGS="-DOMPI_SKIP_MPICXX" \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  ..

# Build
make -j$(nproc)

echo "StreamVorti built successfully"

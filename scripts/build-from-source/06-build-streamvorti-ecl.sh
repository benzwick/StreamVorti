#!/bin/bash
set -e

MFEM_DIR=${1:-"/opt/mfem/mfem-4.8"}
BUILD_TYPE=${2:-"Debug"}
ENABLE_COVERAGE=${3:-"ON"}

echo "Building StreamVorti with ECL support..."
echo "  MFEM_DIR: ${MFEM_DIR}"
echo "  BUILD_TYPE: ${BUILD_TYPE}"
echo "  ENABLE_COVERAGE: ${ENABLE_COVERAGE}"

# Create build directory
mkdir -p build
cd build

# Configure StreamVorti with ECL
cmake \
  -DMFEM_DIR=${MFEM_DIR} \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DSTREAMVORTI_WITH_ECL=ON \
  -DSTREAMVORTI_BUILD_TESTS=ON \
  -DENABLE_COVERAGE=${ENABLE_COVERAGE} \
  -DCMAKE_CXX_FLAGS="-DOMPI_SKIP_MPICXX" \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  ..

# Build
make -j$(nproc)

echo "StreamVorti with ECL built successfully"

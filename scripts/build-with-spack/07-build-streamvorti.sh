#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}
ENV_NAME=${2:-"streamvorti"}
BUILD_TYPE=${3:-"Release"}

echo "Building StreamVorti with Spack environment '${ENV_NAME}'..."

source ${SPACK_DIR}/share/spack/setup-env.sh
spack env activate ${ENV_NAME}

# Get MFEM installation path
MFEM_DIR=$(spack location -i mfem)
echo "MFEM_DIR=$MFEM_DIR"

# Create build directory
mkdir -p build
cd build

# Configure with CMake from Spack environment
# The environment view provides unified paths for all dependencies
cmake \
  -DMFEM_DIR=$MFEM_DIR \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  ..

# Build
make -j$(nproc) VERBOSE=1

echo "StreamVorti built successfully"

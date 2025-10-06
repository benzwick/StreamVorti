#!/bin/bash
set -e

HYPRE_URL=${1:-"https://github.com/hypre-space/hypre/archive"}
HYPRE_ARCHIVE=${2:-"v2.19.0"}
HYPRE_TOP_DIR=${3:-"hypre-2.19.0"}

echo "Building HYPRE ${HYPRE_ARCHIVE}..."

# Download and extract HYPRE
wget ${HYPRE_URL}/${HYPRE_ARCHIVE}.tar.gz
tar -xzf ${HYPRE_ARCHIVE}.tar.gz
cd ${HYPRE_TOP_DIR}/src

# Configure and build HYPRE
./configure --prefix=$PWD/../install CC=mpicc CXX=mpic++
make -j$(nproc)
make install

echo "HYPRE built successfully at ${HYPRE_TOP_DIR}/install"

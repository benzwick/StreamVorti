#!/bin/bash
set -e

MFEM_TOP_DIR=${1:-"mfem-4.8"}
HYPRE_TOP_DIR=${2:-"hypre-2.19.0"}
METIS_TOP_DIR=${3:-"metis-4.0.3"}
INSTALL_PREFIX=${4:-"/opt/mfem/${MFEM_TOP_DIR}"}

echo "Building MFEM ${MFEM_TOP_DIR}..."

# Clone MFEM repository
git clone https://github.com/mfem/mfem.git ${MFEM_TOP_DIR}
cd ${MFEM_TOP_DIR}
git checkout v4.8
cd ..

# Create build directory
mkdir mfem-build
cd mfem-build

# Configure MFEM with CMake
cmake \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMFEM_USE_MPI=YES \
  -DHYPRE_DIR=$PWD/../${HYPRE_TOP_DIR}/install \
  -DMFEM_USE_METIS=YES \
  -DMETIS_DIR=$PWD/../${METIS_TOP_DIR} \
  -DMFEM_USE_LAPACK=YES \
  -DMFEM_USE_SUITESPARSE=YES \
  -DMFEM_USE_ZLIB=YES \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  ../${MFEM_TOP_DIR}

# Build and install
make -j$(nproc)
sudo make install

# Clean up source directories to save space
cd ..
rm -rf ${MFEM_TOP_DIR} mfem-build

echo "MFEM installed successfully at ${INSTALL_PREFIX}"

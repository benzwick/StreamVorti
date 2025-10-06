#!/bin/bash
set -e

METIS_URL=${1:-"https://mfem.github.io/tpls"}
METIS_ARCHIVE=${2:-"metis-4.0.3.tar.gz"}
METIS_TOP_DIR=${3:-"metis-4.0.3"}

echo "Building METIS..."

# Download and extract METIS
wget ${METIS_URL}/${METIS_ARCHIVE}
tar -xzf ${METIS_ARCHIVE}
cd ${METIS_TOP_DIR}

# Build METIS library
make -C Lib CC=mpicc OPTFLAGS="-Wno-error=implicit-function-declaration -O2"

echo "METIS built successfully at ${METIS_TOP_DIR}/libmetis.a"

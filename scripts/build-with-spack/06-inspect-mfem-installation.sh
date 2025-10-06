#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}
ENV_NAME=${2:-"streamvorti"}

echo "Inspecting MFEM installation..."

source ${SPACK_DIR}/share/spack/setup-env.sh
spack env activate ${ENV_NAME}

# Get installation path
MFEM_DIR=$(spack location -i mfem)
echo "========================================="
echo "MFEM Installation Information"
echo "========================================="
echo "MFEM_DIR=$MFEM_DIR"
echo ""

echo "--- Spack MFEM Info ---"
spack find -vl mfem
echo ""

echo "--- Directory Structure (tree-style) ---"
find "$MFEM_DIR" -print | sed -e "s;$MFEM_DIR;.;g;s;[^/]*\/;|__;g;s;__|; |;g"
echo ""

echo "--- Full Recursive Listing ---"
ls -laR "$MFEM_DIR"
echo ""

echo "--- All MFEM-related files ---"
find "$MFEM_DIR" -type f -name "*mfem*" -o -name "*MFEM*"
echo ""

echo "--- All .cmake files ---"
find "$MFEM_DIR" -type f -name "*.cmake"
echo ""

echo "--- All executable files ---"
find "$MFEM_DIR" -type f -executable
echo ""

echo "--- config.mk contents ---"
cat "$MFEM_DIR/share/mfem/config.mk"
echo ""

echo "========================================="

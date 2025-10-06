#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}
ENV_NAME=${2:-"streamvorti"}

echo "Showing Spack environment information for '${ENV_NAME}'..."

source ${SPACK_DIR}/share/spack/setup-env.sh
spack env activate ${ENV_NAME}

echo "=== Spack Environment ==="
spack env status
echo ""
echo "=== Environment Specs ==="
spack find
echo ""
echo "=== MFEM Configuration ==="
spack spec mfem || true

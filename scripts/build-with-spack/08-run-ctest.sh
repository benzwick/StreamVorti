#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}
ENV_NAME=${2:-"streamvorti"}
BUILD_DIR=${3:-"build"}

echo "Running ctest in ${BUILD_DIR} with Spack environment '${ENV_NAME}'..."

source ${SPACK_DIR}/share/spack/setup-env.sh
spack env activate ${ENV_NAME}

cd ${BUILD_DIR}
ctest --output-on-failure || true

echo "ctest completed"

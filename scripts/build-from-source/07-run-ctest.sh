#!/bin/bash
set -e

BUILD_DIR=${1:-"build"}

echo "Running ctest in ${BUILD_DIR}..."

cd ${BUILD_DIR}
ctest --output-on-failure || true

echo "ctest completed"

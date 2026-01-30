#!/bin/bash
set -e

BUILD_DIR=${1:-"build"}

echo "Generating coverage report in ${BUILD_DIR}..."

cd ${BUILD_DIR}

# Capture coverage data
lcov --capture --directory . --output-file coverage.info

# Remove external code from coverage
lcov --remove coverage.info '/usr/*' '*/tests/*' --output-file coverage.info

# Display coverage summary
lcov --list coverage.info

echo "Coverage report generated: ${BUILD_DIR}/coverage.info"

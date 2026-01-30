#!/bin/bash
set -e

BUILD_DIR=${1:-"build"}

echo "Generating coverage report in ${BUILD_DIR}..."

cd ${BUILD_DIR}

# Capture coverage data
# Use --ignore-errors to handle gcov version mismatches
lcov --capture --directory . --output-file coverage.info --ignore-errors mismatch,gcov,source

# Remove external code from coverage
lcov --remove coverage.info '/usr/*' '*/tests/*' --output-file coverage.info --ignore-errors unused

# Display coverage summary
lcov --list coverage.info --ignore-errors unused || true

echo "Coverage report generated: ${BUILD_DIR}/coverage.info"

#!/bin/bash
set -e

BUILD_DIR=${1:-"build"}

echo "Generating coverage report in ${BUILD_DIR}..."

cd ${BUILD_DIR}

# Check lcov version to determine which --ignore-errors options are supported
LCOV_VERSION=$(lcov --version | grep -oP '\d+\.\d+' | head -1)
LCOV_MAJOR=$(echo $LCOV_VERSION | cut -d. -f1)

# Build ignore-errors flags based on lcov version
# Version 2.0+ supports: mismatch,empty,unused
# Older versions only support: gcov,source,graph
if [ "$LCOV_MAJOR" -ge 2 ]; then
    IGNORE_ERRORS="--ignore-errors mismatch,gcov,source,unused"
else
    IGNORE_ERRORS="--ignore-errors gcov,source"
fi

# Capture coverage data
lcov --capture --directory . --output-file coverage.info $IGNORE_ERRORS || true

# Remove external code from coverage
lcov --remove coverage.info '/usr/*' '*/tests/*' --output-file coverage.info $IGNORE_ERRORS || true

# Display coverage summary
lcov --list coverage.info $IGNORE_ERRORS || true

echo "Coverage report generated: ${BUILD_DIR}/coverage.info"

#!/bin/bash
set -e

BUILD_DIR=${1:-"build"}

echo "Checking built executables and libraries in ${BUILD_DIR}..."

ls -la ${BUILD_DIR}/
ls -la ${BUILD_DIR}/lib/StreamVorti/ || true

if [ -f ${BUILD_DIR}/MfemRun ]; then
  echo "MfemRun executable built successfully"
  ldd ${BUILD_DIR}/MfemRun
else
  echo "WARNING: MfemRun executable not found"
fi

if [ -f ${BUILD_DIR}/StreamVorti ]; then
  echo "StreamVorti executable built successfully"
  ldd ${BUILD_DIR}/StreamVorti
else
  echo "WARNING: StreamVorti executable not found"
fi

if [ -f ${BUILD_DIR}/lib/StreamVorti/libStreamVorti_static.a ]; then
  echo "StreamVorti static library built successfully"
  file ${BUILD_DIR}/lib/StreamVorti/libStreamVorti_static.a
else
  echo "WARNING: StreamVorti static library not found"
fi

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
  echo "ERROR: MfemRun executable not found"
  exit 1
fi

if [ -f ${BUILD_DIR}/StreamVorti ]; then
  echo "StreamVorti executable built successfully"
  ldd ${BUILD_DIR}/StreamVorti
else
  echo "ERROR: StreamVorti executable not found"
  exit 1
fi

if [ -f ${BUILD_DIR}/StreamVorti_par ]; then
  echo "StreamVorti_par executable built successfully"
  ldd ${BUILD_DIR}/StreamVorti_par
else
  echo "WARNING: StreamVorti_par executable not found"
fi

if [ -f ${BUILD_DIR}/lib/StreamVorti/libStreamVorti_static.a ]; then
  echo "StreamVorti static library built successfully"
  file ${BUILD_DIR}/lib/StreamVorti/libStreamVorti_static.a
else
  echo "WARNING: StreamVorti static library not found"
fi

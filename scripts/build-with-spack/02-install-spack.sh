#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}

echo "Installing Spack to ${SPACK_DIR}..."

# Check if Spack already exists
if [ -d "${SPACK_DIR}" ]; then
  echo "Spack directory already exists at ${SPACK_DIR}"
  echo "To reinstall, remove the directory first: rm -rf ${SPACK_DIR}"
else
  git clone -c feature.manyFiles=true https://github.com/spack/spack.git ${SPACK_DIR}
  echo "source ${SPACK_DIR}/share/spack/setup-env.sh" >> ~/.bashrc
  echo "Spack installed successfully"
fi

echo "To use Spack in this session: source ${SPACK_DIR}/share/spack/setup-env.sh"

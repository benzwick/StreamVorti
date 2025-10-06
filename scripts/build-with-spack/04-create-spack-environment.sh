#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}
ENV_NAME=${2:-"streamvorti"}
SPACK_YAML=${3:-"spack.yaml"}

echo "Creating Spack environment '${ENV_NAME}' from ${SPACK_YAML}..."

source ${SPACK_DIR}/share/spack/setup-env.sh

# Create environment from spack.yaml in repo
if spack env list | grep -q ${ENV_NAME}; then
  echo "Environment '${ENV_NAME}' already exists, activating..."
  spack env activate ${ENV_NAME}
else
  echo "Creating environment from ${SPACK_YAML}..."
  spack env create ${ENV_NAME} ${SPACK_YAML}
  spack env activate ${ENV_NAME}
fi

# Install all packages in environment
echo "Installing packages in environment..."
spack install

echo "Spack environment '${ENV_NAME}' created and installed successfully"

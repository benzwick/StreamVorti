#!/bin/bash
set -e

SPACK_DIR=${1:-"$HOME/spack"}

echo "Configuring Spack..."

source ${SPACK_DIR}/share/spack/setup-env.sh

# Configure Spack to use system compilers
spack compiler find

# Add binary mirror for faster installs
spack mirror add binary_mirror https://binaries.spack.io/releases/v0.21
spack buildcache keys --install --trust

echo "Spack configured successfully"

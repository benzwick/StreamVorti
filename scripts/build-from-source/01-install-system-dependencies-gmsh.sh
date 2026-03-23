#!/bin/bash
set -e

echo "Installing system dependencies (with ECL + Gmsh support)..."

sudo apt-get update
sudo apt-get install -y \
  build-essential \
  cmake \
  ccache \
  libopenmpi-dev \
  openmpi-bin \
  liblapack-dev \
  libsuitesparse-dev \
  zlib1g-dev \
  libcgal-dev \
  libeigen3-dev \
  ecl \
  libecl-dev \
  sbcl \
  lcov \
  libgtest-dev \
  libfltk1.3-dev \
  libocct-modeling-algorithms-dev \
  libocct-modeling-data-dev \
  libocct-data-exchange-dev \
  libocct-foundation-dev \
  libgl-dev \
  libglu1-mesa-dev

# Install ocicl from GitHub release (not in Ubuntu apt repos)
OCICL_VERSION=2.16.5
wget -q "https://github.com/ocicl/ocicl/releases/download/v${OCICL_VERSION}/ocicl_${OCICL_VERSION}-1_amd64.deb"
sudo apt install -y "./ocicl_${OCICL_VERSION}-1_amd64.deb"
rm -f "ocicl_${OCICL_VERSION}-1_amd64.deb"

# Configure ocicl for SBCL (needed by ocicl install)
echo '(require :asdf)' > ~/.sbclrc
ocicl setup >> ~/.sbclrc

echo "System dependencies (with ECL + Gmsh) installed successfully"

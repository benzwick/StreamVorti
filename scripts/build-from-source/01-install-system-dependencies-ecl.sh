#!/bin/bash
set -e

echo "Installing system dependencies (with ECL support)..."

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
  lcov \
  libgtest-dev

echo "System dependencies (with ECL) installed successfully"

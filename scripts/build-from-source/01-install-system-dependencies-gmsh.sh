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
  ocicl \
  lcov \
  libgtest-dev \
  libfltk1.3-dev \
  libocct-modeling-algorithms-dev \
  libocct-modeling-data-dev \
  libocct-data-exchange-dev \
  libocct-foundation-dev \
  libgl-dev \
  libglu1-mesa-dev

echo "System dependencies (with ECL + Gmsh) installed successfully"

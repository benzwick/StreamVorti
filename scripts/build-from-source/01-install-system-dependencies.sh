#!/bin/bash
set -e

echo "Installing system dependencies..."

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
  libeigen3-dev

echo "System dependencies installed successfully"

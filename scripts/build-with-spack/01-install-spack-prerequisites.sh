#!/bin/bash
set -e

echo "Installing minimal system prerequisites for Spack..."

sudo apt-get update
sudo apt-get install -y \
  build-essential \
  ca-certificates \
  coreutils \
  curl \
  environment-modules \
  gfortran \
  git \
  gpg \
  lsb-release \
  python3 \
  python3-venv \
  unzip \
  zip

echo "Spack prerequisites installed successfully"

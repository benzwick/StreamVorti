#!/bin/bash
set -e

echo "Setting up ccache..."

# Create Eigen symlinks (Debian-specific fix)
sudo ln -sf /usr/include/eigen3/Eigen /usr/include/Eigen 2>/dev/null || true
sudo ln -sf /usr/include/eigen3/unsupported /usr/include/unsupported 2>/dev/null || true

echo "Eigen symlinks created"
echo "ccache setup complete (actual ccache config done by GitHub Actions)"

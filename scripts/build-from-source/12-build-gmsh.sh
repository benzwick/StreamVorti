#!/bin/bash
set -e

echo "Building gmsh-cl and libgmsh..."

# Initialize gmsh-cl submodule
git submodule update --init _reference/gmsh-cl

# Install gmsh-cl CL dependencies via ocicl
cd _reference/gmsh-cl
ocicl install

# Initialize and build Gmsh from gmsh-cl's submodule
git submodule update --init _reference/gmsh
cd _reference/gmsh
mkdir -p build
cd build
cmake .. -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

echo "libgmsh built at: $(pwd)/libgmsh.so"
echo "gmsh-cl and libgmsh built successfully"

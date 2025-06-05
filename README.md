# StreamVorti

# Installing

## Install MFEM

### MacOS

Parallel build:

    cd <mfem-root-dir>

   (download hypre and METIS 4 from above URLs)
   (build METIS 4 in ../metis-4.0 relative to mfem/)
   (build hypre in ../hypre relative to mfem/)

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/opt/mfem/mfem-4.8 -DMFEM_USE_MPI=YES  -DCMAKE_BUILD_TYPE=Debug ..
    make -j 8                  # Build MFEM



    make tests -j 8            # Build unit-tests
    make examples -j 8         # Build examples
    make miniapps -j 8         # Build miniapps
    make test                  # Run the tests

    make install

    Make a note of where MFEM is installed.

    cmake -DMETIS_DIR=/my/METIS


## Install StreamVorti

    cmake -DMFEM_DIR=/opt/mfem/mfem-4.8 -DCMAKE_BUILD_TYPE=Debug ..
    make -j6

Note: On Debian Linux, if you encounter `fatal error: Eigen/Dense: No such file or directory`, create symlinks:

    cd /usr/include
    sudo ln -sf eigen3/Eigen Eigen
    sudo ln -sf eigen3/unsupported unsupported

# Usage

```
MfemRun -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -sm -sn -sd -sdd
```

Read derivatives saved by streamvorti mfem into matlab
check this first:

x: (in MATLAB)
load 'mfem_square10x10.dx.dat'
fx=spconvert(mfem_square10x10_dx)

same for y, xx, yy


TODO:
output the mesh nodes similar to matlab geometry

Read mesh nodes into matlab from file
nodes=load('mfem_square10x10.mesh.geom.dat');
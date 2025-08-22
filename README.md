# StreamVorti

# Installing

## Install MFEM

### MacOS

Parallel build:
```
cd <mfem-root-dir>
```

   (download hypre and METIS 4 from above URLs)
   (build METIS 4 in ../metis-4.0 relative to mfem/)
   (build hypre in ../hypre relative to mfem/)

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/opt/mfem/mfem-4.8 -DMFEM_USE_MPI=YES  -DCMAKE_BUILD_TYPE=Debug ..
make -j 8                  # Build MFEM
```

```
make tests -j 8            # Build unit-tests
make examples -j 8         # Build examples
make miniapps -j 8         # Build miniapps
make test                  # Run the tests

make install
```
Make a note of where MFEM is installed.


## Install StreamVorti

```
mkdir build
cd build
cmake -DMFEM_DIR=/opt/mfem/mfem-4.8 -DCMAKE_BUILD_TYPE=Debug ..
make -j6
```
Note: On Debian Linux, if you encounter `fatal error: Eigen/Dense: No such file or directory`, create symlinks:
```
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
```

Note: On MacOS, you might need to add this option above to find the Eigen headers in the home directory:
```
-DEIGEN3_INCLUDE_DIRS=~/include
```

# Usage

```
MfemRun -dim 2 -sx 1 -sy 1 -nx 40 -ny 40 -sm -sn -sd -sdd
```

Read derivatives saved by streamvorti mfem into matlab
check this first:

x: (in MATLAB)
```
load 'mfem_square10x10.dx.dat'
fx=spconvert(mfem_square10x10_dx)
```
same for y, xx, yy

## OpenMP

For macOS, first do:
```
brew install libomp
```
and set in ~/.zshrc
```
export OpenMP_ROOT=$(brew --prefix)/opt/libomp
```


# Kaya HPC System

## Step 1 Load environment modules using module

https://hpc-wiki.info/hpc/Modules

module avail - show available modules

module list - show modules currently loaded in your environment

module load app/ver - load environment modules for version ver of app

module unload app/ver - unload an application module

Pre-installed modules on Kaya:
```
module load gcc/13.3.0 cmake/3.25.3  boost/1.84.0 eigen/3.4.0 metis/5.1.0 openmpi/4.1.6
```

## Install Spack and packages:
Use Spack to install and manage modules/packages that are not available on Kaya
(Search for package on https://packages.spack.io/)
In HOME directory:
```
git clone --depth=2 https://github.com/spack/spack.git
```

Sourcing Spack setup script file (for bash/zsh/sh)
```
. spack/share/spack/setup-env.sh
```

Create and activate Spack environment:
```
spack env create myenv
spack env list
spack env activate -p myenv
```
(Note: Env activation only works on Login node, activating on Compute node will cause detachment from it)

Add and Install packages to the active env
```
spack add <package_name> <names...>
spack install
```
Note: If error occurs when intalling "diffutils-3.10", try installing on a Compute node in interactive mode(salloc)

cgal 6.0.1
hypre 2.33.0 +amdgpu_target +mpi +gpu-aware-mpi +openmp +superlu-dist
(suite-sparse)

mfem@4.7.0+metis+suite-sparse+openmp


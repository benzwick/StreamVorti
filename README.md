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
cmake -DMFEM_DIR=/home/wli/spack/opt/spack/linux-cascadelake/mfem-4.7.0-undz3pa6cop3pvvmmbczaebgsbo4rqxp -DCMAKE_BUILD_TYPE=Debug ..
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
cmake -DMFEM_DIR=/home/wli/spack/opt/spack/linux-cascadelake/mfem-4.7.0-undz3pa6cop3pvvmmbczaebgsbo4rqxp -DCMAKE_BUILD_TYPE=Debug ..
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

### Step 0 (Optional): Access a Compute node
Allocate a Compute node and access it interactively to avoid later installation process getting killed on a Login node due to its limits. For example of parallel compilation jobs:
```
salloc -p work --ntasks=1 --cpus-per-task=8 --ntasks-per-core=1  --mem=8G --time=2:00:00
# Then can use: make -j8
```

### Step 1: Load environment modules using module

https://hpc-wiki.info/hpc/Modules

```
module avail - show available modules
module list - show modules currently loaded in your environment
module load app/ver - load environment modules for version ver of app
module unload app/ver - unload an application module

```
Available modules on Kaya:
```
module load gcc/13.3.0 cmake/3.25.2 metis/5.1.0 boost/1.84.0 eigen/3.4.0
```


### Step2: Install Spack and packages
- Use Spack to install and manage modules/packages that are not available on Kaya system.
- Search for package on https://packages.spack.io/

```
git clone --depth=2 https://github.com/spack/spack.git
```

- For bash/zsh/sh
```
. spack/share/spack/setup-env.sh
```
- Create and activate Spack environment:
```
spack env create myenv
spack env list
spack env activate -p myenv
spack env status
spack env deactivate
spack env remove myenv
```
(Note: Env activation only works on Login node, activating on Compute node will cause detachment from it)

- Find compiler (gcc) for new env
```
spack compiler find
```
will find and add compilar to the current env.

- Add MFEM package to spec
```
spack add mfem@4.7

```

- Install all added packages (listed in the packages.yaml)
```
spack install
```
Note: If error occurs when intalling "diffutils-3.10", try installing on a Compute node in interactive mode(salloc)

Note:

- If spack install job gets killed on login node, try running on a compute node.

- If error occurs when intalling "diffutils-3.10", try installing on a Compute node in interactive mode(salloc).

- "Error: cmake-3.31.8-b3tbbwlqmj4z5gioxpb24reernudaun6: AttributeError: 'super' object has no attribute 'dag_hash'"

```
cgal 6.0.1
```
- Note: requires Cmake version 3.22 or later. https://doc.cgal.org/latest/Manual/thirdparty.html
-  Module cmake/3.30.2 is the only available version on Kaya. (PS. cmake/3.30.2 is currently broken, use /3.25.2 instead)

```
cgal 5.6 is the next available version on Spack
```





```
hypre 2.33.0 +amdgpu_target +mpi +gpu-aware-mpi +openmp +superlu-dist
(suite-sparse)
```



/home/wli/spack/opt/spack/linux-cascadelake/mfem-4.7.0-undz3pa6cop3pvvmmbczaebgsbo4rqxp

### Using SLURM
Commonly used commands:
```
sinfo
sbatch <batch script>
squeue -u <your username>
srun <resource-parameters>
scancel <jobid>
sacct -j <job id>
```

# StreamVorti

# Installing

    cmake -DMFEM_DIR=/opt/mfem/mfem-4.5 -DCMAKE_BUILD_TYPE=Debug ..
    make -j6

Note: On Debian Linux, if you encounter `fatal error: Eigen/Dense: No such file or directory`, create symlinks:

    cd /usr/include
    sudo ln -sf eigen3/Eigen Eigen
    sudo ln -sf eigen3/unsupported unsupported

# TODO

- restructure following similar structure as MFEM
  examples/poisson
  miniapps/dcpse
  src/efm
  src/general

- replace CGAL with nanoflann or similar
  - https://stackoverflow.com/questions/15124900/why-are-kd-trees-so-damn-slow-for-nearest-neighbor-search-in-point-sets
  - https://pointclouds.org/

- replace config manager with MFEM ArgParser (no more config files)

- rename library and replace headers

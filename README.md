# StreamVorti

# Installing

    cmake -DMFEM_DIR=/opt/mfem/mfem-4.5 -DCMAKE_BUILD_TYPE=Debug ..
    make -j6

Note: On Debian Linux, if you encounter `fatal error: Eigen/Dense: No such file or directory`, create symlinks:

    cd /usr/include
    sudo ln -sf eigen3/Eigen Eigen
    sudo ln -sf eigen3/unsupported unsupported

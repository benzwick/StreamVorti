https://github.com/spack/spack

Installation
To install spack, first make sure you have Python & Git. Then:

    git clone --depth=2 https://github.com/spack/spack.git

    # For bash/zsh/sh
    . spack/share/spack/setup-env.sh

Change to a compute node

    salloc -p work --time=3:00:00 --mem=32G -n16

Install MFEM and other dependencies

    spack install mfem@4.7 cgal@5.6 eigen@3.4.0

Clone StreamVorti into the correct directory

    git clone https://github.com/benzwick/StreamVorti.git

Change to a compute node and compile StreamVorti

    salloc -p work --time=3:00:00 --mem=32G -n16

    mkdir StreamVorti-build
    cd StreamVorti-build
    cmake \
      -DMFEM_DIR=/group/ems018/spack/spack/opt/spack/linux-haswell/mfem-4.7.0-htmjktlx5osy6mkbq2gesy4xizfgdjru/share/mfem/cmake \
      -DCMAKE_BUILD_TYPE=Debug ../StreamVorti

    spack load mfem@4.7 cgal@5.6 eigen@3.4.0
    cmake -DCMAKE_BUILD_TYPE=Debug ../StreamVorti

    make -j6

      -CGAL_DIR=/group/ems018/spack/spack/opt/spack/linux-haswell/cgal-5.6-6yu2sqozhwx5j47qpihtbt3n5mqvrduh/lib64/cmake/CGAL \

# Troubleshooting

1. Get MFEM install path
bash-4.4$ spack location -i mfem@4.7
/group/ems018/spack/spack/opt/spack/linux-haswell/mfem-4.7.0-htmjktlx5osy6mkbq2gesy4xizfgdjru

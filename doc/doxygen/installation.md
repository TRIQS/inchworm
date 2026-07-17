@page installation Installation

[TOC]

@note To guarantee reproducibility in scientific calculations we strongly recommend the use of
a stable [release](https://github.com/TRIQS/triqs/releases) of both TRIQS and its applications.

@section install_prereq Prerequisites

**inchworm** is built on top of the [TRIQS](https://triqs.github.io/triqs/) library; see the
[TRIQS installation instructions](https://triqs.github.io/triqs/latest/install.html). In the
following we assume that TRIQS is installed in the directory `path_to_triqs` and that its
environment has been sourced:

    source path_to_triqs/share/triqs/triqsvars.sh

The `inchworm` version must be compatible with your TRIQS library version; in particular the
major and minor version numbers must match.

@section install_steps Installation steps

1. Clone the `TRIQS/inchworm` repository from GitHub:

        git clone https://github.com/TRIQS/inchworm inchworm.src

2. Create and move to a build directory:

        mkdir inchworm.build && cd inchworm.build

3. Configure, build, test, and install:

        cmake ../inchworm.src -DCMAKE_INSTALL_PREFIX=path_to_inchworm
        make
        make test
        make install

@section install_options Custom CMake options

The compilation of `inchworm` can be configured using CMake options
(`cmake ../inchworm.src -DOPTION1=value1 -DOPTION2=value2 ...`):

| Option                                | Syntax                                     |
|---------------------------------------|--------------------------------------------|
| Specify an installation path          | `-DCMAKE_INSTALL_PREFIX=path_to_inchworm`  |
| Build in debugging mode               | `-DCMAKE_BUILD_TYPE=Debug`                 |
| Disable testing (not recommended)     | `-DBuild_Tests=OFF`                        |
| Build the documentation               | `-DBuild_Documentation=ON`                 |

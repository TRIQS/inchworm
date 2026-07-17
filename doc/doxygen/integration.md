@page integration Integration in C++ projects

[TOC]

**inchworm** installs a CMake package configuration, so it can be used from any CMake project
through `find_package`.

@section integration_cmake CMake

Make the installation visible to CMake — e.g. by adding the install prefix to
`CMAKE_PREFIX_PATH` — and in your `CMakeLists.txt` request the package and link against the
exported target:

    find_package(inchworm CONFIG REQUIRED)
    target_link_libraries(my_target inchworm::inchworm_c)

The imported target `inchworm::inchworm_c` carries the include directories and the transitive
TRIQS dependencies, so no further configuration is required.

@section integration_include Including the headers

The public headers are installed under `inchworm/`. For instance, to use the impurity solver:

    #include <inchworm/solver_core.hpp>

@section integration_env Environment

`inchworm` is built against TRIQS, so make sure the TRIQS environment is set up
(`source path_to_triqs/share/triqs/triqsvars.sh`) before configuring and running your project.

@page examples Examples

[TOC]

The primary user interface to **inchworm** is its Python `Solver`. The runnable example scripts
below are shipped in the `examples/` directory and are also rendered on the
[Python examples page](../examples.html) of the Sphinx documentation:

| Example               | Description                                     |
|-----------------------|-------------------------------------------------|
| `examples/onesite.py` | Single-site spinful Anderson impurity model     |
| `examples/twosite.py` | Two-site spinful Anderson impurity model        |

They can be run in parallel with, e.g.

    mpirun python onesite.py

For direct use of the C++ API, the @ref ex1 walkthrough follows the single-site solve, and
`test/c++/onesite.cpp` provides the corresponding self-contained C++ code.

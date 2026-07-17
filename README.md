[![build](https://github.com/TRIQS/inchworm/workflows/build/badge.svg)](https://github.com/TRIQS/inchworm/actions?query=workflow%3Abuild)

# inchworm — A TRIQS-based quantum impurity solver

> **⚠️ Warning:** This is a research code provided *as is*, with no guarantees, warranties, or support of any kind. It is intended for research purposes and may be incomplete, unstable, or subject to change without notice. Use at your own risk.

`inchworm` implements an equilibrium, imaginary-time Inchworm quantum Monte
Carlo solver for quantum impurity models. It is based on a continuous-time
hybridization expansion and is built on the [TRIQS](https://triqs.github.io/triqs/)
library.

## What is implemented

- Iterative calculation of the imaginary-time impurity propagator
  $\mathcal{U}(\tau)$ and subsequent sampling of the impurity Green's function
  $G(\tau)$.
- General local interactions expressed as TRIQS many-body operators and
  block-structured, matrix-valued hybridization functions.
- Inclusion–exclusion summation of the proper hybridization diagrams, with a
  direct enumeration implementation available for validation.
- MPI-parallel Monte Carlo sampling with automatic calibration, perturbation
  order and sign diagnostics, and optional order-resolved measurements.
- Optional subtraction of a fitted discrete bath before sampling and HDF5
  serialization of solver inputs and results.
- C++ and Python interfaces.

The present implementation is an experimental equilibrium solver; it does not
implement real-time or nonequilibrium Inchworm calculations. The API is not
stable, and users must check convergence with respect to Monte Carlo statistics,
the imaginary-time meshes, and any perturbation-order cutoff.

## Expansion around a fitted discrete bath

By default, the solver expands around the atomic Hamiltonian in the full
hybridization function $\Delta$. The Python interface also supports a shifted
expansion point, enabled by setting the solve parameter `n_bath_sites_ED` to the
number of auxiliary bath sites per hybridization block. It fits a discrete-bath
hybridization $\Delta_{\mathrm{ED}}$, adds the corresponding bath orbitals,
energies, and couplings to the Hamiltonian treated by exact diagonalization, and
samples the residual hybridization

$$
\widetilde{\Delta} = \Delta - \Delta_{\mathrm{ED}}.
$$

This is an add-and-subtract reorganization of the expansion, rather than an ED
approximation to the final answer: when the residual expansion is converged, it
recovers the problem defined by the original $\Delta$. A good bath fit can
reduce the required perturbation order, at the cost of an exponentially larger
local Hilbert space. The default `n_bath_sites_ED=0` selects the ordinary atomic
expansion.

This option was introduced in the original
[August 2021 implementation commit](https://github.com/TRIQS/inchworm/commit/c21ef1e12d82d33a20bbbf3fca77c4a13d8b5a67).

## Build and run

`inchworm` requires a compatible installation of TRIQS 4.0. After loading the
TRIQS environment, configure, build, and test with:

```bash
cmake -S . -B build -GNinja -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build -j 16
ctest --test-dir build -j 16 --output-on-failure
```

The [single-site example](examples/onesite.py) is the smallest complete solver
run. It constructs a spinful Anderson impurity coupled to a discrete bath,
computes $\mathcal{U}(\tau)$ and $G(\tau)$, and writes them to
`onesite.out.h5`:

```bash
PYTHONPATH="$PWD/build/python:$PYTHONPATH" python examples/onesite.py
```

The [two-site example](examples/twosite.py) demonstrates a two-orbital impurity
with matrix-valued hybridization blocks and writes `twosite.out.h5`:

```bash
PYTHONPATH="$PWD/build/python:$PYTHONPATH" python examples/twosite.py
```

These calculations use research-scale Monte Carlo parameters and may take some
time. The corresponding [single-site](test/python/onesite.py) and
[two-site](test/python/twosite.py) regression tests compare against committed
reference results.

## Algorithm and implementation notes

The algorithm is based on the original Inchworm construction, its equilibrium
multiorbital formulation, and the inclusion–exclusion method for summing proper
diagrams:

1. G. Cohen, E. Gull, D. R. Reichman, and A. J. Millis, “Taming the Dynamical
   Sign Problem in Real-Time Evolution of Quantum Many-Body Problems,”
   [Phys. Rev. Lett. **115**, 266802 (2015)](https://doi.org/10.1103/PhysRevLett.115.266802).
2. E. Eidelstein, E. Gull, and G. Cohen, “Multiorbital Quantum Impurity Solver
   for General Interactions and Hybridizations,”
   [Phys. Rev. Lett. **124**, 206405 (2020)](https://doi.org/10.1103/PhysRevLett.124.206405).
3. A. Boag, E. Gull, and G. Cohen, “Inclusion-exclusion principle for many-body
   diagrammatics,”
   [Phys. Rev. B **98**, 115152 (2018)](https://doi.org/10.1103/PhysRevB.98.115152).

The repository's [implementation notes](doc/notes/inchworm_notes.tex) derive the
imaginary-time expansion and document details specific to this code. In
particular, the implementation uses dressed propagators on both sides of the
split time, a corresponding definition of diagram connectivity, and an
additional permutation-sign factor in the adapted inclusion–exclusion formula.

## License

`inchworm` is distributed under the
[GNU General Public License, version 3 or later](LICENSE).

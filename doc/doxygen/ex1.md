@page ex1 Single-site Anderson impurity model

[TOC]

This example solves a single-site, spinful Anderson impurity model with the Inchworm QMC
solver. It corresponds to `test/c++/onesite.cpp` (C++) and `examples/onesite.py` (Python).

@section ex1_setup Setup

Construct the solver from a set of construction parameters — the inverse temperature, the block
structure of the Green's functions, and the imaginary-time mesh sizes:

    inchworm::constr_params_t cp;
    cp.beta        = 2.0;
    cp.gf_struct   = {{"up", 1}, {"dn", 1}};
    cp.n_tau_green = 5;
    cp.n_tau_inch  = 21;
    cp.n_tau       = 10001;

    inchworm::solver_core S{cp};

Set the imaginary-time hybridization function `S.Delta_tau` and define the local interaction
Hamiltonian (here a Hubbard @f$U@f$ on the impurity site) in the solve parameters.

@section ex1_solve Solving

The solver exposes several strategies, all taking the solve parameters `sp`:

    S.solve_cthyb(sp, tau_max);                        // bare hybridization expansion
    S.solve_self_consistently(sp, u_tau, tau_split, tau_max);
    S.solve_inchworm(sp);                              // full inchworm propagation
    S.solve_green(sp);                                 // sample G(tau) from U(tau)

After solving, the impurity propagator is available as `S.u_tau` and the impurity Green's
function as `S.G_tau`.

@section ex1_ref Full source

See `test/c++/onesite.cpp` for the complete, compilable version and `examples/onesite.py` for
the equivalent Python script.

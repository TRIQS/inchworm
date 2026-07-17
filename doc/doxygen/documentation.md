@page documentation API Documentation

[TOC]

This is the detailed API documentation of the **inchworm** C++ source code, generated from the
comments in the headers under `c++/inchworm/`.

The central class is inchworm::solver_core, which drives the Inchworm Monte Carlo solve and
exposes the impurity propagator and Green's function. Its construction and solve parameters are
described by inchworm::constr_params_t and inchworm::solve_params_t, and results are collected in
inchworm::qmc_results_t.

Use the navigation bar above to browse the full list of classes and files.

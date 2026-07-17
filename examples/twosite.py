# Copyright (c) 2019--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""Solve a two-site spinful Anderson impurity model with Inchworm QMC."""

from itertools import product

from h5 import HDFArchive
from inchworm import Solver
import numpy as np
from triqs.operators import n
from triqs.utility import mpi


def one_fermion(tau, eps, beta):
    """Imaginary-time Green's function of one fermionic bath level."""
    if eps >= 0:
        return -np.exp(-tau * eps) / (1.0 + np.exp(-beta * eps))
    return -np.exp((beta - tau) * eps) / (1.0 + np.exp(beta * eps))


beta = 2.0
solver = Solver(
    beta=beta,
    gf_struct=[("up", 2), ("dn", 2)],
    n_tau_green=5,
    n_tau_inch=21,
    n_tau=10001,
)

theta = np.array([[0.1, 0.3, 0.4], [0.1, 0.2, 0.4]])
epsilon = np.array([1.0, -1.0, 1.2])
n_sites, n_bath_sites = theta.shape

for block in ("up", "dn"):
    solver.Delta_tau[block].data[:] = 0
    tau = np.array([point.value for point in solver.Delta_tau[block].mesh])
    for i, j, bath_site in product(
        range(n_sites), range(n_sites), range(n_bath_sites)
    ):
        solver.Delta_tau[block][i, j].data[:] += (
            theta[i, bath_site]
            * theta[j, bath_site]
            * one_fermion(tau, eps=epsilon[bath_site], beta=beta)
        )

U = 1.0
h_imp = U * (
    n("up", 0) * n("dn", 0) + n("up", 1) * n("dn", 1)
) - U / 2 * (n("up", 0) + n("dn", 0) + n("up", 1) + n("dn", 1))
solver.solve(h_imp=h_imp, n_cycles=40000)

if mpi.is_master_node():
    with HDFArchive("twosite.out.h5", "w") as archive:
        archive["u_tau"] = solver.u_tau
        archive["G_tau"] = solver.G_tau

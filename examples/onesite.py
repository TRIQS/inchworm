#!/usr/bin/env python
###############################################################################
#
# inchworm: A TRIQS based impurity solver
#
# Copyright (c) 2019 The Simons foundation
#   authors: Nils Wentzell
#
# inchworm is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# inchworm is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# inchworm. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

"""Solve a single-site spinful Anderson impurity model with Inchworm QMC."""

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
    gf_struct=[("up", 1), ("dn", 1)],
    n_tau_green=5,
    n_tau_inch=21,
    n_tau=10001,
)

theta = np.array([0.1, 0.3, 0.4])
epsilon = np.array([1.0, -1.0, 1.2])

for block in ("up", "dn"):
    solver.Delta_tau[block].data[:] = 0
    tau = np.array([point.value for point in solver.Delta_tau[block].mesh])
    for bath_site in range(len(epsilon)):
        solver.Delta_tau[block][0, 0].data[:] += (
            theta[bath_site] ** 2
            * one_fermion(tau, eps=epsilon[bath_site], beta=beta)
        )

U = 1.0
h_imp = U * n("up", 0) * n("dn", 0) - U / 2 * (
    n("up", 0) + n("dn", 0)
)
solver.solve(h_imp=h_imp, n_cycles=50000)

if mpi.is_master_node():
    with HDFArchive("onesite.out.h5", "w") as archive:
        archive["u_tau"] = solver.u_tau
        archive["G_tau"] = solver.G_tau

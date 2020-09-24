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
#!/usr/bin/env python

import unittest
import numpy as np
from itertools import product

from inchworm import Solver

from triqs.gf import *
from h5 import *
from triqs.operators import *
from triqs.utility.h5diff import h5diff

def one_fermion(tau, eps, beta):
    import math
    if (eps >= 0):
        return -math.exp(-tau * eps) / (1. + math.exp(-beta * eps))
    else:
        return -math.exp((beta - tau) * eps) / (1. + math.exp(beta * eps))

class test_hubbard(unittest.TestCase):

    # System Parameters
    U = 1.0

    # # Construct Parameters
    cp = {}
    cp["beta"]        = 1.0
    cp["gf_struct"]   = [("up", [0]), ("dn", [0])]
    cp["n_tau_green"] = 5
    cp["n_tau_inch"]  = 21
    cp["n_tau"]       = 10001

    theta   = np.array([0.1, 0.3, 0.4])
    epsilon = np.array([1.0, -1.0, 1.2])

    n_bath = len(epsilon)

    # Set up the Solver
    S = Solver(**cp)
    for bl in ["up", "dn"]:
        S.Delta_tau[bl].data[:] = 0
        for mp in S.Delta_tau["up"].mesh:
            for bl, k in product(["up", "dn"], range(n_bath)):
                S.Delta_tau[bl][mp][0, 0] = S.Delta_tau[bl][mp][0, 0] \
                  + theta[k] * theta[k] * one_fermion(tau=mp.value, eps=epsilon[k], beta=cp["beta"])

    # Solve Parameters
    sp = {}
    sp["h_imp"] =  U   *  n("up",0) * n("dn",0) \
                 - U/2 * (n("up",0) + n("dn",0))
    sp["n_cycles"] = 100000
    sp["length_cycle"] = 10
    sp["n_warmup_cycles"] = 20
   
    # Solve the model
    S.solve(**sp)

    # Store the Result
    with HDFArchive("onesite.out.h5", 'w') as arch:
        arch["u_tau"] = S.u_tau
        arch["G_tau"] = S.G_tau

    # -------- Compare ---------
    h5diff("onesite.out.h5", "onesite.ref.h5")

if __name__ == '__main__':
    unittest.main()

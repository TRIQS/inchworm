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
    cp["beta"]        = 2.0
    cp["gf_struct"]   = [("up", 2), ("dn", 2)]
    cp["n_tau_green"] = 5
    cp["n_tau_inch"]  = 21
    cp["n_tau"]       = 10001

    theta   = np.array([[0.1, 0.3, 0.4], [0.1, 0.2, 0.4]])
    epsilon = np.array([1.0, -1.0, 1.2])

    n_site = len(theta)
    n_bath = len(epsilon)

    # Set up the Solver
    S = Solver(**cp)
    for bl in ["up", "dn"]:
        S.Delta_tau[bl].data[:] = 0
        for mp in S.Delta_tau["up"].mesh:
            for i, j, k in product(range(n_site), range(n_site), range(n_bath)):
                S.Delta_tau[bl][mp][i, j] = S.Delta_tau[bl][mp][i, j] \
                  + theta[i, k] * theta[j, k] * one_fermion(tau=mp.value, eps=epsilon[k], beta=cp["beta"])

    # Solve Parameters
    sp = {}
    sp["h_imp"] =  U   * (n("up",0) * n("dn",0) + n("up",1) * n("dn",1)) \
                 - U/2 * (n("up",0) + n("dn",0) + n("up",1) + n("dn",1))
    sp["n_cycles"] = 40000
   
    # Solve the model
    S.solve(**sp)

    # Store the Result
    with HDFArchive("twosite.out.h5", 'w') as arch:
        arch["u_tau"] = S.u_tau
        arch["G_tau"] = S.G_tau

    # -------- Compare ---------
    h5diff("twosite.out.h5", "twosite.ref.h5")

if __name__ == '__main__':
    unittest.main()

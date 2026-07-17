# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

#!/usr/bin/env python

import unittest
import numpy as np
from itertools import product

from inchworm import Solver

from triqs.gfs import *
from h5 import *
from triqs.operators import *
from triqs.utility.h5diff import h5diff

def one_fermion(tau, eps, beta):
    if (eps >= 0):
        return -np.exp(-tau * eps) / (1. + np.exp(-beta * eps))
    else:
        return -np.exp((beta - tau) * eps) / (1. + np.exp(beta * eps))

class test_hubbard(unittest.TestCase):

    # System Parameters
    U = 1.0

    # # Construct Parameters
    cp = {}
    cp["beta"]        = 2.0
    cp["gf_struct"]   = [("up", 1), ("dn", 1)]
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
        taulst = np.array([mp.value for mp in S.Delta_tau[bl].mesh])
        for bl, k in product(["up", "dn"], range(n_bath)):
            S.Delta_tau[bl][0, 0].data[:] = S.Delta_tau[bl][0, 0].data[:] \
              + theta[k] * theta[k] * one_fermion(taulst, eps=epsilon[k], beta=cp["beta"])

    # Solve Parameters
    sp = {}
    sp["h_imp"] =  U   *  n("up",0) * n("dn",0) \
                 - U/2 * (n("up",0) + n("dn",0))
    sp["n_cycles"] = 50000
   
    # Solve the model
    S.solve(**sp)

    # Store the Result
    with HDFArchive("onesite.out.h5", 'w') as arch:
        arch["u_tau"] = S.u_tau
        arch["G_tau"] = S.G_tau

    # -------- Compare ---------
    h5diff("onesite.out.h5", "onesite.ref.h5")

if __name__ == '__main__':
    # exit=False: the h5diff comparison in the class body is the actual gate
    # (it raises on mismatch); unittest.main would otherwise exit 5 on Python
    # 3.12 as the class defines no test_* methods.
    unittest.main(exit=False)

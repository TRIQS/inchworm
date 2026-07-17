# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import sympy as sy
import numpy as np
from sympy.functions import cosh

theta = sy.Symbol('theta')
beta = sy.Symbol('beta')
f = (cosh(beta*theta/2))**4


H=sy.Matrix([[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]])
x=H.charpoly(sy.Symbol('lamda'))
print(sy.factor(x))

zz=H.eigenvects()
sy.pprint(zz)


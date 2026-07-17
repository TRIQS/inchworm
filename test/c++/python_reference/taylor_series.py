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

def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

def taylor(function,theta0,N):
  n = 0
  serie = 0
  while n <= N:
    serie = (function.diff(theta,n).subs(theta,theta0))/(factorial(n))*(theta-theta0)**n + serie
    n += 1
  return serie

x = taylor(f,0,16)

print x
print
sy.pprint(x)
print
print sy.latex(x)

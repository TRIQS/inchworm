# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import sympy as sy
import numpy as np

import re

def propagator_product_for_U0():

 beta = 2.
 tau_split = 0.5*beta
 theta = 3.2
 line = "c(0.052559432586)_0 c^+(1.170881322858)_0 c(1.174520496552)_0 c^+(1.324650712620)_0 c(1.522218051255)_0 c^+(1.564001699808)_0 c(1.914343198276)_0 c^+(1.957088277200)_0"
 size = 14
 
 
 pos = [m.start() for m in re.finditer('\(', line)]
 #[2,17,30,45,58]
 
 times = [0.0]
 print 'times'
 for i in range(len(pos)):
  tau = float(line[pos[i]+1:pos[i]+1+size])
  if(tau_split < tau and tau_split> times[-1]):
    times.append( tau_split )
  times.append( tau )
  #print line[pos[i]+1:pos[i]+1+size]  
 times.append(beta)
 print times
 print
 
 product = 1.0
 print 'product'
 for i in range(len(times)-1):
   cosh = np.cosh( (times[i+1]-times[i]) * theta / 2. )
   product *= cosh*cosh
   print cosh*cosh
 
 print
 print product

propagator_product_for_U0()
exit()






from sympy.functions import cosh

theta = sy.Symbol('theta')
beta  = sy.Symbol('beta')
tau_s = sy.Symbol('tau_s')

B   = sy.Symbol('B')
t_s = sy.Symbol('t_s')

f = (cosh(B-t_s) * cosh(t_s))**2
#f = (cosh(t_s))**2

def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

def taylor(function,t_s0,N):
  return (function.diff(t_s,N).subs(t_s,t_s0))/(factorial(N))*(t_s-t_s0)**N


for n in range(8):

  x = sy.simplify(taylor(f,0,n))
  
  print
  print 'order %d' % n
  print  
  #print x
  #print
  sy.pprint(x)
  print
  print
  print
  print
  print
  #sy.pprint(x.subs(epsilon,0))
  #print
  #print sy.latex(x)

exit()



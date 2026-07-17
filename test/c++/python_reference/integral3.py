# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import sympy as sy
import numpy as np
from sympy.functions import cosh, exp, sinh

epsilon = sy.Symbol('epsilon')
theta = sy.Symbol('theta')
beta = sy.Symbol('beta')
tau1 = sy.Symbol('tau_1')
tau2 = sy.Symbol('tau_2')

taus = sy.Symbol('tau_s')
B = sy.Symbol('Dt')

t = sy.Symbol('t')
ts = 0#sy.Symbol('t_s')

t_1  = sy.Symbol('t_1')
t_2  = sy.Symbol('t_2')
t_3  = sy.Symbol('t_3')
t_4  = sy.Symbol('t_4')
t_5  = sy.Symbol('t_5')


#f = (cosh(beta*theta/2))**2


#f1 = ((exp(ts-t)   + exp(-(ts-t)))   * (exp(t)    + exp(-t) )) **2
#f2 = ((exp(beta-t) + exp(-(beta-t))) * (exp(t-ts) + exp(-t+ts) )) **2

#f1 = (cosh(ts-t_1p) * cosh(t_1p-t_1) * cosh(t_1)) **2
f1 = (cosh(t_2-t_1) * cosh(t_1-ts)) **2
f2 = (cosh(t_3-t_2) * cosh(t_2-t_1) * cosh(t_1-ts)) **2
f3 = (cosh(t_4-t_3) * cosh(t_3-t_2) * cosh(t_2-t_1) * cosh(t_1-ts)) **2
f4 = (cosh(t_5-t_4) * cosh(t_4-t_3) * cosh(t_3-t_2) * cosh(t_2-t_1) * cosh(t_1-ts)) **2
f5 = (cosh(B  -t_5) * cosh(t_5-t_4) * cosh(t_4-t_3) * cosh(t_3-t_2) * cosh(t_2-t_1) * cosh(t_1-ts)) **2

#y1 = sy.simplify(sy.integrate(                      (cosh(t_2-t_1) **2) , (t_1, ts ,t_2) ))
#print 
#print 
#sy.pprint(y1)

#y1 = sy.simplify(sy.integrate( cosh(2.*t_1-2.*ts) * (cosh(t_2-t_1) **2) , (t_1, ts ,t_2) ))
#print 
#print 
#sy.pprint(y1)

#y1 = sy.simplify(sy.integrate( sinh(2.*t_1-2.*ts) * (cosh(t_2-t_1) **2) , (t_1, ts ,t_2) ))
#print 
#print 
#sy.pprint(y1)

#y1 = sy.simplify(sy.integrate(  (t_1-ts) * cosh(2.*t_1-2.*ts) * (cosh(t_2-t_1) **2) , (t_1, ts ,t_2) ))
#print 
#print 
#sy.pprint(y1)


#exit()


#p1 = sy.simplify( sy.integrate( sy.integrate(  f1 , (t_1, 0 ,t_1p) ), (t_1p, 0 , ts)   ))
y1 = sy.simplify(sy.integrate( f1 , (t_1, ts, t_2) ))
print 
sy.pprint(y1)
print 
print y1

y2 = sy.simplify(sy.integrate( cosh(t_3-t_2)**2 * y1 , (t_2, ts, t_3) ))
print 
sy.pprint(y2)
print 
print y2

y3 = sy.simplify(sy.integrate( cosh(t_4-t_3)**2 * y2 , (t_3, ts, t_4) ))
print 
sy.pprint(y3)
print 
print y3

y4 = sy.simplify(sy.integrate( cosh(t_5-t_4)**2 * y3 , (t_4, ts, t_5) ))
print 
sy.pprint(y4)
print 
print y4

y5 = sy.simplify(sy.integrate( cosh(B  -t_5)**2 * y4 , (t_5, ts ,  B) ))
print 
sy.pprint(y5)
print 
print y5


exit()



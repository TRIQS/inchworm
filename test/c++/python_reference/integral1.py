
import sympy as sy
import numpy as np
from sympy.functions import cosh, exp

epsilon = sy.Symbol('epsilon')
theta = sy.Symbol('theta')
beta = sy.Symbol('beta')
tau1 = sy.Symbol('tau_1')
tau2 = sy.Symbol('tau_2')

taus = sy.Symbol('tau_s')
B = sy.Symbol('B')

t = sy.Symbol('t')
ts = sy.Symbol('t_s')

#f = (cosh(beta*theta/2))**2


#f1 = ((exp(ts-t)   + exp(-(ts-t)))   * (exp(t)    + exp(-t) )) **2
#f2 = ((exp(beta-t) + exp(-(beta-t))) * (exp(t-ts) + exp(-t+ts) )) **2

f1 = (cosh(ts-t)   * cosh(t)) **2
f2 = (cosh(B-t) * cosh(t-ts)) **2


sy.pprint(f1)
sy.pprint(f2)
print 

p1 = sy.simplify( sy.integrate(  f1 , (t,0 ,ts) ))
p2 = sy.simplify( sy.integrate(  f2 , (t,ts,B) ))

print
sy.pprint(p1)
print
sy.pprint(p2)

print 
print 
print 
print 
print p1
print 
print p2
#sy.pprint(sy.simplify(p1*p2))


exit()


from sympy import *
x, y = symbols('x, y')
integrate(x + y)

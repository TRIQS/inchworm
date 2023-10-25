
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

t_1  = sy.Symbol('t_1')
t_1p = sy.Symbol('t_1p')
t_2  = sy.Symbol('t_2')
t_2p = sy.Symbol('t_2p')


#f = (cosh(beta*theta/2))**2


#f1 = ((exp(ts-t)   + exp(-(ts-t)))   * (exp(t)    + exp(-t) )) **2
#f2 = ((exp(beta-t) + exp(-(beta-t))) * (exp(t-ts) + exp(-t+ts) )) **2

f1 = (cosh(ts-t_1p) * cosh(t_1p-t_1) * cosh(t_1)) **2
f2 = (cosh(B -t_2p) * cosh(t_2p-t_2) * cosh(t_2-ts)) **2

thet= 0.8
bet = 2.0
ta_s= 0.5*bet

t1  = 0.67376 * thet / 2.
t2  = 1.36633 * thet / 2.
ts1 = 0.5 * bet * thet / 2.
b   = bet * thet / 2.

print 1.0 * (f1.subs([(t, t1), (ts, ts1)])) * (f2.subs([(t, t2), (ts, ts1), (B,b)]))
print
print thet*thet/2.

sy.pprint(f1)
sy.pprint(f2)
#print 
#a1 = 0.08928 / 0.08436757
#a2 = 0.52891 / 0.39499137
#a3 = 8.92976 / 2.88053490
#print a1, a2, a3, a2/a1, a3/a2
#exit()

#p1 = sy.simplify( sy.integrate( sy.integrate(  f1 , (t_1, 0 ,t_1p) ), (t_1p, 0 , ts)   ))
p2 = sy.simplify( sy.integrate( sy.integrate(  f2 , (t_2, ts ,t_2p) ), (t_2p, ts , B)   ))
#p2 = sy.simplify( sy.integrate(  f2 , (t,ts,B) ))

print
#sy.pprint(p1)
print
sy.pprint(p2)

print 
print 
print 
print 
#print p1
print 
print p2
#sy.pprint(sy.simplify(p1*p2))


exit()



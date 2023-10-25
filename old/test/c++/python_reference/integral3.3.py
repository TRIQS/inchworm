
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

t_0  = sy.Symbol('t_0')
t_1  = sy.Symbol('t_1')
t_2  = sy.Symbol('t_2')
t_e  = sy.Symbol('t_e')

a  = sy.Symbol('a')


#f = (cosh(beta*theta/2))**2


#f1 = ((exp(ts-t)   + exp(-(ts-t)))   * (exp(t)    + exp(-t) )) **2
#f2 = ((exp(beta-t) + exp(-(beta-t))) * (exp(t-ts) + exp(-t+ts) )) **2

#f1 = (cosh(ts-t_1p) * cosh(t_1p-t_1) * cosh(t_1)) **2
f1 = sy.simplify( exp(a*t_1) * ( ((exp(t_2-t_1) + exp(t_1-t_2))/2) **2) * (((exp(t_1-t_0) + exp(t_0-t_1))/2) **2)  )

print 
print 
sy.pprint(f1)
print 
print 


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


a#p1 = sy.simplify( sy.integrate( sy.integrate(  f1 , (t_1, 0 ,t_1p) ), (t_1p, 0 , ts)   ))
y1 = sy.simplify(sy.integrate( f1 , (t_1, t_0, t_2) ))
print 
sy.pprint(y1)
print 
print y1


exit()



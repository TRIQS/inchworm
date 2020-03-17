
import sympy as sy
import numpy as np
from sympy.functions import cosh, exp

epsilon = sy.Symbol('epsilon')
theta = sy.Symbol('theta')
beta = sy.Symbol('beta')
tau1 = sy.Symbol('tau_1')
tau2 = sy.Symbol('tau_2')
f = (cosh(beta*theta/2))**2

U = sy.Function('U')(beta)

#sy.pprint(sy.integrate(x/(x**2+2*x+1), x))
print U

sy.pprint(sy.integrate(cosh((tau1)*theta/2)**2, tau1))

#sy.pprint(sy.integrate( sy.integrate( (cosh((beta-tau1)*theta/2))**2  * (cosh((tau1-tau2)*theta/2))**2  * (cosh((tau2)*theta/2))**2   , tau1)  , tau2 ))



exit()


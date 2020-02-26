
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


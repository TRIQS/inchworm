
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

x = taylor(f,0,10)

print x
print
sy.pprint(x)
print
print sy.latex(x)

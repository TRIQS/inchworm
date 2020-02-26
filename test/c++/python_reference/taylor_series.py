
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
    serie = serie + (function.diff(theta,n).subs(theta,theta0))/(factorial(n))*(theta-theta0)**n
    n += 1
  return serie
    
print taylor(f,0,10)


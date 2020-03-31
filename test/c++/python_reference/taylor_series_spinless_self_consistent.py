
import sympy as sy
import numpy as np
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



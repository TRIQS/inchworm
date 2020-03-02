
import sympy as sy
import numpy as np
from sympy.functions import cosh, exp

epsilon = sy.Symbol('epsilon')
theta = sy.Symbol('theta')
beta = sy.Symbol('beta')
f = (cosh(beta*theta/2))**2

def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

def taylor(function,theta0,N):
  return (function.diff(theta,N).subs(theta,theta0))/(factorial(N))*(theta-theta0)**N
  



H=sy.Matrix([[0    ,theta  ],
             [theta,epsilon]])

'''
[T, D] = H.diagonalize(normalize=True)
print
sy.pprint(sy.simplify(T))
T2 = sy.simplify(T)
print
sy.pprint(D)


norm1 = (T2*T2.transpose())
norm2 = (T2.transpose()*T2)

print
sy.pprint(sy.simplify(norm1))
print
sy.pprint(sy.simplify(norm2))

exit()

H2 = sy.simplify(norm1*T*D*T.transpose())
D2 = sy.simplify(norm2*T.transpose()*H*T)

I1 = sy.simplify(T*T.transpose())
I2 = sy.simplify(T.transpose()*T)

print
sy.pprint(H2)
print
sy.pprint(D2)
print
sy.pprint(I1)
print 
sy.pprint(I2)
print 


exit()
'''


epsilon = 0

H=sy.Matrix([[0,0,0,0],
             [0,0,theta,0],
             [0,theta,epsilon,0],
             [0,0,0,epsilon]])

[T, D] = H.diagonalize(normalize=True)
zz=H.eigenvects()
sy.pprint(zz)

#sy.pprint(H)
#sy.pprint(T)
#sy.pprint(D)

#print dir(T)
norm1 = (T*T.transpose())**(-1)
norm2 = (T.transpose()*T)**(-1)

H2 = sy.simplify(norm1*T*D*T.transpose())
D2 = sy.simplify(norm2*T.transpose()*H*T)

sy.pprint(H2)
sy.pprint(D2)


Ud=sy.Matrix([[0,0,0,0],
              [0,0,0,0],
              [0,0,0,0],
              [0,0,0,0]])

for ii in range(4):
  Ud[ii,ii] += exp(-beta*D[ii,ii])


sy.pprint(Ud)
#exit()

U = sy.simplify(norm1* T*Ud*T.transpose())

U_traced = sy.Matrix([[0,0],
                      [0,0]])

for ii in range(2):
  for jj in range(2):
    for kk in range(2):
      U_traced[ii,jj] += U[ii+2*kk,jj+2*kk]


sy.pprint(U)
sy.pprint(U_traced)


zb = 1 + exp(-beta*epsilon)

print '\n'
sy.pprint(zb)

print '\n'
sy.pprint(U_traced[0,0]/zb)
print '\n'
sy.pprint(U_traced[1,1]/zb)

for n in range(4):

  x = sy.simplify(taylor(U_traced[0,0]/zb,0,2*n))
  
  print
  print
  print 'order %d' % n
  print  
  print x
  print
  sy.pprint(x)
  print
  #sy.pprint(x.subs(epsilon,0))
  print
  print sy.latex(x)

exit()





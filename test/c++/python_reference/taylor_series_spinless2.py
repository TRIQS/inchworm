# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import sympy as sy
import numpy as np
from sympy.functions import cosh, exp


t1 = 0.2
e1 = -10.5
b1 = 2.0
ts = b1*0.95

B = e1*b1

eB = np.exp(B)
emB= np.exp(-B)

es = np.exp(ts*e1)
ems= np.exp(-ts*e1)

c_ = t1*t1/(e1*e1*(1+eB))
cm = t1*t1/(e1*e1*(1+emB))

print eB, emB, es, ems, c_, cm
print
U0t0 = (1-c_ * (  2-eB*ems + B - es ) )
U0t1 = (1+cm * ( -2+emB*es + B + ems ) )
print U0t0
print U0t1
print 
U1t0 = (-c_ * (eB - es) * (ems-1)) 
U1t1 = (-cm * (emB - ems) * (es-1) )
print U1t1
print U1t1
print 
#U1U1t0 = c_*c_*(1-eB*ems + B-ts*e1) * (1-es + ts*e1) 
#U1U1t1 = cm*cm*(1-emB*es + B-ts*e1) * (1-ems + ts*e1) 
#print U1U1t0
#print U1U1t1
Utot0 = U0t0 + U1t0
Utot1 = U0t1 + U1t1
print
print Utot0
print Utot1

print 1-c_ * (1. - eB + B)
print 1-cm * (1. - emB - B)

print
print (1-c_ * (1. - eB*ems + B-ts*e1)) * (1-c_ * (1. - es  + ts*e1))
print (1-cm * (1. - emB*es - B+ts*e1)) * (1-cm * (1. - ems - ts*e1))



c_Bs = t1*t1/(e1*e1*(1+eB*ems))
cmBs = t1*t1/(e1*e1*(1+emB*es))

c_s = t1*t1/(e1*e1*(1+es))
cms = t1*t1/(e1*e1*(1+ems))

print
print (1-c_Bs * (1. - eB*ems + B-ts*e1)) * (1-c_s * (1. - es  + ts*e1))
print (1-cmBs * (1. - emB*es - B+ts*e1)) * (1-cms * (1. - ems - ts*e1))


#exit()




epsilon = sy.Symbol('epsilon')
theta = sy.Symbol('theta')
beta = sy.Symbol('beta')

tau_m = sy.Symbol('tau_m')


e = sy.Symbol('e')
b = sy.Symbol('b')
f = (cosh(beta*theta/2))**2

def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

def taylor(function,theta0,N):
  return (function.diff(theta,N).subs(theta,theta0))/(factorial(N))*(theta-theta0)**N
  


'''
H=np.array([[0, 1],
            [1,10 ]])


[T, D] = np.linalg.eigh(H)
print T
print
print D

exit()
'''



H=sy.Matrix([[ 0,theta ],
             [ theta,epsilon ]])


[T, D] = H.diagonalize(normalize=True)
print
sy.pprint(sy.simplify(T))
T2 = sy.simplify(T)
print
sy.pprint(D)


H2 = sy.simplify(T*D*T.transpose())
D2 = sy.simplify(T.transpose()*H*T)

#I1 = sy.simplify(T*T.transpose())
#I2 = sy.simplify(T.transpose()*T)



print
sy.pprint(H)
print '='
sy.pprint(H2)
print

print
sy.pprint(D)
print '='
sy.pprint(D2)
print
#sy.pprint(I1)
#print 
#sy.pprint(I2)
print 
print 
print 
print 
print 
print 
print 



#epsilon = 0

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
  Ud[ii,ii] += exp(-tau_m*D[ii,ii])


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



for ii in range(2):
 for n in range(4):

  x = sy.simplify(taylor(U_traced[ii,ii]/zb,0,n))
  
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
  #print sy.latex(x)
  print x.subs([(epsilon,e1),(theta,t1),(beta, b1),(tau_m, ts)])
  print 

for ii in range(2):
  print '\nU_traced[%d]\n' %ii
  print (U_traced[ii,ii]/zb).subs([(epsilon,e1),(theta,t1),(beta, b1),(tau_m, b1)])
  print 
  print (U_traced[ii,ii]/zb).subs([(epsilon,e1),(theta,t1),(beta, ts),(tau_m, ts)]) * (U_traced[ii,ii]/zb).subs([(epsilon,e1),(theta,t1),(beta, b1 - ts),(tau_m, b1 - ts)])
  print 
  print (U_traced[ii,ii]/zb).subs([(epsilon,e1),(theta,t1),(beta, b1),(tau_m, ts)]) * (U_traced[ii,ii]/zb).subs([(epsilon,e1),(theta,t1),(beta, b1),(tau_m, b1-ts)])
  
exit()





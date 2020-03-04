import numpy as np


def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

#print factorial(4)

theta = 0.02
U = 0.
beta = 2.

tau1 =1.90054
tau1p=0.21846

print np.exp((beta-tau1+tau1p)*U/2.)
print np.exp((tau1-tau1p)*U/2.)

print
print 1.*(np.cosh(beta*theta/2.))
print 1.*(np.cosh(beta*theta/2.))**2
print 1.*(np.cosh(beta*theta/2.))**4
print 1.*(np.cosh(beta*theta/2.))**8
# c^+(0.21846)_0 c(1.90054)_0
print
print "spin"
print 1.
print (beta*theta)**2
print 5./12.*(beta*theta)**4
print 17./180.*(beta*theta)**6
print 13./1008.*(beta*theta)**8
print 1. + (beta*theta)**2 + 5./12.*(beta*theta)**4 + 17./180.*(beta*theta)**6 + 13./1008.*(beta*theta)**8

print
print "spinless"
print 1.
print 1./2.*(beta*theta)**2
print 1./12.*(beta*theta)**4
print 1./180.*(beta*theta)**6
print 1./5040.*(beta*theta)**8
print 1. + 1./2.*(beta*theta)**2  +1./12.*(beta*theta)**4   +1./180.*(beta*theta)**6  +1./5040.*(beta*theta)**8


import numpy as np


def factorial(n):
  if n <= 0:
    return 1
  else:
    return n*factorial(n-1)

#print factorial(4)

theta = 0.02
U = 10.
beta = 2.

tau1 =1.90054
tau1p=0.21846

print np.exp((beta-tau1+tau1p)*U/2.)
print np.exp((tau1-tau1p)*U/2.)

# c^+(0.21846)_0 c(1.90054)_0


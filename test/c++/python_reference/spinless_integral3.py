# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import numpy as np

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
fig = plt.figure()
ax = plt.axes()

plt.ylabel('U-U_{approx}')
plt.xlabel('theta')

plt.yscale("log")
plt.xscale("log")





  
def y1(Dt):
  s2 = np.sinh(2 * Dt)
  c2 = np.cosh(2 * Dt)
  return 1. / 8 * Dt * (c2 + 2.) + 5. / 16 * s2

def y2(Dt):
  s3 = np.sinh(2 * Dt)
  c3 = np.cosh(2 * Dt)
  return 1. / 64 * Dt * Dt * (c3 + 4.) + 15. / 128 * Dt * s3 + 3. / 32 * (c3 - 1.)

def y3(Dt):
  s4 = np.sinh(2 * Dt)
  c4 = np.cosh(2 * Dt)
  return 1. / 768 * Dt * Dt * Dt * (c4 + 8.) + 5. / 256 * Dt * Dt * s4 + 1. / 1024 * Dt * (57. * c4 - 64.) + 7. / 2048 * s4

def y4(Dt):
  s5 = np.sinh(2 * Dt)
  c5 = np.cosh(2 * Dt)
  return 1. / 12288 * Dt * Dt * Dt * Dt * (c5 + 16.) + 25. / 12288 * Dt * Dt * Dt * s5 + 1. / 16384 * Dt * Dt * (205. * c5 - 320.) + 435. / 32768 * Dt * s5 - 5. / 512. * (c5 - 1.)

def y5(Dt):
  s6 = np.sinh(2 * Dt)
  c6 = np.cosh(2 * Dt)
  return 1. / 245760 * Dt * Dt * Dt * Dt * Dt * (c6 + 32.) + 5. / 32768 * Dt * Dt * Dt * Dt * s6 + 1. / 65536 * Dt * Dt * Dt * (107. * c6 - 256.) + 316. / 65536 * Dt * Dt * s6 - 1. / 262144 * Dt * (279. * c6 - 2304.) - 2025. / 524288 * s6



def calculate_orders(beta,theta,t_s_scale):
  tau_split = t_s_scale * beta;
  tau_max   = 1.00 * beta;

  B   = tau_max * theta / 2.;
  t_s = tau_split * theta / 2.;

  tmp = np.cosh(B - t_s) * np.cosh(t_s - 0.0)
  order0 = tmp * tmp
  order1 = 2 * y1(t_s) * y1(B - t_s)
  order2 = 4 * (y2(t_s) * y2(B - t_s))
  order3 = -8 * (y1(t_s) * y5(B - t_s) + y2(t_s) * y4(B - t_s) + y4(t_s) * y2(B - t_s) + y5(t_s) * y1(B - t_s))
  full = np.cosh(B) * np.cosh(B)

  
  return full, order0, order1, order2, order3



beta = 2.0
theta = 3.0
t_s_scale = 0.5

N = 100

v_x = []

v_full = []
v_full0 = []
v_full1 = []
v_full2 = []
v_full3 = []

for ii in range(N):
  theta_ii = theta*(ii+1)/float(N+2)
  #t_s_scale_ii = t_s_scale#*(ii+1)/float(N+2)
  
  full, order0, order1, order2, order3 = calculate_orders(beta,theta_ii,t_s_scale)
  full0 = full  - order0
  full1 = full0 - order1
  full2 = theta_ii**6
  full3 = order3
  print full, full0, full1, full2, full3
  v_full.append(np.abs(full))
  v_full0.append(np.abs(full0))
  v_full1.append(np.abs(full1))
  v_full2.append(np.abs(full2))
  v_full3.append(np.abs(full3))
  
  v_x.append(theta_ii)

#print v_full3

print v_x 

#ax.plot(v_x, v_full0, label='order 0');
#ax.plot(v_x, v_full1, label='order 0+1');
ax.plot(v_x, v_full2, label='order 0+1+2');
ax.plot(v_x, v_full3, label='order 0+1+2+3');
ax.legend(loc='lower right')

plt.savefig("analytical_spinless_orders.pdf")
plt.show()

print "order 0:       % 4.8f " % order0
print "order 1:       % 4.8f " % order1
print "order 2:       % 4.8f " % order2
print "order 3:       % 4.8f " % order3
print "order 0+1+2+3: % 4.8f " % (order0 + order1 + order2 + order3)
print "\n\ncosh^2(B): % 4.8f \n" % (full)



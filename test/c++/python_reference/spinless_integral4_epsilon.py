# Copyright (c) 2020--present, The Simons Foundation
# This file is part of inchworm and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

import numpy as np

def y1(a,t0,t2):
  def e(a2):
    return (np.exp(a2*t2) - np.exp(a2*t0)) / a2
  
  e2 = np.exp(2.*t2)
  e0 = np.exp(2.*t0)
  
  em2 = np.exp(-2.*t2)
  em0 = np.exp(-2.*t0)
  
  return (1/16.) * ( e(a) * (e2*em0 + e0*em2 + 1)  + e(a-4)*e2*e0 + e(a+4)*em2*em0 + e(a-2)*(e2+e0) + e(a+2)*(em2+em0) )


def calculate_orders(beta, theta, t_s_scale, epsilon):
  tau_split = t_s_scale * beta;
  tau_max   = 1. * beta;

  B   = tau_max * theta / 2.;
  t_s = tau_split * theta / 2.;
  a   = 2. * epsilon / theta;
  prefactor = -4. / (1.+np.exp(-epsilon * beta))

  tmp = np.cosh(B - t_s) * np.cosh(t_s)
  order0 = tmp * tmp
  order1 = prefactor * y1(a,0,t_s) * y1(-a,t_s,B)
  full = np.cosh(B) * np.cosh(B)
  
  return full, order0, order1

beta = 2.0
theta = 4.0
t_s_scale = 0.5
epsilon = 20.

full, order0, order1 = calculate_orders(beta,theta,t_s_scale,epsilon)
full0 = full  - order0
full1 = full0 - order1
print full, full0, full1

print "order 0:       % 4.8f " % order0
print "order 1:       % 4.8f " % order1
print "order 0+1:     % 4.8f " % (order0 + order1)
print "\n\ncosh^2(B): % 4.8f \n" % (full)
exit()




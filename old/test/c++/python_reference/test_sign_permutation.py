#!/usr/bin/python
# zlib license:

# Copyright (c) 2019 Maxime Charlebois

# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.

# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:

# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import numpy as np
import numpy.linalg as la
import os,re,sys
from copy import deepcopy
from itertools import permutations 


t1 = [0.4,0.701,0.9,0.4,0.5]  #o
t2 = [0.05,0.305,0.8,0.55,0.23] #x

t1 = [0.4,0.701,0.9,0.4]  #o
t2 = [0.05,0.305,0.8,0.23] #x

#t1 = [0.2,0.3,0.4]  #o
#t2 = [0.22,0.25,0.5] #x

print len(t1), len(t2)

verbose1=1

def hyb(dtau):
  #return tau +0.5
  return (2.2 + dtau + 0.7*dtau*dtau + 0.1*dtau*dtau*dtau)

class DiagramTimes:
  def __init__(self, t1,t2,subrange=()):
    assert(len(t1)==len(t2))
    
    diag0 = len(t1)*'o' + len(t1)*'x'
    
    t1.sort()
    t2.sort()
    t = t1+t2
    t_indices = 2* range(len(t1))
    
    tmp_order = ''.join([x for _,x in sorted(zip(t,diag0))])
    tmp_t_indices = np.array([x for _,x in sorted(zip(t,t_indices))])
    #t.sort()
    t.sort()

    if subrange == ():
      self.order = tmp_order
      self.t_indices = tmp_t_indices
      #t.sort()
      self.t = t
      self.t1 = t1
      self.t2 = t2
    else:
      if verbose1: print "subrange ", subrange
      self.order = tmp_order[subrange[0]:subrange[1]+1]
      self.t_indices = tmp_t_indices[subrange[0]:subrange[1]+1]
      #t.sort()
      self.t = t[subrange[0]:subrange[1]+1]
    
      t1sub = []
      t2sub = []
      
      for index in range(len(self.t_indices)):
        time_index = self.t_indices[index]
        if self.order[index] == 'o':
          t1sub.append(t1[time_index])
          pos_cdag.append(index)
        else:
          t2sub.append(t2[time_index])
          pos_c.append(index)
      if len(t1sub)!=len(t2sub):
        print 'error - diagram does not have same number of x and o'
        exit()
      self.t1 = t1sub
      self.t2 = t2sub

    pos_c = []
    pos_cdag = []
    for index in range(len(self.t)):
      if self.order[index] == 'o':
        pos_cdag.append(index)
      else:
        pos_c.append(index)

      self.pos_c = pos_c
      self.pos_cdag = pos_cdag
    self.kOrder = len(self.t)/2



diagram = DiagramTimes(t1,t2)

def printDiag():
 string1 = ''
 for j in range(2*diagram.kOrder):
   string1 += diagram.order[j]+'-'
 s1=string1[:-1]
 print s1

printDiag()
print diagram.pos_cdag
    
print  np.array(diagram.pos_c).sum()%2
print  np.array(diagram.pos_c).sum()
print  np.array(diagram.pos_cdag).sum()
exit()



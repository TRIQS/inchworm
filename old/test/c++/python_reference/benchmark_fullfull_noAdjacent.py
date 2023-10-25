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


t1 = [0.4,0.801,0.9,0.4,0.5]  #o
t2 = [0.05,0.305,0.8,0.55,0.23] #x

#t1 = [0.4,0.501,0.9]  #o
#t2 = [0.05,0.305,0.2] #x

#t1 = [0.1,0.3,0.5,0.7,0.9] #o
#t2 = [0.21,0.45,0.69,0.81,1.0] #x

#t1 = [0.1,0.3,0.5,1.7,1.9] #o
#t2 = [1.21,1.45,1.69,1.81,2.0] #x

#t1 = [0.0,0.2,0.4,1.0,1.4] #o
#t2 = [1.1,1.3,1.5,1.9,2.0] #x



#t1 = [0.2,0.4,1.0,1.5] #o
#t2 = [0.0,0.5,1.9,2.0] #x


#t1 = [0.4344,0.1,0.3,0.5,0.55,0.566,0.33,0.99839,0.333] #o
#t2 = [0.0,0.3452,0.21,0.45,0.81,0.122,0.833,0.513,0.2314] #x

#t1 = [0.0,0.05,4.0,2.3]  #o
#t2 = [0.2,0.5,0.8,1.2] #x

#t1 = [0.1,1.3,1.7,1.9] #o
#t2 = [1.21,1.69,1.81,2.0] #x

#t1 = [0.1,1.3,1.7,1.9,2.05] #o
#t2 = [1.21,1.69,1.81,2.0,2.1] #x

#t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33] #o
#t2 = [0.0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833] #x

#t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45] #o
#t2 = [0.0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4] #x

t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45,2.3,1.6,4.3,1.222,2.98,3.1244] #o
t2 = [0.0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4,3.4,2.34,3.11,2.9,1.99,3.098] #x


t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941] #o
t2 = [0.0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111] #x


t1 = [0.4344,0.1,0.3,0.5,0.75,0.9]
t2 = [0,0.3452,0.21,0.45,0.69,0.81]

print len(t1), len(t2)
#t1 = [0.2,0.4,1.0,1.6] #o
#t2 = [1.1,1.3,1.5,1.9] #x

#t1 = [0.9,0.4,0.5]  #o
#t2 = [0.8,0.6,0.2] #x

#t1 = [0.55,0.33]     #o
#t2 = [0.81,0.122] #x

#t1 = [0.3,0.5,0.9]     #o
#t2 = [0.0,0.4,0.8] #x

#t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566]     #o
#t2 = [0.0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122] #x

#t1 = [0.4344,0.1,0.3 , 0.9]     #o
#t2 = [0.0,0.3452,0.21, 0.7] #x



#t1 = [0.4344,0.1,0.3,0.5,0.75,0.9,0.55,0.566,0.33,.4959594,.4494929,.12349512,.62343,0.123412444,0.2134444,.99949941,1.1,1.23,1.45,2.3,1.6,4.3,1.222,2.98,3.1244]
#t2 = [0,0.3452,0.21,0.45,0.69,0.81,0.998,0.122,0.833,0.4934,.210342134,.210343,.02134,.0030404,.02142430,0.1111,1.2,1.3,1.4,3.4,2.34,3.11,2.9,1.99,3.098]

#t1 = [1.0,2.0,4.0,4.5]
#t2 = [0.8,1.8,3.6,4.7]


t1 = [0.1, 0.2, 0.3, 0.5, 1.3, 1.4, 1.5, 3.6, 3.78]
t2 = [0.0, 0.53, 0., 1.2, 1.8, 2.2, 2.3, 4.6, 4.7]

t1 = [0.1, 0.2, 0.3, 0.5, 1.3, 1.4, 1.5, 3.6, 3.78, 3.88]
t2 = [0.0, 0.53, 0., 1.2, 1.8, 2.2, 2.3, 4.6, 4.7, 4.9]

#t1 = [0.1, 0.2, 0.3, 0.5, 1.3, 1.4, 1.5]
#t2 = [0.0, 0.53, 0., 1.2, 1.8, 2.2, 2.3]

#t1=[1.0,2.0,3.0,4.5];
#t2=[1.2,1.8,4.6,4.7];

doProper = 1
split_point = 5
verbose1=1
SMALLEST_SEGMENT=4

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
        else:
          t2sub.append(t2[time_index])
      if len(t1sub)!=len(t2sub):
        print 'error - diagram does not have same number of x and o'
        exit()
      self.t1 = t1sub
      self.t2 = t2sub
    self.kOrder = len(self.t)/2


class Delta:
  def __init__(self, diagram):
    delta_matrix = np.zeros([diagram.kOrder,diagram.kOrder],dtype='float')
    for it1 in range(len(diagram.t1)):
      for it2 in range(len(diagram.t2)):
        delta_matrix[it1,it2] = hyb(diagram.t2[it2]-diagram.t1[it1])
    #print delta_matrix
    if SMALLEST_SEGMENT==4:
      for ii in range(len(diagram.t)-1):
       if(ii != len(diagram.t)-split_point-1):
        if (diagram.order[ii] == 'x') and (diagram.order[ii+1] == 'o'): #segment of length 2
            it1 = diagram.t_indices[ii]
            it2 = diagram.t_indices[ii+1]
            delta_matrix[it1,it2] = 0.0
        elif (diagram.order[ii] == 'o') and (diagram.order[ii+1] == 'x'): #segment of length 2
            it1 = diagram.t_indices[ii+1]
            it2 = diagram.t_indices[ii]
            delta_matrix[it1,it2] = 0.0
    self.matrix = delta_matrix
    self.diagram = diagram

  def prints(self,diagram):
    for i in range(diagram.kOrder):
      for j in range(diagram.kOrder):
        print "% 5.3f" % self.matrix[i,j],
      print 
    print 

  def det(self,list_of_indices): #very dirty, need to find a more clever way
    
    it1 = []
    it2 = []
    
    for index in list_of_indices:
      time_index = self.diagram.t_indices[index]
      if self.diagram.order[index] == 'x':
        #print 'o'
        it1.append(time_index)
      else:
        #print 'x'
        it2.append(time_index)
    if len(it1)!=len(it2):
      print 'error - diagram does not have same number of x and o'
      exit()
    
    value = la.det(self.matrix[np.ix_(it1,it2)])
    return value














def inclusion_exclusion(t1,t2,split_point):
  diagramTime = DiagramTimes(t1,t2)
  delta = Delta(diagramTime)
  
  print diagramTime.order
  delta.prints(diagramTime)
  #exit()

  kOrder = diagramTime.kOrder 
  diagram = diagramTime.order
  #print kOrder

  def printPos():
    string1 = ''
    for j in range(2*kOrder):
      if(2*kOrder - split_point -1 == j):
        string1 += ' |'
      else:
        string1 += '  '
      
    print string1

  def printDiag():
    printPos()
    string1 = ''
    for j in range(2*kOrder):
      string1 += diagram[j]+'-'
    s1=string1[:-1]
    print s1


  def printLine(segmentarray,chars='-=', additionnal_string = ''):
    current_segment=0
    string1=''
    for ii in range(len(segmentarray)):
     if((current_segment != segmentarray[ii]) or (ii==0)):
       current_segment = segmentarray[ii]
       if(current_segment==-1):
         char1=' '
       elif(current_segment==0):
         char1=chars[0]
       else:
         char1=chars[1]  
       if len(string1)>0: string1=string1[:-1]+' '
     string1 += char1+char1
    s1=string1[:-1]
    print s1+additionnal_string 


  def removeSubList_fromList(list1,sublist1):
    return [x for x in list1 if x not in sublist1]

  class Segment:
    def __init__(self, size,position):
      self.sz =   size
      self.pos1 = position
      self.pos2 = position+size
      self.value = 0.0
      self.value_without_cuts = 0.0
      self.calculated = False
      self.cuts = []   #Combination_of_Segment number
      self.subs = []   #Combination_of_Segment number
      self.number = 0   #to be compatible with Combination_of_Segment


    def calculate(self,seg_J,combinaison,verbose=0):
      #print 'calculating diag ', self.pos1, self.pos2
      if verbose: printLine(formatted_output(seg_J, combinaison[self.number]),'-=')
              
      self.calculated = True
      rangeOfVertex = range(self.pos1,self.pos2)
      
      self.value = delta.det(rangeOfVertex)
      if verbose: print rangeOfVertex, self.value
      
      if len(self.subs) >0:
        for sub in self.subs:
          if verbose: printLine(formatted_output(seg_J, combinaison[sub]),'-=')
              
          value = 1.0
          rangeOfSubVertex = deepcopy(rangeOfVertex)
          signe_de_parcollet_charlebois = 1.0
          for segNo in combinaison[sub].list:
            if seg_J[segNo].calculated == False:
              print 'error: sub not calculated?'
              exit()
            value *= -seg_J[segNo].value
            #print segNo, -seg_J[segNo].value
            #print rangeOfSubVertex, range(seg_J[segNo].pos1,seg_J[segNo].pos2)
            rangeOfSubVertex = removeSubList_fromList(rangeOfSubVertex,range(seg_J[segNo].pos1,seg_J[segNo].pos2))
            #print rangeOfSubVertex, range(seg_J[segNo].pos1,seg_J[segNo].pos2)
            if (seg_J[segNo].sz % 4 !=0): 
              if (seg_J[segNo].pos2-self.pos1) %2 ==1:
                signe_de_parcollet_charlebois *= -1.0
          
          
          #if len(combinaison[sub].list) >0:
          #  last_segment_end = seg_J[segNo].pos2
            
            
          if rangeOfSubVertex !=[]:
            det1 = delta.det(rangeOfSubVertex)
            if verbose: print rangeOfSubVertex, det1
            value *= signe_de_parcollet_charlebois * det1
            #print 'rest', signe_de_parcollet_charlebois * det1
            #if signe_de_parcollet_charlebois <0:
            #  exit()
          if verbose: print "value %f\n" % value
          self.value += value
        if verbose: print "self.value %f " % self.value
      
      self.value_without_cuts = self.value
      if len(self.cuts) >0:
        for cut in self.cuts:
          if verbose: printLine(formatted_output(seg_J, combinaison[cut]),'-=')
          #print cut
          #print combinaison[cut].list
          value = 1.0
          
          for segNo in combinaison[cut].list:
            if seg_J[segNo].calculated == False:
              print 'error: sub not calculated?'
              exit()
            value *= -seg_J[segNo].value_without_cuts
            #print 'seg.value',seg_J[segNo].value
          #print 'value ', value
          #value *= delta.det(rangeOfSubVertex)
          self.value -= value
        #print "self.value %f" % self.value
        
      if verbose: print 'self.value %f \n' % self.value
      if verbose: print "\n\n\n\n" 
      
          
  class Combination_of_Segment:
    def __init__(self, segment_number,seg_J):
      segment0  = seg_J[segment_number]
      self.sz   = segment0.pos2-segment0.pos1
      self.pos1 = segment0.pos1
      self.pos2 = segment0.pos2
      self.value = 0.0
      self.fullyAdjacent = True
#      self.fullyDisjoint = True
      self.adjacent = True
      self.disjoint = True
      self.special = False
      self.nSegment = 1 
      self.list = [segment_number]
      self.calculated = False

    def append(self, segment_number, seg_J):
      segment0  = seg_J[segment_number]
      #if self.adjacent: # test if it is still adjacent or not
      if segment0.pos1 != self.pos2:
        self.fullyAdjacent = False
      elif (self.pos2 == 2*kOrder-split_point):
        self.fullyAdjacent = False
#      elif segment0.pos1 == split_point:
#        self.fullyAdjacent = False
      else: 
        self.disjoint = False
          
      
      self.sz   = segment0.pos2-self.pos1
      self.pos2 = segment0.pos2
      
      self.list.append(segment_number)
      self.nSegment = len(self.list)
      if not not_segment_cross_p(2*kOrder-split_point, self.pos1,self.pos2):
        self.special = True
      

  start_index = np.zeros([2*kOrder],dtype='int')

  def not_segment_cross_p(point, a,b):
    return (a-(point))*(b-(point))>=0

  def determine_indep_segments(split_point,diagram1):
    seg_J=[]
    assert(split_point<len(diagram1))
      #print aa,
    for ii in range(len(diagram1)-1): # starting position of segment j
      start_index[ii] = len(seg_J)
      for aa in range(len(diagram1)-ii+1): # lenght of segments
        if aa >= SMALLEST_SEGMENT: #because we will find another way to filter out segment of size 2.
          #print ii,
          if not_segment_cross_p(2*kOrder-split_point, ii,ii+aa):
           if(diagram1[ii:(ii+aa)].count('o') == diagram1[ii:(ii+aa)].count('x')):
            seg_tmp = Segment(aa,ii)
            
            seg_J.append(Segment(aa,ii))
            seg_J[len(seg_J)-1].number = len(seg_J)-1 
    start_index[ii+1] = len(seg_J)
    return seg_J

  def combine_segments(seg_J):
    nCombination = 0
    nSingleSeg = len(seg_J)
    
    start_index_number_j = [0,0]
    combination_of_seg = []
    
    
    for ii in range(nSingleSeg):
      combination_of_seg.append(Combination_of_Segment(ii,seg_J))
      for jj in range(len(seg_J)):
        if((seg_J[ii].pos1 >= seg_J[jj].pos1) and (seg_J[ii].pos2 < seg_J[jj].pos2)):
          if nCombination!=jj:
            seg_J[jj].subs.append(nCombination)
      
      nCombination +=1

    for number_of_segment in range(kOrder+1): #we know maximum of potential segments of length 4 or more.
      if number_of_segment > 1:
        start_index_number_j.append(len(combination_of_seg))
        for previous_combination_of_seg_index in range(start_index_number_j[-2],start_index_number_j[-1]):
          
          previous_combination_of_seg = combination_of_seg[previous_combination_of_seg_index]
          #print previous_combination_of_seg,combination_of_seg
          
          last_segment = seg_J[previous_combination_of_seg.list[-1]]
          if last_segment.pos2  < len(start_index):
            for additionnal_segment in range(start_index[last_segment.pos2],start_index[-1]):
              comb_of_seg = deepcopy(previous_combination_of_seg)
              comb_of_seg.append(additionnal_segment,seg_J)
              #print previous_combination_of_seg
              combination_of_seg.append(comb_of_seg)

              #print comb_of_seg.pos1,  comb_of_seg.pos2
              #if not comb_of_seg.fullyAdjacent:
#              for ii in range(len(seg_J)):
#                #print '' , seg_J[ii].pos1, seg_J[ii].pos2
#                if((seg_J[ii].pos1 <= comb_of_seg.pos1) and (seg_J[ii].pos2 > comb_of_seg.pos2)):
#                  #print 'accepted'
#                  seg_J[ii].subs.append(nCombination)

              
              ##analyse this segment:
              if comb_of_seg.fullyAdjacent:
                for ii in range(len(seg_J)):
                  if((seg_J[ii].pos1 == comb_of_seg.pos1) and (seg_J[ii].pos2 == comb_of_seg.pos2)):
                    seg_J[ii].cuts.append(nCombination)
                    
              elif(comb_of_seg.disjoint):
                for ii in range(len(seg_J)):
                  if((seg_J[ii].pos1 <= comb_of_seg.pos1) and (seg_J[ii].pos2 > comb_of_seg.pos2)):
                    seg_J[ii].subs.append(nCombination)
              
              nCombination +=1
    #exit()
    return combination_of_seg


  def formatted_output(seg_J, comb_of_seg):
    segmentarray = np.zeros(len(diagram),dtype='int')
    NN = 0
    for segment_no in comb_of_seg.list:
      NN+=1
      seg = seg_J[segment_no]
      for a in range(seg.sz):
        segmentarray[a+seg.pos1] = NN
    return segmentarray
    #printLine(segmentarray,'-=',adjacent)
  

  def enumerate_segments(seg_J):
    printDiag()
    for seg in seg_J:
      #print seg.pos1, seg.pos2
      segmentarray = np.zeros(len(diagram),dtype='int')
      for a in range(seg.sz):
        segmentarray[a+seg.pos1] = 1
      value_string = '  '
      if seg.calculated:
        value_string += ' val=% 15.5f   no_cuts=% 15.5f ' % (seg.value, seg.value_without_cuts)
      printLine(segmentarray,'-=',value_string)



  def enumerate_combination_of_segments(seg_J,combination_of_seg):
    printDiag()
    segNo=0
    for comb_of_seg in combination_of_seg:
      segmentarray = formatted_output(seg_J, comb_of_seg)
      
      adjacent = '  '+str(segNo)
      if comb_of_seg.fullyAdjacent:
        adjacent += '  A'
#      if comb_of_seg.fullyDisjoint:
#        adjacent += '  D'
#      if comb_of_seg.adjacent:
#        adjacent += '  a'
      if comb_of_seg.disjoint:
        adjacent += '  d'
#      if comb_of_seg.special:
#        adjacent += '  special'
      
        printLine(segmentarray,'-=',adjacent)
        segNo+=1


    


  def enumerate_subdiagram_of_single_segments(seg_J,combination_of_seg,use_cuts_instead):
    printDiag()
    for ii in range(len(seg_J)):
      seg = seg_J[ii]

      segmentarray = np.zeros(len(diagram),dtype='int')
      for a in range(seg.sz):
        segmentarray[a+seg.pos1] = 1
      adjacent = '  '+str(ii)
      printLine(segmentarray,'-=',adjacent)
      
      list_to_sparse = seg_J[ii].subs
      if use_cuts_instead:
        list_to_sparse = seg_J[ii].cuts
        
      
      for comb_of_seg_number in list_to_sparse:
        comb_of_seg = combination_of_seg[comb_of_seg_number]
        
        segmentarray = np.zeros(len(diagram),dtype='int')-1
        for jj in range(seg_J[ii].pos1,seg_J[ii].pos2):
          segmentarray[jj] = 0
          
        NN = 0
        value = 1.0
        values = ''
        for segment_no in comb_of_seg.list:
          NN+=1
          seg = seg_J[segment_no]
          for a in range(seg.sz):
            segmentarray[a+seg.pos1] = NN
          
          values += ' %10.3f' % (-seg.value_without_cuts)
          
          if seg.calculated == True:
            value *= -seg.value_without_cuts
        
        values = '  %10.5f =' % (value) + values
        #print segmentarray
        printLine(segmentarray,'-=',values)
        #segNo+=1        
      print
      


  seg_J = determine_indep_segments(split_point,diagram)
  print '\ndifferent segments:'
  enumerate_segments(seg_J)
  
  combination_of_seg = combine_segments(seg_J)

  if verbose1>1: 
    print '\ndifferent combination of segments:'
    enumerate_combination_of_segments(seg_J,combination_of_seg)
    print '\n\n\ndifferent subdiagram of single segments:'
    enumerate_subdiagram_of_single_segments(seg_J,combination_of_seg,False)

    print '\n\n\ndifferent cuts of single segments:'
    enumerate_subdiagram_of_single_segments(seg_J,combination_of_seg,True)

    print
    print "## diagram = '%s'" % diagram 
    print "kOrder = %d" % kOrder
    print "number of segments = %d" % len(seg_J)
    print "number of combinations = %d" % len(combination_of_seg)

  #calculate segment starting from the smallest
  for length in range(2,2*kOrder,2):
    if verbose1: print '\n\n############\nsegment length=', length

    for seg in seg_J:
      if seg.sz == length:
        seg.calculate(seg_J,combination_of_seg)
  
  full_seg = Segment(2*kOrder,0)
  for ii in range(len(combination_of_seg)):
    comb_of_seg = combination_of_seg[ii]
    if comb_of_seg.disjoint:# or (comb_of_seg.special):
      if verbose1: printLine(formatted_output(seg_J, comb_of_seg),'-=')
      full_seg.subs.append(ii)
      
  full_seg.calculate(seg_J,combination_of_seg)
  print 'c_k = %f' %full_seg.value
  

  if verbose1: 

    print '\ndifferent segments:'
    enumerate_segments(seg_J)
    
    print '\ndifferent combination of segments:'
    enumerate_combination_of_segments(seg_J,combination_of_seg)


  print
  print "## diagram = '%s'" % diagram 
  print "kOrder = %d" % kOrder
  print "number of segments = %d" % len(seg_J)
  print "number of combinations = %d" % len(combination_of_seg)

  
  return full_seg.value



























def proper_enum(t1,t2,split_point1,subrange=()):

  diagramTime = DiagramTimes(t1,t2,subrange)
  delta = Delta(diagramTime)

  #print diagramTime.order
  #print delta.matrix

  kOrder1 = diagramTime.kOrder 
  diagram1 = diagramTime.order

  #print kOrder1
  #print diagram1

  kRange = range(kOrder1)

  pos_o = [m.start() for m in re.finditer('o', diagram1)]
  pos_x = [m.start() for m in re.finditer('x', diagram1)]
  assert(len(pos_o) == len(pos_x))
  #print pos_o, pos_x

  def printPos():
    string1 = ''
    for j in range(2*kOrder1):
      if(2*kOrder1 - split_point -1 == j):
        string1 += ' |'
      else:
        string1 += '  '
      
    print string1

  def printDiag():
    string1 = ''
    for j in range(2*kOrder1):
      string1 += diagram1[j]+'-'
    s1=string1[:-1]
    print s1


  def order(a,b):
    if a>b:
      tmp=b
      b=a
      a=tmp
    return a,b



  def printLine(a,b,char1='.'):
    assert(a!=b)
    a,b = order(a,b)
    #print a,b
    string1 = ''
    
    if(a<=b):
     for ii in range(2*kOrder1):
      if((ii>=a) and (ii<=b)):
        string1 += char1 + char1
      else:
        string1 += '  '
    s1=string1.replace(char1+' ','  ')[:-1]
    
    print s1


  #for jj in range(6):
  #  for ii in range(6):
  #    if(ii!=jj):
  #      printLine(jj,ii)
  #exit()


  def point_overlap(point,a,b):
    return (a-point)*(b-point)>0

  def segment_cross(a1,b1,a2,b2):
    return (a1-a2)*(b1-a2)*(b1-b2)*(a1-b2)<0

  def segment_before_point_overlap(point,a,b):
    #print min(a,b)
    return not((min(a,b)-point)>0)

  def segment_after_point_overlap(point,a,b):
    return (max(a,b)-point)>0

  def segment_cross_p(point, a,b):
    return (a-(point+0.5))*(b-(point+0.5))<0

  
  def print_graph(permutation,cross_split_point,visited):
    
    printPos()
    printDiag()

    for arch in range(kOrder1): #as much arches than pair of vertices
      char1 = '.'
      if(cross_split_point[arch]):
        char1='_'
      if(visited[arch]):
        char1='='
      #elif()
      printLine(pos_o[arch],pos_x[permutation[arch]],char1)
            
    print ''

  NN = 0



  #print
  #print permutation0
  #print


  #print cross_split_point_pile
  def grow_pile(arch,permutation,connexion_pile,visited):
    #print_graph()
    connexion_pile.append(arch)
    o = pos_o[arch]
    x = pos_x[permutation[arch]]
    for arch2 in range(kOrder1):
      if not visited[arch2]:      
        if segment_cross(o,x,pos_o[arch2],pos_x[permutation[arch2]]):
          visited[arch2]=True
          grow_pile(arch2,permutation,connexion_pile,visited)
    connexion_pile.pop()
    visited[arch]=True

  def test_diagram_connection(permutation,verbose=1):

    cross_split_point = np.zeros(kOrder1,dtype='bool')
    visited = np.zeros(kOrder1,dtype='bool')
    connexion_pile = []  
    cross_split_point_pile = []

    for arch in range(kOrder1):
      #print arch,kOrder1
      o = pos_o[arch]
      x = pos_x[permutation[arch]]
      #if(segment_before_point_overlap(split_point,o,x)):
      #if(segment_after_point_overlap(2*kOrder1-split_point1-1,o,x)):
      if segment_cross_p(2*kOrder1-split_point1-1,o,x):
        cross_split_point[arch] = True
        #visited[arch] = True
        cross_split_point_pile.append(arch)
    
    if verbose>1: 
      print_graph(permutation,cross_split_point,visited)
    
    for arch in cross_split_point_pile:
      if(not visited[arch]):
        visited[arch] = True
        grow_pile(arch,permutation,connexion_pile,visited)

    if verbose>1: print_graph(permutation,cross_split_point,visited)
    #if verbose>0: print visited.all()  
    return visited.all()

  #print test_diagram_connection(permutation0,verbose=1)

  # Get all permutations of [1, 2, 3] 
  perm = permutations(range(kOrder1))








  def parity(perm):
    perm0 = list(perm)
    par = 1
    for i in range(len(perm0)):
      if perm0[i] != i:
        par *= -1
        swap_i = perm0.index(i)
        perm0[i],perm0[swap_i] = perm0[swap_i],perm0[i]
    return par    


  #print_graph()
  NN=0
  N_proper = 0
  value = 0.0


  for perm0 in perm:
    #print perm0
    #continue
    #print '%d: ' % (NN+1)
    NN +=1
    #print i
    #printPos()
    #printDiag()
      
    #for j in range(kOrder):
    #  printLine(pos_o[j],pos_x[perm0[j]]) 
    test = test_diagram_connection(perm0,verbose=verbose1)
    value1 = 1.0
    for k in range(kOrder1):
      value1 *= delta.matrix[perm0[k],k]
    
    pari1 = parity(perm0)
  
    if test:
      N_proper +=1
      value += pari1*value1
      #print '%d: ' % (NN+1),
      #print perm0
      
      #printPos()
      #printDiag()
      #test = test_diagram_connection(perm0,verbose=2)
      
      if verbose1 > 1: print 'proper, value=%f, parity=%d\n' %(value1,pari1)
    else:
      if verbose1 > 1: print 'IMPROPER, value=%f, parity=%d\n' %(value1,pari1)
    if verbose1 > 1: print '\n\n'

  print
  print "## diagram = '%s'" % diagram1
  print "kOrder = %d" % kOrder1
  print "number of diagram = %d" % NN
  print "number of proper diagram = %d" % N_proper
  print
  print 'c_k = %f' %value
  
  return value
  
  
if (doProper): val1 = proper_enum(t1,t2,split_point)
val2 = inclusion_exclusion(t1,t2,split_point)

if (doProper): print 'proper_enum():         c_k = %f' %val1
print 'inclusion_exclusion(): c_k = %f' %val2

if (doProper): print 'difference = %e' %(val2-val1) 
    
    
exit()



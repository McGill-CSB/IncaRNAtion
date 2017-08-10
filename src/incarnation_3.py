###############################################################################
#Copyright (c) 2013, Vladimir Reinharz, Yann Ponty & Jerome Waldispuhl        #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#* Redistributions of source code must retain the above copyright             #
#notice, this list of conditions and the following disclaimer.                #
#* Redistributions in binary form must reproduce the above copyright          #
#notice, this list of conditions and the following disclaimer in the          #
#documentation and/or other materials provided with the distribution.         #
#* Neither the name of the <organization> nor the                             #
#names of its contributors may be used to endorse or promote products         #
#derived from this software without specific prior written permission.        #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
#ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY       #
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT   #
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS#
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE                  #
###############################################################################
import itertools
import random
import math
import sys
import os

def MPMATH_MISSING():
  print("""The module `mpmath` was not found. This might impeed the
  processing of long RNAs.
  More information can be found on http://code.google.com/p/mpmath/
  We also recommand installation of `gmpy` (http://code.google.com/p/gmpy/),
  automatically leveraged by `mpmath` to increase the speed of computations""")

try: #For infinite precision
  from mpmath import mpf
except ImportError:
  MPMATH_MISSING()
  def mpf(n):
    return n

sys.setrecursionlimit(10000)

IUPACBASES = {
  'A':['A'],
  'C':['C'],
  'G':['G'],
  'U':['U'],
  'N':['A','C','G','U'],
  'R':['A','G'],
  'Y':['C','U'],
  'S':['G','C'],
  'W':['A','U'],
  'K':['G','U'],
  'M':['A','C'],
  'B':['C','G','U'],
  'D':['A','G','U'],
  'H':['A','C','U'],
  'V':['A','C','G']
  }
BASES = []

BOLTZMANN = 0.0019872041
global T
T = 310.15

#Just generate all possible combination with maxint
global STACKING_ENERGY
STACKING_ENERGY = {k:sys.maxsize for k in itertools.product(
                  IUPACBASES['N'], repeat=4)}
#Adjust with turner04
#The order of the nucleotides is from 5' -> 3'
STACKING_ENERGY.update({('A', 'A', 'U', 'U'):-0.9,
                        ('A', 'C', 'G', 'U'):-2.2,
                        ('A', 'G', 'C', 'U'):-2.1,
                        ('A', 'G', 'U', 'U'):-0.6,
                        ('A', 'U', 'A', 'U'):-1.1,
                        ('A', 'U', 'G', 'U'):-1.4,
                        ('C', 'A', 'U', 'G'):-2.1,
                        ('C', 'G', 'U', 'G'):-1.4,
                        ('C', 'U', 'A', 'G'):-2.1,
                        ('C', 'U', 'G', 'G'):-2.1,
                        ('G', 'A', 'U', 'C'):-2.4,
                        ('G', 'C', 'G', 'C'):-3.4,
                        ('G', 'G', 'C', 'C'):-3.3,
                        ('G', 'G', 'U', 'C'):-1.5,
                        ('G', 'U', 'A', 'C'):-2.2,
                        ('G', 'U', 'G', 'C'):-2.5,
                        ('G', 'A', 'U', 'U'):-1.3,
                        ('G', 'C', 'G', 'U'):-2.5,
                        ('G', 'G', 'C', 'U'):-2.1,
                        ('G', 'G', 'U', 'U'):-0.5,
                        ('G', 'U', 'A', 'U'):-1.4,
                        ('G', 'U', 'G', 'U'):1.3,
                        ('U', 'A', 'U', 'A'):-1.3,
                        ('U', 'C', 'G', 'A'):-2.4,
                        ('U', 'G', 'C', 'A'):-2.1,
                        ('U', 'G', 'U', 'A'):-1.0,
                        ('U', 'U', 'A', 'A'):-0.9,
                        ('U', 'U', 'G', 'A'):-1.3,
                        ('U', 'A', 'U', 'G'):-1.0,
                        ('U', 'C', 'G', 'G'):-1.5,
                        ('U', 'G', 'C', 'G'):-1.4,
                        ('U', 'G', 'U', 'G'):0.3,
                        ('U', 'U', 'A', 'G'):-0.6,
                        ('U', 'U', 'G', 'G'):-0.5,
                        ('C', 'C', 'G', 'G'):-3.3,
                        ('C', 'G', 'C', 'G'):-2.4})

#All isostericity values. GG->anything or anything->GG is 10
ISO = {(('A', 'C'), ('C', 'A')):4.93,
        (('G', 'C'), ('A', 'G')):3.5,
        (('G', 'C'), ('G', 'C')):0.0,
        (('U', 'G'), ('U', 'G')):0.0,
        (('A', 'C'), ('G', 'U')):0.8,
        (('A', 'G'), ('U', 'U')):8.18,
        (('G', 'A'), ('A', 'G')):2.25,
        (('A', 'A'), ('A', 'C')):4.58,
        (('A', 'C'), ('U', 'G')):4.76,
        (('U', 'G'), ('C', 'G')):2.14,
        (('U', 'G'), ('G', 'G')):10.0,
        (('U', 'C'), ('U', 'U')):4.36,
        (('G', 'G'), ('U', 'A')):10.0,
        (('C', 'G'), ('G', 'C')):0.26,
        (('A', 'U'), ('C', 'A')):2.75,
        (('A', 'U'), ('U', 'U')):3.8,
        (('U', 'U'), ('G', 'C')):3.94,
        (('U', 'G'), ('G', 'C')):2.39,
        (('C', 'A'), ('C', 'A')):0.0,
        (('G', 'A'), ('U', 'G')):3.8,
        (('A', 'G'), ('G', 'C')):4.38,
        (('A', 'A'), ('C', 'G')):3.44,
        (('G', 'G'), ('C', 'U')):10.0,
        (('G', 'C'), ('U', 'A')):0.34,
        (('U', 'A'), ('A', 'A')):4.66,
        (('U', 'U'), ('U', 'C')):6.46,
        (('U', 'U'), ('C', 'G')):3.8,
        (('C', 'C'), ('U', 'C')):8.25,
        (('G', 'A'), ('C', 'U')):0.0,
        (('C', 'U'), ('G', 'U')):5.06,
        (('C', 'G'), ('U', 'C')):3.44,
        (('U', 'C'), ('G', 'U')):2.89,
        (('U', 'U'), ('U', 'G')):2.89,
        (('G', 'A'), ('U', 'U')):5.97,
        (('G', 'U'), ('A', 'G')):4.44,
        (('A', 'C'), ('G', 'A')):5.05,
        (('G', 'C'), ('U', 'G')):2.39,
        (('G', 'G'), ('A', 'A')):10.0,
        (('U', 'U'), ('C', 'A')):2.39,
        (('G', 'G'), ('C', 'A')):10.0,
        (('A', 'A'), ('C', 'U')):1.53,
        (('G', 'U'), ('C', 'C')):6.25,
        (('C', 'C'), ('U', 'U')):2.37,
        (('C', 'U'), ('U', 'C')):7.97,
        (('C', 'C'), ('G', 'G')):10.0,
        (('G', 'G'), ('A', 'U')):10.0,
        (('C', 'G'), ('A', 'U')):0.34,
        (('C', 'U'), ('A', 'U')):5.3,
        (('U', 'A'), ('C', 'A')):2.47,
        (('U', 'U'), ('A', 'U')):3.8,
        (('G', 'G'), ('U', 'U')):10.0,
        (('U', 'A'), ('C', 'C')):5.3,
        (('G', 'G'), ('A', 'G')):10.0,
        (('A', 'G'), ('A', 'U')):4.52,
        (('A', 'A'), ('A', 'G')):2.33,
        (('U', 'C'), ('U', 'A')):3.8,
        (('G', 'U'), ('U', 'U')):5.27,
        (('A', 'A'), ('C', 'A')):5.3,
        (('U', 'G'), ('C', 'U')):3.8,
        (('G', 'U'), ('G', 'U')):0.0,
        (('G', 'A'), ('A', 'U')):3.57,
        (('U', 'A'), ('A', 'U')):0.31,
        (('U', 'A'), ('U', 'C')):3.57,
        (('U', 'G'), ('A', 'G')):4.33,
        (('C', 'G'), ('C', 'A')):2.55,
        (('U', 'A'), ('U', 'A')):0.0,
        (('G', 'G'), ('G', 'A')):10.0,
        (('U', 'U'), ('C', 'U')):5.97,
        (('A', 'C'), ('A', 'G')):5.14,
        (('U', 'C'), ('A', 'A')):6.91,
        (('A', 'C'), ('A', 'U')):2.47,
        (('G', 'C'), ('A', 'C')):2.55,
        (('C', 'C'), ('A', 'C')):5.91,
        (('A', 'U'), ('A', 'U')):0.0,
        (('C', 'C'), ('C', 'C')):0.0,
        (('U', 'C'), ('A', 'C')):2.39,
        (('G', 'A'), ('C', 'C')):7.97,
        (('A', 'A'), ('G', 'C')):3.39,
        (('U', 'C'), ('A', 'G')):6.96,
        (('G', 'C'), ('G', 'A')):3.49,
        (('U', 'U'), ('U', 'A')):3.63,
        (('G', 'U'), ('G', 'G')):10.0,
        (('C', 'G'), ('U', 'G')):2.14,
        (('A', 'G'), ('U', 'A')):4.66,
        (('A', 'A'), ('U', 'A')):3.57,
        (('U', 'G'), ('U', 'A')):2.11,
        (('C', 'A'), ('U', 'C')):5.3,
        (('G', 'A'), ('C', 'G')):3.39,
        (('U', 'A'), ('C', 'G')):0.21,
        (('A', 'U'), ('C', 'U')):3.57,
        (('G', 'C'), ('U', 'U')):3.94,
        (('U', 'U'), ('C', 'C')):2.37,
        (('G', 'U'), ('A', 'A')):4.1,
        (('A', 'C'), ('C', 'U')):5.3,
        (('A', 'G'), ('C', 'C')):9.77,
        (('U', 'A'), ('G', 'C')):0.34,
        (('A', 'U'), ('A', 'C')):2.47,
        (('U', 'G'), ('G', 'A')):4.44,
        (('G', 'U'), ('A', 'U')):2.11,
        (('A', 'U'), ('A', 'A')):4.52,
        (('C', 'U'), ('C', 'C')):3.02,
        (('C', 'U'), ('U', 'A')):5.39,
        (('A', 'C'), ('U', 'C')):4.58,
        (('C', 'G'), ('G', 'A')):3.5,
        (('C', 'A'), ('G', 'G')):10.0,
        (('G', 'A'), ('G', 'G')):10.0,
        (('C', 'G'), ('U', 'U')):3.8,
        (('U', 'A'), ('A', 'C')):2.75,
        (('G', 'A'), ('C', 'A')):4.58,
        (('G', 'A'), ('A', 'C')):5.3,
        (('A', 'C'), ('C', 'G')):2.78,
        (('U', 'C'), ('G', 'G')):10.0,
        (('A', 'A'), ('A', 'U')):3.5,
        (('C', 'C'), ('C', 'U')):7.97,
        (('A', 'G'), ('U', 'C')):2.71,
        (('A', 'C'), ('A', 'C')):0.0,
        (('U', 'A'), ('G', 'A')):3.67,
        (('C', 'C'), ('G', 'U')):6.25,
        (('A', 'U'), ('U', 'G')):2.4,
        (('U', 'C'), ('C', 'G')):3.94,
        (('U', 'C'), ('U', 'C')):5.97,
        (('A', 'G'), ('A', 'A')):0.0,
        (('A', 'A'), ('U', 'G')):4.59,
        (('C', 'C'), ('A', 'U')):5.39,
        (('A', 'G'), ('C', 'U')):3.77,
        (('C', 'C'), ('U', 'G')):5.06,
        (('G', 'U'), ('C', 'G')):2.39,
        (('G', 'C'), ('A', 'A')):4.38,
        (('U', 'C'), ('U', 'G')):5.27,
        (('U', 'G'), ('C', 'A')):0.8,
        (('A', 'A'), ('G', 'U')):3.8,
        (('C', 'U'), ('U', 'U')):4.31,
        (('C', 'C'), ('G', 'C')):5.56,
        (('G', 'C'), ('G', 'U')):2.14,
        (('U', 'U'), ('U', 'U')):0.0,
        (('U', 'C'), ('G', 'A')):6.91,
        (('C', 'U'), ('C', 'U')):8.25,
        (('C', 'G'), ('A', 'C')):2.78,
        (('G', 'U'), ('U', 'C')):3.8,
        (('A', 'G'), ('C', 'G')):4.5,
        (('G', 'G'), ('U', 'C')):10.0,
        (('U', 'G'), ('U', 'C')):4.59,
        (('C', 'A'), ('C', 'U')):4.58,
        (('G', 'A'), ('U', 'C')):1.53,
        (('G', 'G'), ('C', 'G')):10.0,
        (('A', 'U'), ('G', 'U')):2.11,
        (('G', 'U'), ('C', 'A')):4.76,
        (('A', 'C'), ('U', 'U')):5.21,
        (('A', 'G'), ('C', 'A')):6.7,
        (('C', 'A'), ('U', 'A')):2.47,
        (('A', 'G'), ('U', 'G')):6.07,
        (('C', 'G'), ('U', 'A')):0.21,
        (('U', 'G'), ('A', 'U')):2.4,
        (('G', 'G'), ('C', 'C')):10.0,
        (('C', 'G'), ('A', 'A')):4.5,
        (('C', 'U'), ('C', 'A')):5.91,
        (('G', 'U'), ('U', 'G')):4.48,
        (('G', 'C'), ('C', 'A')):2.78,
        (('A', 'A'), ('U', 'C')):0.0,
        (('G', 'G'), ('G', 'C')):10.0,
        (('A', 'U'), ('G', 'C')):0.21,
        (('C', 'A'), ('C', 'G')):2.55,
        (('U', 'C'), ('C', 'A')):5.21,
        (('U', 'G'), ('U', 'U')):2.89,
        (('G', 'G'), ('G', 'U')):10.0,
        (('C', 'A'), ('A', 'G')):5.05,
        (('A', 'C'), ('U', 'A')):2.75,
        (('A', 'U'), ('C', 'G')):0.34,
        (('U', 'G'), ('C', 'C')):5.06,
        (('A', 'C'), ('G', 'C')):2.55,
        (('G', 'U'), ('A', 'C')):0.8,
        (('C', 'A'), ('C', 'C')):4.49,
        (('A', 'U'), ('G', 'G')):10.0,
        (('A', 'A'), ('U', 'U')):6.46,
        (('G', 'G'), ('G', 'G')):0.0,
        (('A', 'C'), ('G', 'G')):10.0,
        (('A', 'C'), ('A', 'A')):4.8,
        (('A', 'A'), ('G', 'G')):10.0,
        (('G', 'C'), ('U', 'C')):3.39,
        (('U', 'C'), ('G', 'C')):3.8,
        (('U', 'A'), ('G', 'U')):2.4,
        (('G', 'A'), ('G', 'C')):3.44,
        (('G', 'G'), ('A', 'C')):10.0,
        (('A', 'A'), ('A', 'A')):2.71,
        (('C', 'C'), ('A', 'G')):8.82,
        (('C', 'G'), ('C', 'C')):5.49,
        (('U', 'A'), ('U', 'G')):2.11,
        (('U', 'G'), ('A', 'C')):4.76,
        (('C', 'U'), ('U', 'G')):6.25,
        (('C', 'C'), ('C', 'G')):5.49,
        (('G', 'U'), ('G', 'A')):4.33,
        (('G', 'C'), ('C', 'G')):0.26,
        (('A', 'G'), ('A', 'G')):2.41,
        (('A', 'A'), ('C', 'C')):8.25,
        (('A', 'G'), ('G', 'A')):2.18,
        (('C', 'G'), ('G', 'U')):2.39,
        (('A', 'C'), ('C', 'C')):5.91,
        (('A', 'G'), ('G', 'G')):10.0,
        (('U', 'U'), ('G', 'A')):6.96,
        (('U', 'C'), ('C', 'U')):6.46,
        (('A', 'U'), ('U', 'A')):0.31,
        (('G', 'A'), ('G', 'A')):2.33,
        (('G', 'C'), ('C', 'U')):3.44,
        (('G', 'A'), ('U', 'A')):3.5,
        (('U', 'G'), ('A', 'A')):6.07,
        (('C', 'U'), ('A', 'G')):8.86,
        (('U', 'U'), ('G', 'U')):5.27,
        (('A', 'U'), ('A', 'G')):3.67,
        (('G', 'C'), ('C', 'C')):5.56,
        (('A', 'U'), ('C', 'C')):5.39,
        (('U', 'C'), ('A', 'U')):3.63,
        (('C', 'A'), ('G', 'C')):2.78,
        (('C', 'G'), ('C', 'G')):0.0,
        (('C', 'C'), ('A', 'A')):9.77,
        (('G', 'C'), ('A', 'U')):0.21,
        (('A', 'U'), ('G', 'A')):3.67,
        (('U', 'U'), ('A', 'A')):8.18,
        (('U', 'C'), ('C', 'C')):4.31,
        (('C', 'U'), ('G', 'A')):8.82,
        (('C', 'A'), ('U', 'G')):0.8,
        (('C', 'U'), ('A', 'A')):9.05,
        (('C', 'A'), ('A', 'A')):6.7,
        (('C', 'C'), ('U', 'A')):5.3,
        (('G', 'U'), ('C', 'U')):4.59,
        (('C', 'A'), ('A', 'C')):4.93,
        (('C', 'A'), ('G', 'U')):4.76,
        (('G', 'C'), ('G', 'G')):10.0,
        (('U', 'A'), ('G', 'G')):10.0,
        (('G', 'A'), ('A', 'A')):3.77,
        (('U', 'U'), ('A', 'G')):6.91,
        (('C', 'G'), ('C', 'U')):3.39,
        (('A', 'G'), ('A', 'C')):4.8,
        (('U', 'U'), ('A', 'C')):5.21,
        (('U', 'A'), ('U', 'U')):3.63,
        (('G', 'U'), ('G', 'C')):2.14,
        (('C', 'U'), ('A', 'C')):4.49,
        (('C', 'A'), ('G', 'A')):5.14,
        (('C', 'G'), ('G', 'G')):10.0,
        (('C', 'C'), ('G', 'A')):8.86,
        (('G', 'A'), ('G', 'U')):4.59,
        (('U', 'A'), ('A', 'G')):3.67,
        (('U', 'U'), ('G', 'G')):10.0,
        (('U', 'A'), ('C', 'U')):3.5,
        (('C', 'C'), ('C', 'A')):4.49,
        (('C', 'U'), ('G', 'C')):5.49,
        (('C', 'U'), ('C', 'G')):5.56,
        (('A', 'A'), ('G', 'A')):2.25,
        (('C', 'U'), ('G', 'G')):10.0,
        (('C', 'A'), ('U', 'U')):2.39,
        (('A', 'G'), ('G', 'U')):4.1,
        (('G', 'U'), ('U', 'A')):2.4,
        (('U', 'G'), ('G', 'U')):4.48,
        (('C', 'A'), ('A', 'U')):2.75,
        (('C', 'G'), ('A', 'G')):3.49,
        (('G', 'G'), ('U', 'G')):10.0,
        (('A', 'U'), ('U', 'C')):3.5}


class memoize(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,profile,ref_seq,struct,*args):
      nargs = (args)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(profile,ref_seq,struct,*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()

class memoize_iso(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,ref_seq,*args):
      nargs = (args)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(ref_seq,*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()

def energy(xxx_todo_changeme6, xxx_todo_changeme7,alpha):
  #stacking energy of base pair (a,b) around base pair (a2,b2)
  (a,b) = xxx_todo_changeme6
  (a2,b2) = xxx_todo_changeme7
  E = STACKING_ENERGY[a,a2,b2,b]
  return  math.exp(-(alpha*E)/(BOLTZMANN*T))

@memoize_iso
def isostericity(ref_seq, xxx_todo_changeme, xxx_todo_changeme1, alpha):
  #isostericity of going from original base pair to (a,b)
  (i,j) = xxx_todo_changeme
  (a,b) = xxx_todo_changeme1
  if not ref_seq:
    return 1
  iso = sum(ISO[(ref[i],ref[j]),(a,b)] for ref in ref_seq)/len(ref_seq)
  #iso_start = sum(ISO[(ref[i],ref[j]),(seq[i],seq[j])] for ref in ref_seq)
  #iso = mpf(iso_mut-iso_start)/len(ref_seq)
  return  math.exp(-((1-alpha)*iso)/(BOLTZMANN*T))

@memoize
def forward(profile,ref_seq,struct, xxx_todo_changeme2, xxx_todo_changeme3,alpha):
  #alpha gives the weight energy vs isostericity
  (i,j) = xxx_todo_changeme2
  (a,b) = xxx_todo_changeme3
  result = 0.
  if i > j :
    result=1.
  else:
    k = struct[i]
    if k==-1:
      for a2 in BASES[i]:
        pro = profile[i][a2]
        f = forward(profile,ref_seq,struct,
                    (i+1,j),
                    (a2,b),
                    alpha)
        result +=  pro*f
    elif i < k <= j: #If k > j we return 0
      for a2 in BASES[i]:
        for b2 in BASES[k]:
          pro = profile[i][a2]*profile[k][b2]
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(struct)-1 or not (j==k and struct[i-1]==j+1):
            f1 = forward(profile,ref_seq,struct,
                        (i+1,k-1),
                        (a2,b2),
                        alpha)
            f2 = forward(profile,ref_seq,struct,
                         (k+1,j),
                         (b2,b),
                         alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*f1*f2*iso
          #if stack, we add energy
          else :
            f = forward(profile,ref_seq,struct,
                       (i+1,k-1),
                       (a2,b2),
                       alpha)
            e = energy((a,b),
                       (a2,b2),
                       alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*f*e*iso
  return result

@memoize
def backward(profile,ref_seq,struct, xxx_todo_changeme4, xxx_todo_changeme5,alpha):
  (i,j) = xxx_todo_changeme4
  (a,b) = xxx_todo_changeme5
  result = 0.
  if i<0:
    result=forward(profile,ref_seq,struct,
                   (j,len(struct)-1),
                   ('X','X'),
                   alpha)
  else:
    k = struct[i]
    if k==-1:
      for a2 in BASES[i]:
        pro = profile[i][a2]
        back = backward(profile,ref_seq,struct,
                        (i-1,j),
                        (a2,b),
                        alpha)
        result += pro*back
    #BP to the left
    elif k<i:
      for a2 in BASES[k]:
        for b2 in BASES[i]:
          pro = profile[k][a2]*profile[i][b2]
          back = backward(profile,ref_seq,struct,
                       (k-1,j),
                       (a2,b),
                       alpha)
          forw = forward(profile,ref_seq,struct,
                      (k+1,i-1),
                      (a2,b2),
                      alpha)
          iso = isostericity(ref_seq,
                             (k,i),
                             (a2,b2),
                             alpha)
          result += pro*back*forw*iso
    #BP to the right
    elif k>=j:
      for a2 in BASES[i]:
        for b2 in BASES[k]:
          pro = profile[i][a2]*profile[k][b2]
          #No stack
          if not (j==k and struct[i+1]==j-1):
            back = backward(profile,ref_seq,struct,
                         (i-1,k+1),
                         (a2,b2),
                         alpha)
            forw = forward(profile,ref_seq,struct,
                        (j,k-1),
                        (b,b2),
                        alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*back*forw*iso
          #Stack
          else:
            back = backward(profile,ref_seq,struct,
                         (i-1,k+1),
                         (a2,b2),
                         alpha)
            e = energy((a2,b2),
                       (a,b),
                       alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*back*e*iso
  return result

def parseStruct(dbn):
  p = []
  result = [-1 for c in dbn]
  for i in range(len(dbn)):
    c = dbn[i]
    if c=='(':
      p.append(i)
    elif c==')':
      j = p.pop()
      result[j] = i
      result[i] = j
  return result

def product_given_i(profile,ref_seq,struct,i,a,alpha):
  """Will compute the sum of boltzmann weights of structures 
  where the 'i-th' nucleotide is 'a'.
  """
  n = len(struct)
  tot = forward(profile,ref_seq,struct,(0,n-1),('X', 'X'),alpha)
  k = struct[i]
  result = mpf(0)
  if k == -1:
    pro = profile[i][a]
    result += pro*backward(profile,ref_seq,struct,(i-1,i+1),(a,a),alpha)
  elif k < i:
    for c in BASES[k]:
      pro = profile[k][c]*profile[i][a]
      f = forward(profile,ref_seq,struct,(k+1,i-1),(c,a),alpha) 
      b = backward(profile,ref_seq,struct,(k-1,i+1),(c,a),alpha)
      iso = isostericity(ref_seq,(k,i),(c,a),alpha)
      result += pro*f*b*iso
  else:
    for c in BASES[k]:
      pro = profile[i][a]*profile[k][c]
      f = forward(profile,ref_seq,struct,(i+1,k-1),(a,c),alpha) 
      b = backward(profile,ref_seq,struct,(i-1,k+1),(a,c),alpha)
      iso = isostericity(ref_seq,(i,k),(a,c),alpha)
      result += pro*f*b*iso
  return result

def random_weighted_sampling(l_samples):
  tot = sum(x[1] for x in l_samples)
  scaled_weights = [x[1]/tot for x in l_samples]
  rand_nb = random.random()
  accumulation = 0
  for i,x in enumerate(scaled_weights):
    accumulation += x
    if accumulation > rand_nb:
      return l_samples[i][0]
  return l_samples[-1][0]


def backtrack(profile,ref_seq,struct, xxx_todo_changeme8, xxx_todo_changeme9,alpha):
  #alpha gives the weight energy vs isostericity
  (i,j) = xxx_todo_changeme8
  (a,b) = xxx_todo_changeme9
  result_list = []
  #max_seq_list = []
  max_seq = ''
  if i > j :
    return '' 
  else:
    k = struct[i]
    if k==-1:
      l_samples = []
      for a2 in BASES[i]:
        pro = profile[i][a2]
        f = forward(profile,ref_seq,struct,
                    (i+1,j),
                    (a2,b),
                    alpha)
        result =  pro*f
        l_samples.append((a2,result))
      a2 = random_weighted_sampling(l_samples)
      max_seq = a2 + backtrack(profile,ref_seq,struct,(i+1,j),(a2,b),alpha)
      
    elif i < k <= j: #If k > j we return 0
      l_samples = []
      for a2 in BASES[i]:
        for b2 in BASES[k]:
          pro = profile[i][a2]*profile[k][b2]
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(struct)-1 or not (j==k and struct[i-1]==j+1):
            f1 = forward(profile,ref_seq,struct,
                        (i+1,k-1),
                        (a2,b2),
                        alpha)
            f2 = forward(profile,ref_seq,struct,
                         (k+1,j),
                         (b2,b),
                         alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result = pro*f1*f2*iso
          #if stack, we add energy
          else :
            f = forward(profile,ref_seq,struct,
                       (i+1,k-1),
                       (a2,b2),
                       alpha)
            e = energy((a,b),
                       (a2,b2),
                       alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result = pro*f*e*iso
          l_samples.append(((a2,b2),result))
      """
      for a2,b2 in result_list:
        for best_1 in backtrack(profile,ref_seq,struct,(i+1,k-1),(a2,b2),alpha):
          for best_2 in backtrack(profile,ref_seq,struct,(k+1,j),(b2,b),alpha): 
            max_seq_list.append(a2+best_1+b2+best_2)
      """
      a2,b2 = random_weighted_sampling(l_samples)
      best_1 = backtrack(profile,ref_seq,struct,(i+1,k-1),(a2,b2),alpha)
      best_2 = backtrack(profile,ref_seq,struct,(k+1,j),(b2,b),alpha)
      max_seq = a2+best_1+b2+best_2
  return max_seq

def probability_given_i(profile,ref_seq,struct,i,a,alpha):
  """Will compute the probability that the 'i-th' nucleotide
  is 'a' over all sequences at 'm' mutations from seq
  """
  n = len(struct)
  tot = forward(profile,ref_seq,struct,(0,n-1),('X', 'X'),alpha)
  result = product_given_i(profile,ref_seq,struct,i,a,alpha)
  if tot == 0:
    print("""The total partition function is 0, you might want to increase
    the number of mutations allowed""")
    sys.exit(1)
  return result/tot

def testSingleSequence(profile,ref_seq,struct,alpha):
  forward.resetCache()
  backward.resetCache()
  n = len(struct)

  i = 29
  #print "Given i=%s and m=%s the probability of sequence is:" %(i,m)
  #print "\t", probability_given_i_most_m(ref_seq,struct,i,'C',m,alpha=1.0)
  print("  Forward: \t",forward(profile,ref_seq,struct,(0,n-1),('A','G'),alpha))
  i = n-1
  res = 0
  for j in BASES[i]:
    res += profile[i][j]*backward(profile,ref_seq,struct,(i-1,i+1),(j,j),alpha)
  print("  Backward of nuc %s:\t" % i, res)

def test():
  
  seq = "UAUAUAUGUAAUACAACAAACAAUAAAAGGUG"
  dbn = "................................"
  ref_seq = seq
  alpha = 0.2 

  ref_seq = ref_seq*1
  seq = seq*1
  dbn = dbn*1
  struct = parseStruct(dbn)
  profile = tuple({'A':0.25,
                   'C':0.25,
                   'G':0.25,
                   'U':0.25} for x in range(len(seq)))



  testSingleSequence(profile,ref_seq,struct,alpha)

def parse_fasta(file_name):
  #Sequences in the MFE and the target secondary structure
  seq = []
  struct = []
  with open(file_name) as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      if all(x in 'AUGC' for x in line):
        seq.append(line)
        continue
      if all(x in '(.)' for x in line):
        struct.append(line)
        continue
  return seq, parseStruct(struct[0])

def parse_profile(file_name):
  #Parse a profile in the format of the output of RNApyro
  profile = []
  with open(file_name) as f:
    for l in f:
      l = l.strip().split()
      prob = {'A':mpf(l[0]),
              'C':mpf(l[1]),
              'G':mpf(l[2]),
              'U':mpf(l[3])}
      profile.append(prob)
  return tuple(profile)

def equiprob_profile(n):
      prob = {'A':0.25,
              'C':0.25,
              'G':0.25,
              'U':0.25}
      profile = tuple(prob for _ in range(n))
      return profile

def all_probabilities(profile,ref_seq, stuct, alpha):
  n = len(struct)
  results = []
  for i in range(n):
    results.append([])
    for a in IUPACBASES['N']:
      if a not in BASES[i]:
        results[-1].append(0)
      else:
        results[-1].append(
          probability_given_i(profile,ref_seq,struct,i,a,alpha))
  return results

def display_all_probabilities(results):
  for x in results:
    print(x[0],x[1],x[2],x[3])

def gc_content(sequence,structure=None):
  if not structure:
    nb = len([x for x in sequence if x in 'GC'])
    return float(nb)/len(sequence)
  else:
    nb = len([x for i,x in enumerate(sequence) if x in 'GC' and 
             structure[i] > -1])
    nb_bp = len([x for x in structure if x > -1])
    return float(nb) / nb_bp

def update_profile(profile,max_bound,min_bound,increase=True):
  """Given a profile, will increase (resp. decrease) by half the probability
  of GC content at all positions. Will equaly distribute the remaining between
  AU
  """
  new_profile = []
  for position in profile:
    if increase:
      up_G = (max_bound-position['G'])/2
      up_C = (max_bound-position['C'])/2
    else:
      up_G = -(position['G']-min_bound)/2
      up_C = -(position['C']-min_bound) /2
    to_remove = up_G+up_C
    new_profile.append(
      {'G':position['G']+up_G,
       'C':position['C']+up_C,
       'A':position['A']-to_remove/2})
    #make sure sums to 1
    new_profile[-1]['U'] = 1.0 - sum(new_profile[-1][x] for x in 'GCA')
  return new_profile
     
def sample_gc_target(profile,ref_seq,struct,alpha,nb_gc_sample,gc_target,
                     file_gc_data,f_gc_only_structure,
                     max_err=0.1,sample_before_update=1000):
  """Will sample the required number of sequences, inside the gc_target

  """
  n = len(struct)
  lower_gc = gc_target-max_err
  upper_gc = gc_target+max_err
  max_bound = 1.0
  min_bound = 0.0
  l_all_sample = []
  l_correct_gc = []

  if file_gc_data:
    l_contents = []
    f = open(file_gc_data,'w')

  while len(l_correct_gc) < nb_gc_sample:
    l_all_sample.append([backtrack(profile,ref_seq,struct,(0,n-1),('',''),alpha)
                         for _ in range(sample_before_update)])
    over_under = 0 
    if file_gc_data:
      l_contents[:] = []
    for sample in l_all_sample[-1]:
      if f_gc_only_structure:
        content = gc_content(sample,structure=struct)
      else:
        content = gc_content(sample)
      if file_gc_data:
        l_contents.append(content)
      over_under += (content - gc_target)
      if lower_gc <= content <= upper_gc: 
        l_correct_gc.append(sample)

    if file_gc_data:
      f.write('%s\n' % profile[0]['C']) 
      for c in l_contents[:-1]:
        f.write('%s\t' % c)
      f.write('%s\n' % l_contents[-1])
    if over_under < 0:
      min_bound = profile[0]['C']
      profile = update_profile(profile,max_bound,min_bound) 
      forward.clear()
    elif over_under > 0:
      max_bound = profile[0]['C']
      profile = update_profile(profile,max_bound,min_bound,increase=False)
      forward.clear()

  if file_gc_data:
    f.close()

  return l_correct_gc

def sub_seq_structure(sequence,struct):
  return ''.join(x for i,x in enumerate(sequence) if struct[i] > -1)

def diversity_seq(l_sequence,struct):
  n = len(l_sequence)
  n_set = float(len(set(l_sequence)))
  n_struct_set = float(len(set(sub_seq_structure(x,struct)
                               for x in l_sequence)))
  return n,n_set/n,n_struct_set/n

def help():
  print("""
  Required:
    -d <file_path> 
      A file containing the target secondary structure and an optional MSA
    -a <float> 
      The value of alpha, between 0 and 1.
      1 takes only into account the secondary structure, 0 only the MSA

  Optional:
    -m <int> 
      The max penality for an invalid base pair, -1 for infinity
    -b <int> 
      print 'n' stochasticly backtracked sequences
    -no_profile <> 
      Doesn't output a profile
    -s_gc <target_gc> <nb_samples> 
      Sampling sequences with a 0<=target_gc<=1 and a given
      number of samples 
    -gc_sec_struct <>
      Only the nucleotides with an interaction in the secondary structure
      will be considered for the GC content
    -gc_max_err <float>
      Max error from GC target allowed in sample, default 0.1
    -gc_data <file_name>
      Create a file with the given name, will print to it a first line
      the weight of 'C' and on the next the list of sampled GC content
    -t <float>
      The temperature (default 310.5K)
    -c <IUPAC sequence>
      An IUPAC sequence to constrain the outputed sequences
    -p <file_path> 
      A data file with  starting RNA profile (i.e. every lines contains
       nucleotides probability in order: 'ACGU')

  e.g.
    python IncaRNAtion -d data.txt -a 0.5 -m 20
    python IncaRNAtion -d data.txt -a 1 -m 20 -no_profile -s_gc 0.5 100
    python IncaRNAtion -d data.txt -a 0.5 -m 20 -b 5 -no_profile
    """)

if __name__ == "__main__":
  opts = sys.argv
  l = len(opts)
  i = 1

  #Action Flags!!
  f_no_profile = False#Don't do profile
  f_backtrack = False#Do backtrack
  f_sample_gc = False#Sample with given gc content
  f_gc_only_structure = False#Only take into account GC with interactions
  f_gc_max_error = False#Custom max error for gc content, default 0.1
  file_gc_data=None#If a file given, print to it for each iter. weight C, GC content

  while i < l:
    cmd = opts[i]
    if not cmd.startswith('-'):
      i += 1
    else:
      if cmd == '-d':  
        file_name = opts[i+1]
        i += 2
        if not os.path.isfile(file_name):
          help()
          sys.exit(1)
      elif cmd == '-p':
        profile_path = opts[i+1]
        i += 2
        if not os.path.isfile(profile_path):
          help()
          sys.exit(1)
      #alpha
      elif cmd == '-a':
        alpha = opts[i+1]
        i += 2
        try:
          alpha = float(alpha)
        except ValueError:
          help()
          sys.exit(1)
        if not 0 <= alpha <= 1:
          help()
          sys.exit(1)
      #Max BP penality (def inf.)
      elif cmd == '-m':
        z = opts[i+1]
        i += 2
        try:
          z = float(z)
          if z >= 0:
            for x in STACKING_ENERGY:
              if STACKING_ENERGY[x] == sys.maxsize:
                STACKING_ENERGY[x] = z
        except ValueError:
          pass
      #Don't output profile
      elif cmd == '-no_profile':
        f_no_profile = True
        i += 1
      #Backtrack N optimal sequences
      elif cmd == '-b':
        f_backtrack = True
        nb_backtrack = opts[i+1]
        i += 2
        try:
          nb_backtrack = int(nb_backtrack)
        except ValueError:
          print("Error -b")
          help()
          sys.exit(1)
        if nb_backtrack < 1:
          print("Error -b")
          help()
          sys.exit(1)
      #Sampling Sequences given GC target
      elif cmd == '-s_gc':
        f_sample_gc = True
        gc_target = opts[i+1]
        nb_gc_sample = opts[i+2]
        i += 3
        try:
          nb_gc_sample = int(nb_gc_sample)
          gc_target = float(gc_target)
        except ValueError:
          print("Error s_gc")
          help()
          sys.exit(1)
        if not (0 <= gc_target <= 1) or not nb_gc_sample > 0:
          print("Error s_gc")
          help()
          sys.exit(1)
      #Only check sec struct pos for GC content
      elif cmd == '-gc_sec_struct':
        f_gc_only_structure = True
        i += 1
      #Temperature of the system
      elif cmd == '-t':
        T = opts[i+1]
        i += 2
        try:
          T = float(T)
        except ValueError:
          help()
          sys.exit(1)
      #Max error from GC_Target allowed
      elif cmd == '-gc_max_err':
        f_gc_max_error = True
        max_error = opts[i+1]
        i += 2
        try:
          max_error = float(max_error)
        except ValueError:
          print("Error gc_max_err")
          help()
          sys.exit(1)
        if not 0 <= max_error <= 1:
          print("Error gc_max_err")
          help()
          sys.exit(1)
      #If want data from GC sampling
      elif cmd == '-gc_data':
        file_gc_data = opts[i+1]
        i += 2
      elif cmd == '-c':
        iupac = opts[i+1]
        i+=2 
      else:
        print("Unrecognized arg")
        help()
        sys.exit(1)



  try:
    ref_seq,struct = parse_fasta(file_name)
  except NameError:
    help()
    sys.exit(1)

  n = len(struct)

  #iupac
  try:
    if len(iupac) != n:
      print("IUPAC code too short")
      help()
      sys.exit(1)
    for x in iupac.upper():
      try:
        BASES.append(IUPACBASES[x])
      except KeyError:
        print("Unrecognized IUPAC symbol")
        help()
        sys.exit(1)
  except NameError:
    for _ in range(n):
      BASES.append(IUPACBASES['N'])



  try:
    profile = parse_profile(profile_path)
  except NameError:
    profile = equiprob_profile(n)

  #Action!

  if not f_no_profile:
    results = all_probabilities(profile,ref_seq,struct,alpha)
    display_all_probabilities(results)

  if f_backtrack:
    for _ in range(nb_backtrack):
      res = backtrack(profile,ref_seq,struct,(0,n-1),('',''),alpha)
      print(res, gc_content(res,struct))

  if f_sample_gc:
    if f_gc_max_error:
      res = sample_gc_target(profile,ref_seq,struct,alpha,nb_gc_sample,
                             gc_target,file_gc_data,f_gc_only_structure,
                             max_err=max_error)
    else:
      res = sample_gc_target(profile,ref_seq,struct,alpha,
                             nb_gc_sample,gc_target,file_gc_data,
                             f_gc_only_structure)
      print('\n'.join(res))

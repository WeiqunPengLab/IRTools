#!/usr/bin/env python
'''
This python3 module was converted from python2.7 code using 2to3
'''

#import built-in modules
import re
import sys
from random import shuffle
from heapq import nlargest
#import third-party modules

#changes to the paths

#changing history to this module



__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="3.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def S_diff(lst):
    '''Given a list of int or float, calculate S_diff and S_point'''
    
    S_avg = sum(lst) / len(lst)
    S_dist = [i-S_avg for i in lst] #distance to average
    S_cum=[]    #list of cumulative sum
    S_cum.append(0)
    for i in range(0,len(S_dist)):
        S_cum.append(S_cum[i] + S_dist[i])
    return [nlargest(1,list(range(0,len(S_cum))),key=lambda i: S_cum[i]),(max(S_cum) - min(S_cum))]
    #return the index of maximum_diff index, and maximum_diff

def bootstrap(lst,obs,rep=1000):
    '''Given a list of int or float (lst) and an observation value(obs). calcualte the chance (pvalue) 
    of getting this observation through bootstrapping.'''
    
    shuffled_diff=[]
    count=0
    tmp=lst
    for i in range(0,rep):
        shuffle(tmp)
        shuffled_diff.append(S_diff(tmp))
    
    for i in sorted(shuffled_diff):
        if (i>=obs):
            count += 1 
    if count/rep <0.5:
        return count/rep
    else:
        return 1- count/rep

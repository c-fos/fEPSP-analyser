'''
Created on 11.03.2012

@author: pilat
'''
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from numpy import array

def neuroCheck(maxdiff,length,front,frequ):
    importr("nnet")
    ro.r['load'](".RneuroModel2")
    tmpMatrix=ro.r['matrix'](ro.FloatVector([maxdiff,length/frequ,front*frequ]),nrow = 1,ncol = 3)
    result=ro.r['predict'](ro.r['neuroModel'],tmpMatrix)
    #print(([maxdiff,length/frequ,front*frequ],result))
    return array(result)[0][0]
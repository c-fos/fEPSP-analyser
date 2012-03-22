'''
Created on 11.03.2012

@author: pilat
'''
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from numpy import array

def neuroCheck(length,front,frequ):
    importr("nnet")
    ro.r['load'](".RneuroModel")
    tmpMatrix=ro.r['matrix'](ro.FloatVector([length/frequ,front*frequ]),nrow = 1,ncol = 2)
    result=ro.r['predict'](ro.r['neuroModel'],tmpMatrix)
    return array(result)[0][0]
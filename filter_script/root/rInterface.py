'''
Created on 11.03.2012

@author: pilat
'''
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.rlike import container
from numpy import array
def neuroCheck(maxtomin,delay,length,angle1,angle2,frequ):
    angles=angle1/angle2
    importr("nnet")
    ro.r['load']("/home/pilat/workspace/fEPSP-analyser/filter_script/root/.RneuroModel4")
    tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length/frequ,delay,maxtomin,angle1,angles],['length','delay','maxtomin','angle1','angles'])))
    result=ro.r['predict'](ro.r['neuroModel'],tmpMatrix)
    #print(([maxtomin,length/frequ,angles],result))
    return array(result)[0][0]
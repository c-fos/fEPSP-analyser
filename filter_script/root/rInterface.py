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
    print("fibrecheck 0")
    ro.r('.libPaths("/home/pilat/R/i686-pc-linux-gnu-library/2.15/")')
    nnet=importr("nnet")
    print(1)
    ro.r.load("/home/pilat/workspace/fEPSP-analyser/filter_script/root/.RneuroModel3")
    print(2)
    tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length/frequ,delay,maxtomin,angle1,angles],['length','delay','maxtomin','angle1','angles'])))
    print(3)
    result=ro.r.predict(ro.r['neuroModel'],tmpMatrix)
    #print(([maxtomin,length/frequ,angles],result))
    return array(result)[0][0]
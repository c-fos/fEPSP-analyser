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
    ro.r('.libPaths("/home/pilat/R/i686-pc-linux-gnu-library/2.15/")')
    nnet=importr("nnet")
    ro.r.load("/home/pilat/workspace/fEPSP-analyser/filter_script/root/.RneuroModel8")
    tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length,delay,maxtomin,angle1,angles],['length','delay','maxtomin','angle1','angles'])))
    #tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length,maxtomin,angle1,angles],['length','maxtomin','angle1','angles'])))
    result=ro.r.predict(ro.r['neuroModel'],tmpMatrix)
    return array(result)[0][0]

def stimNeuroCheck(length,ptp,base_ptp,std,median,mean,base_mean,base_median,base_std,diff_median,base_diff_median,sampleStdDiff,samplePtpDiff):
    median_baseMedian=median-base_median
    if base_ptp==0:
        base_ptp=10**-16
    ptp_base_ptp=ptp/base_ptp
    med_ptp=abs(median_baseMedian)/base_ptp
    ro.r('.libPaths("/home/pilat/R/i686-pc-linux-gnu-library/2.15/")')
    nnet=importr("nnet")
    ro.r.load("/home/pilat/workspace/fEPSP-analyser/filter_script/root/.RStimNeuroModel11")
    tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length,ptp,base_ptp,std,median,mean,base_mean,base_median,base_std,diff_median,base_diff_median,median_baseMedian,med_ptp,ptp_base_ptp,sampleStdDiff,samplePtpDiff],['length','ptp','base_ptp','std','median','mean','base_mean','base_median','base_std','diff_median','base_diff_median','median_baseMedian','med_ptp','ptp_base_ptp','sampleStdDiff','samplePtpDiff'])))
    #tmpMatrix=ro.r['data.frame'](ro.DataFrame(container.TaggedList([length,maxtomin,angle1,angles],['length','maxtomin','angle1','angles'])))
    result=ro.r.predict(ro.r['neuroModel'],tmpMatrix)
    return array(result)[0][0]
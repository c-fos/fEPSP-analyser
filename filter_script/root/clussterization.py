'''
Created on 15.02.2012

@author: pilat
'''
from hcluster import fclusterdata
from numpy import array, unique, where, zeros
import sys

def clusterization(fromObject,spikeDict,stimuli,debug,isClusterOn):
    if len(spikeDict)>1:
        dictValues=array(spikeDict.values())
        listOfSpikes=[]
        for i in dictValues:
            tmpObject=getattr(fromObject,i)
            listOfSpikes.append([tmpObject.spikeMin])
        ndarrayOfSpikes=array(listOfSpikes)
        rightClasterOrder=zeros(ndarrayOfSpikes.size,dtype=int)
        if isClusterOn==1:
            if debug==1:
                print("clusterization is on")
            try:
                clusteredSpikes=fclusterdata(ndarrayOfSpikes,1.1,depth=4,method='average')
            except:
                print "Unexpected error in fclusterdata:", sys.exc_info()        
            clusterNumbers=unique(clusteredSpikes)
            for i in clusterNumbers:
                k=where(rightClasterOrder==0)[0][0]
                mask=clusteredSpikes==clusteredSpikes[k]
                rightClasterOrder+=mask*i
            if debug==1:
                print((clusteredSpikes,rightClasterOrder))
            return(rightClasterOrder)
        else:
            if debug==1:
                print("clusterization is off")
            ndarrayOfSpikes.shape=(1,ndarrayOfSpikes.size)
            for i in range(len(stimuli[0])):
                rightClasterOrder[ndarrayOfSpikes[0]>=stimuli[0][i]]=i+1
            return (rightClasterOrder)
    else:
        return(array([1]))
        
def clusterAnalyser(fromObject,spikeDict,clusters):
    dictValues=array(spikeDict.values())
    for i in range(len(dictValues)):
        try:
            tmpObject=getattr(fromObject,dictValues[i])
            tmpObject.responsNumber=int(clusters[i])
            k=where(clusters==int(clusters[i]))[0][0]
            tmpObject.spikeNumber=int(i)-k
        except:
            print "Unexpected error in cluster analysis:", sys.exc_info()
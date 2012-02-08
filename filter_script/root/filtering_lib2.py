#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''
#library for filtering
from mimetypes import guess_type
from numpy import math, zeros, arange, loadtxt, std, fromfile, array, int16, unique, where, mean#arange,zeros,array,math,fromfile,int16, r_,convolve, hanning, hamming, bartlett, blackman, ones,loadtxt,roll,diff,sign,nonzero
#from externalFunctions import iswt, smooth, extrema
#import sys
from matplotlib.transforms import Bbox
from matplotlib.path import Path
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
from externalFunctions import *
from dbAccess import *
from objects import *
from hcluster import fclusterdata

class dataSample:
    def __init__(self,filename,dbobject,arguments):        
        #variables
        self.spikeObjectList=[]
        self.spikeDict=dict()
        self.responseList=[]
        self.responseDict=dict()
        self.fileName = str(filename)
        self.mysql_writer=dbobject
        self.arguments=arguments
       
        try:
            self.argReading()#0
        except:
            print("argReading() complete with error")
        if self.debug==1:
            print(self.fileName)
        ##
        try:
            self.dataLoading()
        except:
            print("dataLoading() complete with error")
        try:
            self.tresholdCreating() 
        except:
            print("tresholdCreating() complete with error")
        try:
            self.cleanData=self.cutStimuli(self.data)
        except:
            print("cutStimuli() complete with error")  
        try: 
            self.snr=self.snrFinding(self.cleanData,self.defaultFrame)
        except:
            print("snrFinding() complete with error")  
        try:
            self.mainLevelFinding()
        except:
            print("mainLevelFinding() complete with error")
        try:
            self.filtering()
        except:
            print("filtering() complete with error")
            
        try:
            self.spikeFinding()
        except:
            print("spikeFinding() complete with error")
        try:
            self.clusters=self.clusterization()
        except:
            print("clusterization() complete with error") 
        try:
            self.clusterAnalyser()
        except:
            print("clusterAnalyser() complete with error")   
        try:
            self.responsMatrix=self.responsLength()
        except:
            print("responsLength() complete with error")
        try:
            self.responsAnalysis()
        except:
            print("responsAnalysis() complete with error")
            
        """     
        try:
            self.spikesIntoResponses()
        except:
            print("spikesIntoResponses() complete with error")

        """
        try:
            self.plotData()
        except:
            print("plotData() complete with error") 
        try:
            self.writeData()
        except:
            print("writeData() complete with error")
        
    #0-reading command line arguments


    def argReading(self):
        self.wavelet='sym7'#'bior3.5
        self.frequency = int(self.arguments[2])
        self.destination = str(self.arguments[3])
        self.coeffTreshold = int(self.arguments[4])
        self.debug = int(self.arguments[6])
        self.write = int(self.arguments[7])
        self.isClusterOn = int(self.arguments[8])
                      
    #1-data Loading

    def dataLoading(self):
        if guess_type(self.fileName)[0] == 'text/plain' or guess_type(self.fileName)[0] == 'chemical/x-mopac-input':
            self.data = loadtxt(self.fileName)
            self.data = self.dataFitting(self.data,self.frequency)
        else:
            self.data=fromfile(self.fileName,int16)
            self.data=self.dataFitting(self.data,self.frequency)
    
    #1.1-data fitting


    def dataFitting(self,data,frequency):
        tmp=2**(math.ceil(math.log(len(data),2)))
        delay=2*frequency/1000#asume that first 2 msec doesn`t contain any signal 
        newData=zeros(tmp)
        deltaLen=tmp-len(data)
        for j in range(int(math.ceil(deltaLen/delay))):
            for k in range(delay):
                newData[delay*j+k]=data[k]
        for i in range(len(data)):
            newData[deltaLen+i]=data[i]
        return newData


    def tresholdCreating(self):
        msec=self.frequency/1000 # 1 msec =  frequency/1000 points        
        self.defaultFrame=5*msec #frame size for mean() and std() finding must depend on frequency. assume it equeal to 5msec
        self.stimulyDuration=int(1.5*msec) # treshold for stimuli filtering ~20points==2msec==2*self.msec
        self.level=int((math.log(100.0/self.frequency,0.5)-1)) #wavelet decomposition level. level 6 to 10kHz signal.
        self.baseFrequency=150 #we must separate the levels of wavelet decomposition wich contains most part of the signal
        highNoiseFrequency=14000.0
        self.highNoiseLevel=int((math.log(highNoiseFrequency/self.frequency,0.5)))
        self.localDelay=1*msec #time delay before spike(filtering of local responses)
        

    def snrFinding(self,data,frameSize):
        minSD=data[:frameSize].std()
        maxSD=max(data)-min(data)
        snr=int(maxSD/minSD)
        return snr


    def stdFinder(self,data,frameSize,mean=False,maxStd=False):
        dataSample=[]
        minSD=200
        maxSD=0    
        i=frameSize
        for j in range(int(len(data)/frameSize)):
            dataSample+=[[i,j,data[j*i:i*(j+1)].mean(),data[j*i:i*(j+1)].std()]]#1/3
            dataSample+=[[i,j,data[j*i+j/3.0:i*(j+1)+j/3.0].mean(),data[j*i+j/3.0:i*(j+1)+j/3.0].std()]]#2/3
            dataSample+=[[i,j,data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].mean(),data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].std()]]#3/3
        if len(dataSample)>0:
            for k in range(len(dataSample)):
                if dataSample[k][3]< minSD and dataSample[k][3]!=0:
                    minSD=dataSample[k][3]
                    index=k
                if dataSample[k][3]>maxSD and dataSample[k][3]!=0:
                    maxSD=dataSample[k][3]
            if mean==True:
                return minSD,dataSample[index][2]
            else:
                if maxStd==True:
                    if self.debug==1:
                        print((maxSD,len(dataSample)))
                    return maxSD
                else:
                    return minSD
        else:
            if self.debug==1:
                print("#")
            if mean==True:
                return data.std(),data.mean()
            else:
                return data.std()        


    def findStimuli(self,data):
        pwr=pywt.swt(data, 'db2', 3)
        self.pwr2=array(pwr[2][1])
        treshold=abs(max(self.pwr2))/4
        self.pwr2[self.pwr2<treshold]=0
        self.pwr5=smooth(self.pwr2,10)
        self.pwr5[self.pwr5>0]=treshold 
        self.dpwr=where((diff(self.pwr5)>treshold/2)==True)[0]
        self.dpwrM=where((diff(self.pwr5)<-treshold/2)==True)[0]
        stimList=[[],[]]
        for i in range(len(self.dpwr)):
            start=self.dpwr[i]
            length=(self.dpwrM[i]-self.dpwr[i])*3
            stop=start+length
            stimList[0]+=[start]
            stimList[1]+=[stop]
        self.stimuli=stimList

    def cutStimuli(self,data):
        self.findStimuli(data)
        if len(self.stimuli[0])!=0:
            for i in range(len(self.stimuli[0])):
                pathValue=data[self.stimuli[0][i]-(self.stimuli[1][i]-self.stimuli[0][i]):self.stimuli[0][i]].mean()
                data[self.stimuli[0][i]:self.stimuli[1][i]]=pathValue
                
            return data
        else:
            return data


    def mainLevelFinding(self):
        if self.debug==1:
            print(self.snr)
        self.mainLevel=int((math.log((self.baseFrequency-4.0*(self.baseFrequency/self.snr))/self.frequency,0.5))-2)#
        if self.debug==1:
            print(self.mainLevel)
    

    def filtering(self):
        self.coeffs=pywt.swt(self.cleanData, self.wavelet, level=self.mainLevel+1)
        for i in range(len(self.coeffs)):
            cA, cD = self.coeffs[i]
            if i>=(len(self.coeffs)-self.highNoiseLevel):
                cD=zeros(len(cA))
                if self.debug==1:
                    print(("noisLevel",i))
            else:
                minSD=self.stdFinder(cD,self.defaultFrame)
                cD=pywt.thresholding.soft(cD,minSD*(self.coeffTreshold+i**4))
            self.coeffs[i]=cA, cD
        self.result=iswt(self.coeffs,self.wavelet)
   

    def spikeFinding(self):
        resultData=self.result
        start=self.defaultFrame
        stop=-self.defaultFrame
        minimum,minimumValue = extrema(resultData[start:stop],_max = False, _min = True, strict = False, withend = True)
        maximum,maximumValue = extrema(resultData[start:stop],_max = True, _min = False, strict = False, withend = True)
        std=self.stdFinder(self.cleanData,self.defaultFrame)
        SD=std*2*(self.coeffTreshold-5*self.coeffTreshold/self.snr)#? maybe we must add the snr check?
        spikePoints=[]
        if minimum[0]<maximum[0]:
            minimum=minimum[1:]
            minimumValue=minimumValue[1:]
        for i in range(len(minimum)):
            tmpMaximum1=0
            for j in range(len(maximum)):
                if maximum[j]<minimum[i]:
                    tmpMaximum1=j
                if maximum[j]>minimum[i]:
                    tmpMaximum2=j
                    break
            if maximumValue[tmpMaximum1]-minimumValue[i]>SD and maximumValue[tmpMaximum2]-minimumValue[i]>SD:
                spikePoints.append([start+maximum[tmpMaximum1],start+minimum[i],start+maximum[tmpMaximum2]])
        for i in range(len(spikePoints)):
            ampl=round(resultData[spikePoints[i][0]]-resultData[spikePoints[i][1]]+(resultData[spikePoints[i][2]]-resultData[spikePoints[i][0]])/(spikePoints[i][2]-spikePoints[i][0])*(spikePoints[i][1]-spikePoints[i][0]),1)
            if ampl>SD:
                index=len(self.spikeDict)
                self.spikeDict[index]="n"+str(i)
                setattr(self,self.spikeDict[index],Spike(self.frequency))
                tmpObject=getattr(self,self.spikeDict[index])                
                tmpObject.responseStart=start
                tmpObject.responseEnd=stop
                tmpObject.spikeNumber=i
                tmpObject.spikeMax1=spikePoints[i][0]
                tmpObject.spikeMax1Val=self.result[spikePoints[i][0]]
                tmpObject.spikeMin=spikePoints[i][1]
                tmpObject.spikeMinVal=self.result[spikePoints[i][1]]
                tmpObject.spikeMax2=spikePoints[i][2]
                tmpObject.spikeMax2Val=self.result[spikePoints[i][2]]
                tmpObject.spikeAmpl=ampl
                tmpObject.calculate()
        print((self.fileName,self.spikeDict))
     
        
    def clusterization(self):
        if len(self.spikeDict)>1:
            dictValues=array(self.spikeDict.values())
            listOfSpikes=[]
            for i in dictValues:
                tmpObject=getattr(self,i)
                listOfSpikes.append([tmpObject.spikeMin])
            ndarrayOfSpikes=array(listOfSpikes)
            rightClasterOrder=zeros(ndarrayOfSpikes.size,dtype=int)
            if self.isClusterOn==1:
                print("clusterization is on")
                try:
                    clusteredSpikes=fclusterdata(ndarrayOfSpikes,1.1,depth=4,method='average')
                except:
                    print("fclusterdata error")                
                clusterNumbers=unique(clusteredSpikes)
                for i in clusterNumbers:
                    k=where(rightClasterOrder==0)[0][0]
                    mask=clusteredSpikes==clusteredSpikes[k]
                    rightClasterOrder+=mask*i
                if self.debug==1:
                    print((clusteredSpikes,rightClasterOrder))
                return(rightClasterOrder)
            else:
                print("clusterization is off")
                ndarrayOfSpikes.shape=(1,ndarrayOfSpikes.size)
                for i in range(len(self.stimuli[0])):
                    rightClasterOrder[ndarrayOfSpikes[0]>=self.stimuli[0][i]]=i+1
                return (rightClasterOrder)
        else:
            return(array([1]))
        

    def clusterAnalyser(self):
        dictValues=array(self.spikeDict.values())
        clusters=self.clusters
        for i in range(len(dictValues)):
            try:
                tmpObject=getattr(self,dictValues[i])
                tmpObject.responsNumber=int(clusters[i])
                k=where(clusters==int(clusters[i]))[0][0]
                tmpObject.spikeNumber=int(i)-k
                #print((i,clusters[i],tmpObject.responsNumber,tmpObject.spikeNumber))
            except:
                pass
            
    def responsLength(self):
        responsMatrix=zeros((max(self.clusters),2),dtype=int)#[[start1,stop1],[start2,stop2]]
        length=len(self.result)
        smallFrame=self.defaultFrame/5
        print(len(unique(self.clusters)),len(self.clusters))
        for i in unique(self.clusters):
            print("0")
            firstSpike=self.spikeDict.values()[list(self.clusters).index(i)]
            try:
                lastSpike=self.spikeDict.values()[list(self.clusters).index(i+1)-1]
            except:
                lastSpike=self.spikeDict.values()[-1]
            print("1")
            if self.isClusterOn==1:
                tmpObject=getattr(self,firstSpike)
                firstMax=tmpObject.spikeMax1
                k=firstMax
                baseLevel=self.result[k-self.defaultFrame:k-smallFrame].mean()
                std2=self.result[k-self.defaultFrame:k-smallFrame].std()
                while(abs(self.result[k-smallFrame:k].mean()-baseLevel)>std2/2 or self.result[k-smallFrame:k].std()>std2/2):
                    k-=smallFrame/2
                start=k
                tmpObject=getattr(self,lastSpike)
                lastMax=tmpObject.spikeMax2
                k=lastMax
                while((abs(self.result[k:k+smallFrame*3].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*3].std()>std2/2) and k<length-smallFrame*4):
                    k+=smallFrame
                stop=k
                responsMatrix[i-1]=start,stop
            else:
                print("2")
                start=self.stimuli[0][i-1]
                baseLevel=self.result[start-smallFrame*2:start].mean()
                tmpObject=getattr(self,lastSpike)
                lastMax=tmpObject.spikeMax2
                k=lastMax
                std2=self.result[start-smallFrame*2:start].std()
                try: 
                    while((abs(self.result[k:k+smallFrame*2].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*2].std()>std2/2) and (k<length-smallFrame*4 and k<self.stimuli[0][i])):
                        k+=smallFrame
                        print((i,k,k+smallFrame*2,self.stimuli[0][i]))
                except:
                    print("3")
                    k=lastMax
                    std2=self.result[start-smallFrame*2:start].std()
                    while((abs(self.result[k:k+smallFrame*2].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*2].std()>std2/2) and k<length-smallFrame*4):
                        k+=smallFrame
                        
                stop=k
                print(("Stop--",stop))
                responsMatrix[i-1]=start,stop
                print("4")
        return responsMatrix
    
    def responsAnalysis(self):
        rMatrix=self.responsMatrix
        for i in unique(self.clusters):
            index=len(self.responseDict)
            self.responseDict[index]="r"+str(i)
            setattr(self,self.responseDict[index],Response())
            tmpObject=getattr(self,self.responseDict[index])
            tmpObject.responseStart=rMatrix[i-1][0]
            tmpObject.responseEnd=rMatrix[i-1][1]
            tmpObject.response_top=max(self.result[rMatrix[i-1][0]:rMatrix[i-1][1]])
            tmpObject.response_bottom=min(self.result[rMatrix[i-1][0]:rMatrix[i-1][1]])
            tmpObject.baselevel=self.result[rMatrix[i-1][0]-self.defaultFrame/2:rMatrix[i-1][0]].mean()
            tmpObject.responsNumber=i
            tmpObject.vpsp=round(tmpObject.response_top-tmpObject.baselevel,1)
            tmpObject.spikes=array(self.spikeDict.values())[self.clusters==i]
        print(self.responseDict)
        

    def plotData(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.pwr5,'r')
        ax.plot(self.data,'y')
        ax.plot(self.result,'b')
        ax.grid(color='k', linestyle='-', linewidth=0.4)
        try:
            for i in self.responseDict.values():
                print(i)
                tmpObject=getattr(self,i)
                rect = Rectangle((tmpObject.responseStart, tmpObject.response_bottom), tmpObject.responseEnd-tmpObject.responseStart, tmpObject.response_top-tmpObject.response_bottom, facecolor="#aaaaaa", alpha=0.3)
                ax.add_patch(rect)   
                #ax.axhline(y=tmpObject.baselevel,color='g')
                ax.text(tmpObject.responseStart,tmpObject.response_top+20, "VPSP="+str(tmpObject.vpsp), fontsize=12, va='bottom')
                try:
                    for j in tmpObject.spikes:
                        print(j)
                        tmpObject2=getattr(self,j)
                        tex = str((tmpObject2.responsNumber,tmpObject2.spikeNumber,tmpObject2.spikeAmpl))
                        ax.plot(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin],'or')
                        ax.plot(tmpObject2.spikeMax1,self.result[tmpObject2.spikeMax1],'og')
                        ax.plot(tmpObject2.spikeMax2,self.result[tmpObject2.spikeMax2],'og')
                        ax.vlines(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin], self.result[tmpObject2.spikeMin]+tmpObject2.spikeAmpl, color='k', linestyles='dashed')
                        ax.text(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin]-15, tex, fontsize=12, va='bottom')
                except:
                    pass
        except:
            pass
        for i in range(len(self.stimuli[0])):
            ax.axvline(x=self.stimuli[0][i],color='g')
        plt.savefig(self.fileName+"_graph.png")
        #plt.show()
        plt.close()# very important to stop memory leak
        del fig
        
    def writeData(self):
        if self.write:
            for i in self.responseDict.values():
                tmpObject=getattr(self,i)
                self.mysql_writer.dbWriteResponse(tmpObject)
                for j in tmpObject.spikes:
                    tmpObject2=getattr(self,j)
                    self.mysql_writer.dbWriteSpike(tmpObject2)
                    #del tmpObject2
            #del tmpObject
            
            
            
            
            
            
            
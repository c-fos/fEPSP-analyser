#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011
master 
@author: pilat
'''
#library for filtering
from mimetypes import guess_type
from numpy import math, zeros, arange, loadtxt, fromfile, array, int16#arange,zeros,array,math,fromfile,int16, r_,convolve, hanning, hamming, bartlett, blackman, ones,loadtxt,roll,diff,sign,nonzero
#from externalFunctions import iswt, smooth, extrema
#import sys
from matplotlib.transforms import Bbox
from matplotlib.path import Path
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
from externalFunctions import *
from dbAccess import *
from objects import *

class dataSample:
    def __init__(self,filename,dbobject,spikeListObject,arguments):        
        #variables
        self.spikeObjectList=[]
        self.spikeDict=dict()
        self.fileName = str(filename)
        self.mysql_writer=dbobject
        self.spikeList=spikeListObject
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
            self.responsesFind()
        except:
            print("responsesFind() complete with error")      
        try:
            self.spikesIntoResponses()
        except:
            print("spikesIntoResponses() complete with error")
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
            dataSample=dataSample+[[i,j,data[j*i:i*(j+1)].mean(),data[j*i:i*(j+1)].std()]]#1/3
            dataSample=dataSample+[[i,j,data[j*i+j/3.0:i*(j+1)+j/3.0].mean(),data[j*i+j/3.0:i*(j+1)+j/3.0].std()]]#2/3
            dataSample=dataSample+[[i,j,data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].mean(),data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].std()]]#3/3
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


    def findStarts(self,data):
        pwr=pywt.swt(data, 'db2', 2, start_level=0)
        pwr2=abs(array(pwr[0][1]))
        treshold=int(max(pwr2)/2)
        pwr3=pywt.thresholding.greater(pwr2, treshold)
        pwr4=smooth(pwr3,self.stimulyDuration)
        if self.debug==1:
            print((type(extrema),max(pwr4)))
        self.stimuli=extrema(pwr4,_min=False)[0]
        for i in range(len(self.stimuli)):
            if self.stimuli[i] > self.defaultFrame and self.stimuli[i] < len(data) - self.defaultFrame:
                self.stimuli[i]=self.stimuli[i]-self.stimulyDuration/4
            else:
                pass
        self.stimuli=self.stimuli[self.stimuli.__gt__(self.defaultFrame)*self.stimuli.__lt__(len(data)-self.defaultFrame)]


    def cutStimuli(self,data):
        self.findStarts(data)
        if len(self.stimuli)!=0:
            tmpData=zeros(len(data))
            for i in range(len(data)):
                tmpData[i]=data[i]
            for i in range(len(self.stimuli)):
                start=self.stimuli[i]-self.stimulyDuration/2
                summ=0    
                for k in range(self.stimulyDuration):
                    summ=summ+data[start+k-self.stimulyDuration]
                averege=summ/self.stimulyDuration
                for j in range(self.stimulyDuration):
                    if start+j< len(tmpData)-1:
                        tmpData[start+j]=averege
        return tmpData


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

    
    def responsesFind(self):
        self.responsesList=[]
        for i in range(len(self.stimuli)):
            if len(self.result[self.stimuli[i]:]) > self.defaultFrame:
                if i<len(self.stimuli)-1 and self.stimuli[i+1]< self.stimuli[i]+self.responsLength(self.result[self.stimuli[i]:]):
                    self.responsesList=self.responsesList+[[self.stimuli[i],self.stimuli[i+1]]]
                else:
                    self.responsesList=self.responsesList+[[self.stimuli[i],self.stimuli[i]+self.responsLength(self.result[self.stimuli[i]:])]]
            


    def responsLength(self,sample):
        dataStd,dataMean=self.stdFinder(sample,self.defaultFrame,True)
        dataSample=[]
        data=sample
        i=self.defaultFrame
        if len(sample)<self.defaultFrame*7:
            sampleLength=len(sample)
        else:
            for j in arange(2,int(len(data)/i),1):
                dataSample=dataSample+[[i,j,data[j*i:i*(j+1)].mean(),data[j*i:i*(j+1)].std(),0]]#1/3
                dataSample=dataSample+[[i,j,data[j*i+j/3.0:i*(j+1)+j/3.0].mean(),data[j*i+j/3.0:i*(j+1)+j/3.0].std(),0.33]]#2/3
                dataSample=dataSample+[[i,j,data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].mean(),data[j*i+j*2/3.0:i*(j+1)+j*2/3.0].std(),0.66]]#3/3
            sampleLength=len(sample)
            if dataSample>0:        
                for k in range(len(dataSample)):
                    if abs(dataSample[k][2]-dataMean)<=dataStd*5 and dataSample[k][3]<dataStd*3:
                        sampleLength=(dataSample[k][1]+dataSample[k][4])*dataSample[k][0]
                        break
        if self.debug==1:
            print(("sampleLength",sampleLength,len(sample),dataStd,dataMean))
        return sampleLength
    


    def spikesIntoResponses(self):        
        self.resultList=[]
        if self.debug==1:
            print("len(self.responsesList)",len(self.responsesList))
        for i in range(len(self.responsesList)):
            self.extremums(i,self.responsesList[i][0]+self.localDelay,self.responsesList[i][1])


    def extremums(self,responsNumber,start,stop):
        minimum,minimumValue = extrema(self.result[start:stop],_max = False, _min = True, strict = False, withend = True)
        maximum,maximumValue = extrema(self.result[start:stop],_max = True, _min = False, strict = False, withend = True)
        std,baseLevel=self.stdFinder(self.cleanData,self.defaultFrame,True)
        SD=std*2*(self.coeffTreshold-5*self.coeffTreshold/self.snr)#? maybe we must add the snr check?
        #extremum filtering
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
            ampl=self.result[spikePoints[i][0]]-self.result[spikePoints[i][1]]+(self.result[spikePoints[i][2]]-self.result[spikePoints[i][0]])/(spikePoints[i][2]-spikePoints[i][0])*(spikePoints[i][1]-spikePoints[i][0])
            width=spikePoints[i][2]-spikePoints[i][0]
            response_top=max(self.result[start:stop])
            response_bottom=min(self.result[start:stop])     
            if width>self.stimulyDuration and ampl>SD and spikePoints[i][2]-start>self.localDelay:
                index=len(self.spikeDict)
                self.spikeDict[index]="r"+str(responsNumber)+"n"+str(i)
                setattr(self.spikeList,self.spikeDict[index],Spike(self.frequency))
                tmpObject=getattr(self.spikeList,self.spikeDict[index])
                tmpObject.responsNumber=responsNumber
                tmpObject.responseStart=start
                tmpObject.responseEnd=stop
                tmpObject.spikeNumber=i
                tmpObject.allSpikes=len(spikePoints)
                tmpObject.spikeMax1=spikePoints[i][0]
                tmpObject.spikeMax1Val=self.result[spikePoints[i][0]]
                tmpObject.spikeMin=spikePoints[i][1]
                tmpObject.spikeMinVal=self.result[spikePoints[i][1]]
                tmpObject.spikeMax2=spikePoints[i][2]
                tmpObject.spikeMax2Val=self.result[spikePoints[i][2]]
                tmpObject.vpsp=response_top-baseLevel
                tmpObject.response_top=response_top
                tmpObject.response_bottom=response_bottom
                tmpObject.baselevel=baseLevel
                tmpObject.spikeAmpl=ampl
                tmpObject.calculate()
        print((self.fileName,self.spikeDict))


    def plotData(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.cleanData,'r')
        ax.plot(self.data,'y')
        ax.plot(self.result,'b')
       # ax.plot(self.spikeLevel,'r')
        ax.grid(color='k', linestyle='-', linewidth=0.4)
        
        
        for i in range(len(self.spikeDict)):
            tmpObject=getattr(self.spikeList,self.spikeDict[i])
            tex = str((tmpObject.responsNumber,tmpObject.spikeNumber,tmpObject.spikeAmpl))
            ax.plot(tmpObject.spikeMin,self.result[tmpObject.spikeMin],'or')
            ax.plot(tmpObject.spikeMax1,self.result[tmpObject.spikeMax1],'og')
            ax.plot(tmpObject.spikeMax2,self.result[tmpObject.spikeMax2],'og')
            ax.vlines(tmpObject.spikeMin,self.result[tmpObject.spikeMin], self.result[tmpObject.spikeMin]+tmpObject.spikeAmpl, color='k', linestyles='dashed')
            ax.text(tmpObject.spikeMin,self.result[tmpObject.spikeMin]-15, tex, fontsize=12, va='bottom')
            rect = Rectangle((tmpObject.responseStart, tmpObject.response_bottom), tmpObject.responseEnd, tmpObject.response_top-tmpObject.response_bottom, facecolor="#aaaaaa", alpha=0.1)
            ax.add_patch(rect)   
            ax.axhline(y=tmpObject.baselevel,color='g')
            ax.text(tmpObject.responseStart,tmpObject.response_top+20, "VPSP="+str(tmpObject.vpsp), fontsize=12, va='bottom')
        for i in range(len(self.stimuli)):
            ax.axvline(x=self.stimuli[i],color='g')
        plt.savefig(self.fileName+"_graph.png")
        #plt.show()
        plt.close()# very important to stop memory leak
        del fig
        
    def writeData(self):
        tmpResponse=-1                
        for i in range(len(self.spikeDict)):
            tmpObject=getattr(self.spikeList,self.spikeDict[i])#self.spikeDict[0])
            if self.write:    
                self.mysql_writer.variables_local(tmpObject)
                if tmpResponse!=self.spikeDict[i].split('r')[1].split('n')[0]:
                    tmpResponse=self.spikeDict[i].split('r')[1].split('n')[0]
                    self.mysql_writer.dbWriteResponse()            
                self.mysql_writer.dbWriteSpike()
            del tmpObject
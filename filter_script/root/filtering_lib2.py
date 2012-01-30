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
        """
        try:
            self.writeData()
        except:
            print("writeData() complete with error")
        """
        
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
        pwr=pywt.swt(data, 'db2', 2, start_level=0)
        self.pwr2=abs(array(pwr[0][1]))
        treshold=pwr[0][1].std()*10
        pwr5=pywt.thresholding.greater(self.pwr2, treshold)
        self.dpwr=diff(pwr5)>=treshold*3
        print(("dpwr",where(self.dpwr==True)))
        treshold=int(max(self.pwr2)/2)
        pwr3=pywt.thresholding.greater(self.pwr2, treshold)
        self.stimulMask=self.pwr2>=treshold
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
        self.findStimuli(data)
        if len(self.stimuli)!=0:
            """
            smallFrame=int(self.defaultFrame/50)
            print((self.stimuli,smallFrame))
            data=array(data)
            for i in self.stimuli:
                iSample=data[i-self.defaultFrame/2:i]
                print((iSample,i-self.defaultFrame/2,i,len(data[i-self.defaultFrame/2:i])))
                iStd,iMean=1,100
                try:
                    iMean=mean(iSample)
                except:
                    pass
                try:
                    iStd=iSample.std()
                except:
                    pass
                print(("*",iStd,iMean))
                iMask=(iSample<(iMean+iStd))-(iSample<(iMean-iStd/2))#?
                print(("**",iMask))
                baselevel=mean(iSample[iMask])
                std1=std(iSample[iMask])
                k=smallFrame
                sample=data[(i-k)-smallFrame/2:(i-k)+smallFrame/2]
                print((std1,iStd,iMean,baselevel,sample))
                while(sample.std()>std1 or abs(sample.mean()-baselevel)>std1*3):
                    k+=smallFrame/2
                    sample=data[(i-k)-smallFrame/2:(i-k)+smallFrame/2]
                start=i-k
                k=i-smallFrame
                while(abs(data[k]-baselevel)>std1/4):
                    k+=1
                stop=k
                """
            mask=self.stimulMask
            patchValues=data[:-100]*mask[100:]
            data=data*(-mask)
            #data[100:]+=patchValues
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
            try:
                clusteredSpikes=fclusterdata(ndarrayOfSpikes,1.1,depth=4,method='average')
            except:
                print("fclusterdata error")
            rightClasterOrder=zeros(clusteredSpikes.size,dtype=int)
            clusterNumbers=unique(clusteredSpikes)
            for i in clusterNumbers:
                k=where(rightClasterOrder==0)[0][0]
                mask=clusteredSpikes==clusteredSpikes[k]
                rightClasterOrder+=mask*i
            if self.debug==1:
                print((clusteredSpikes,rightClasterOrder))
            return(rightClasterOrder)
        else:
            return(array([1]))


    def clusterAnalyser(self):
        dictValues=array(self.spikeDict.values())
        clusters=self.clusters
        #print((clusters,dictValues))
        for i in range(len(dictValues)):
            #print(i)
            try:
                tmpObject=getattr(self,dictValues[i])
                tmpObject.responsNumber=int(clusters[i])
                k=where(clusters==int(clusters[i]))[0][0]
                tmpObject.spikeNumber=int(i)-k
                #print((i,clusters[i],tmpObject.responsNumber,tmpObject.spikeNumber))
            except:
                pass
            
    def responsLength(self):
        responsMatrix=zeros((len(unique(self.clusters)),2))#[[start1,stop1],[start2,stop2]]
        length=len(self.result)
        smallFrame=self.defaultFrame/4
        #print((responsMatrix,length))
        for i in unique(self.clusters):
            firstSpike=self.spikeDict.values()[list(self.clusters).index(i)]
            try:
                lastSpike=self.spikeDict.values()[list(self.clusters).index(i+1)-1]
            except:
                lastSpike=self.spikeDict.values()[-1]
            #print("*")
            tmpObject=getattr(self,firstSpike)
            firstMax=tmpObject.spikeMax1
            k=firstMax
            baseLevel=self.result[k-self.defaultFrame:k].mean()
            std2=self.result[k-self.defaultFrame:k].std()
            #print((baseLevel,std2))
            while(abs(self.result[k-smallFrame:k].mean()-baseLevel)>std2/2 or self.result[k-smallFrame:k].std()>std2/2):
                #print("start")
                k-=smallFrame/2
            start=k
            tmpObject=getattr(self,lastSpike)
            lastMax=tmpObject.spikeMax2
            k=lastMax
            #print(("#",baseLevel,std2))
            while((abs(self.result[k:k+smallFrame*4].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*4].std()>std2/2) and k<length):
                k+=smallFrame
                #print("stop")
            stop=k
            #print((start,firstMax,lastMax,stop))
            responsMatrix[i-1]=start,stop
        return responsMatrix
    
    def responsAnalysis(self):
        rMatrix=self.responsMatrix
        #print((len(unique(self.clusters)),self.clusters))
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
            tmpObject.vpsp=tmpObject.response_top-tmpObject.baselevel
            tmpObject.spikes=array(self.spikeDict.values())[self.clusters==i]
            #print(tmpObject.spikes)
        print(self.responseDict)
        

    def plotData(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[self.dpwr],'r')
        ax.plot(self.data,'y')
        ax.plot(self.result,'b')
        #ax.plot(self.spikeLevel,'r')
        ax.grid(color='k', linestyle='-', linewidth=0.4)
        try:
            for i in self.responseDict.values():
                print(i)
                tmpObject=getattr(self,i)
                rect = Rectangle((tmpObject.responseStart, tmpObject.response_bottom), tmpObject.responseEnd-tmpObject.responseStart, tmpObject.response_top-tmpObject.response_bottom, facecolor="#aaaaaa", alpha=0.3)
                ax.add_patch(rect)   
                ax.axhline(y=tmpObject.baselevel,color='g')
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
        for i in range(len(self.stimuli)):
            ax.axvline(x=self.stimuli[i],color='g')
        #plt.savefig(self.fileName+"_graph.png")
        plt.show()
        #plt.close()# very important to stop memory leak
        del fig
        
    def writeData(self):
        tmpResponse=-1                
        for i in range(len(self.spikeDict)):
            tmpObject=getattr(self,self.spikeDict[i])#self.spikeDict[0])
            if self.write:    
                self.mysql_writer.variables_local(tmpObject)
                if tmpResponse!=self.spikeDict[i].split('r')[1].split('n')[0]:
                    tmpResponse=self.spikeDict[i].split('r')[1].split('n')[0]
                    self.mysql_writer.dbWriteResponse()            
                self.mysql_writer.dbWriteSpike()
            del tmpObject
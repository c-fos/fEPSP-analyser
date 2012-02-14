#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''
#library for filtering
import sys
from mimetypes import guess_type
from numpy import zeros, append, math, ones, diff, loadtxt, fromfile, array, int16, unique, where,float16,float32
import pywt
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from externalFunctions import iswt,extrema,smooth
from objects import Spike,Response
from hcluster import fclusterdata
from scipy.interpolate import spline, Rbf
from scipy import polyval, polyfit

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
        self.errorState=0
    
    def dataProcessing(self):       
        try:
            self.argReading()#0
        except:
            print "argReading() error:", sys.exc_info()
            self.errorState=1
            ##raise
        if self.debug==1:
            print(self.fileName)
        ##
        try:
            self.dataLoading()
        except:
            print "dataLoading() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.tresholdCreating() 
        except:
            print "tresholdCreating() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.cleanData=self.cutStimuli(self.data)
        except:
            print "cutStimuli() error:", sys.exc_info()
            self.errorState=1
            ##raise 
        try: 
            self.snr=self.snrFinding(self.cleanData,self.defaultFrame)
        except:
            print "snrFinding() error:", sys.exc_info()
            self.errorState=1
            ##raise 
        try:
            self.mainLevelFinding()
        except:
            print "mainLevelFinding() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.filtering()
        except:
            print "filtering() error:", sys.exc_info()
            self.errorState=1
            ##raise
            
        try:
            self.spikeFinding()
        except:
            print "spikeFinding() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.clusters=self.clusterization()
        except:
            print "clusterization() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.clusterAnalyser()
        except:
            print "clusterAnalyser() error:", sys.exc_info()
            self.errorState=1
            ##raise   
        try:
            self.responsMatrix=self.responsLength()
        except:
            print "responsLength() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.responsAnalysis()
        except:
            print "responsAnalysis() error:", sys.exc_info()
            self.errorState=1
            ##raise
        try:
            self.plotData()
        except:
            print "plotData() error:", sys.exc_info()
            self.errorState=1
            ###raise
        if self.mysql_writer!="pass":
            try:
                self.writeData()
            except:
                print "writeData() error:", sys.exc_info()
                self.errorState=1
                ###raise
        
        
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
            self.data = loadtxt(self.fileName,dtype='float32')
            self.data = self.dataFitting(self.data,self.frequency)
        else:
            self.data=fromfile(self.fileName,int16)
            self.data=self.dataFitting(self.data,self.frequency)
    
    #1.1-data fitting


    def dataFitting(self,data,frequency):
        tmp=2**(math.ceil(math.log(len(data),2)))
        delay=2*frequency/1000#asume that first 2 msec doesn`t contain any signal 
        newData=zeros(tmp, dtype='float32')
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
        self.stimulyDuration=int(1*msec) # treshold for stimuli filtering ~20points==2msec==2*self.msec
        self.level=int((math.log(100.0/self.frequency,0.5)-1)) #wavelet decomposition level. level 6 to 10kHz signal.
        self.baseFrequency=150 #we must separate the levels of wavelet decomposition wich contains most part of the signal
        highNoiseFrequency=14000.0
        self.highNoiseLevel=int((math.log(highNoiseFrequency/self.frequency,0.5)))
        self.localDelay=1*msec #time delay before spike(filtering of local responses)
        

    def snrFinding(self,data,frameSize):
        minSD=data[:frameSize].std()
        maxSD=max(data)-min(data)
        snr=float16(maxSD/minSD)
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
        pwr=pywt.swt(data, 'haar', 2)
        pwr2=array(pwr[0][1])
        treshold=max(abs(pwr2))/4
        pwr2[pwr2<treshold]=0
        pwr2[pwr2>0]=treshold 
        dpwr=where((diff(pwr2)>treshold/2)==True)[0]
        dpwrMask=ones(len(dpwr),dtype='bool')
        for i in range(len(dpwr)-1):
            if dpwr[i+1]-dpwr[i]<self.stimulyDuration:
                dpwrMask[i+1]=0
        dpwr=dpwr[dpwrMask]            
        stimList=[[],[]]
        for i in range(len(dpwr)):
            start=dpwr[i]
            length=self.stimulyDuration/10
            baseline=data[start-30:start].mean()
            baseStd=data[start-30:start].std()
            tmpStop=start+length
            if tmpStop+self.defaultFrame+5>len(data):
                tmpStop=len(data)-self.defaultFrame-5
            stimMean=data[start:tmpStop].mean()
            if self.debug==1:
                print((baseline,baseStd,stimMean,"baseline,baseStd,stimMean"))
            sample=smooth(data[tmpStop-5:tmpStop+self.defaultFrame+5],10)
            firstArray=abs(diff(sample))<=baseStd
            secondArray=abs(sample[1:]-baseline)<baseStd
            thirdArray=abs(sample[1:]-baseline)<abs(baseline-stimMean)/10
            try:
                shift=where((firstArray+secondArray)*thirdArray==True)[0]
                if len(shift)>0:
                    realStop=tmpStop+where((firstArray+secondArray)*thirdArray==True)[0][0]
                else:
                    realStop=tmpStop
            except:
                print "Unexpected error in finding of stimuli end:", sys.exc_info()
                realStop=tmpStop
            if self.debug==1:
                print((start,tmpStop,realStop,"start,tmpStop and stop of stimule"))
            stimList[0]+=[start]
            stimList[1]+=[realStop]
        self.stimuli=stimList

    def cutStimuli(self,data):
        try:
            self.findStimuli(data)
        except:
            print "Unexpected error in findStimuli", sys.exc_info()
            ###raise
        if len(self.stimuli[0])!=0:
            if self.debug==1:
                print("cut stimuli")
            for i in range(len(self.stimuli[0])):
                pathValue=data[self.stimuli[0][i]-(self.stimuli[1][i]-self.stimuli[0][i]):self.stimuli[0][i]].mean()
                data[self.stimuli[0][i]:self.stimuli[1][i]]=pathValue    
            return data
        else:
            return data


    def mainLevelFinding(self):
        if self.debug==1:
            print(self.snr)
        self.mainLevel=int((math.log((self.baseFrequency-4*(self.baseFrequency/self.snr))/self.frequency,0.5))-2)#
        if self.debug==1:
            print(self.mainLevel)
    

    def filtering(self):
        self.coeffs=pywt.swt(self.cleanData, self.wavelet, level=self.mainLevel+1)
        for i in range(len(self.coeffs)):
            cA, cD = self.coeffs[i]
            if i>=(len(self.coeffs)-self.highNoiseLevel):
                cD=zeros(len(cA),dtype='float32')
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
        if self.debug==1:
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
                if self.debug==1:
                    print("clusterization is on")
                try:
                    clusteredSpikes=fclusterdata(ndarrayOfSpikes,1.1,depth=4,method='average')
                except:
                    print "Unexpected error in fclusterdata:", sys.exc_info()
                    ###raise           
                clusterNumbers=unique(clusteredSpikes)
                for i in clusterNumbers:
                    k=where(rightClasterOrder==0)[0][0]
                    mask=clusteredSpikes==clusteredSpikes[k]
                    rightClasterOrder+=mask*i
                if self.debug==1:
                    print((clusteredSpikes,rightClasterOrder))
                return(rightClasterOrder)
            else:
                if self.debug==1:
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
            except:
                print "Unexpected error in cluster analisys:", sys.exc_info()
                ###raise
            
    def responsLength(self):
        responsMatrix=zeros((max(self.clusters),2),dtype=int)#[[start1,stop1],[start2,stop2]]
        length=len(self.result)
        smallFrame=self.defaultFrame/5
        if self.debug==1:
            print(len(unique(self.clusters)),len(self.clusters))
        for i in unique(self.clusters):
            firstSpike=self.spikeDict.values()[list(self.clusters).index(i)]
            try:
                lastSpike=self.spikeDict.values()[list(self.clusters).index(i+1)-1]
            except:
                lastSpike=self.spikeDict.values()[-1]
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
                start=self.stimuli[0][i-1]
                baseLevel=self.result[start-smallFrame*2:start].mean()
                tmpObject=getattr(self,lastSpike)
                lastMax=tmpObject.spikeMax2
                k=lastMax
                std2=self.result[start-smallFrame*2:start].std()
                try: 
                    while((abs(self.result[k:k+smallFrame*2].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*2].std()>std2/2) and (k<length-smallFrame*4 and k<self.stimuli[0][i])):
                        k+=smallFrame
                        if self.debug==1:
                            print((i,k,k+smallFrame*2,self.stimuli[0][i]))
                except:
                    k=lastMax
                    std2=self.result[start-smallFrame*2:start].std()
                    while((abs(self.result[k:k+smallFrame*2].mean()-baseLevel)>std2/2 or self.result[k:k+smallFrame*2].std()>std2/2) and k<length-smallFrame*4):
                        k+=smallFrame
                        
                stop=k
                try:
                    sampleLen=stop-lastMax
                    if sampleLen>0:
                        ar=polyfit(array(range(sampleLen)),self.result[lastMax:stop],2)
                        sample=polyval(ar,array(range(sampleLen)))
                        shift=where(diff(sample)>=-0.001)[0]
                        if len(shift)>0:
                            realStop=lastMax+shift[0]
                        else:
                            realStop=lastMax
                        responsMatrix[i-1]=start,realStop
                    else:
                        responsMatrix[i-1]=start,stop
                except:
                    responsMatrix[i-1]=start,stop
                    print "Unexpected error in response length finding:", sys.exc_info()
                    ###raise
        return responsMatrix
    
    def epspReconstructor(self,tmpObject):
        tmpObject2=getattr(self,tmpObject.spikes[-1])
        #before lastMax of lastSpike
        sample1=self.result[tmpObject.responseStart:tmpObject2.spikeMax2]     
        mask=zeros(tmpObject.responseEnd-tmpObject.responseStart,dtype='bool')
        tmpObject2=getattr(self,tmpObject.spikes[0])
        mask[0]=1
        mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
        mask[tmpObject2.spikeMax1-tmpObject.responseStart]=1
        if len(tmpObject.spikes)>=2:
            try:
                tmpObject2=getattr(self,tmpObject.spikes[1])
                mask[tmpObject2.spikeMax1-tmpObject.responseStart]=1
                mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
            except:
                print "there is no second spike?:", sys.exc_info()
                ###raise
        #after lastMax of lastSpike (curve reconstruction)
        tmpObject2=getattr(self,tmpObject.spikes[-1])
        sample2=self.result[tmpObject2.spikeMax2:tmpObject.responseEnd]
        #mask2=zeros(len(sample2),dtype='bool')
        sample2Points=array(range(len(sample2)))
        ar1=Rbf(sample2Points[sample2Points%200==0],sample2[sample2Points%200==0],smooth=0.01)#,function='gaussian')
        xr1=ar1(sample2Points)
        step=len(sample2Points)/8
        for i in range(8):
            mask[tmpObject2.spikeMax2-tmpObject.responseStart+step*i]=1
        mask[-1]=1
        #append
        #mask=append(mask1,mask2)
        sample=append(sample1,xr1)
        if self.debug==1:
            print((len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd]),"len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd])"))
        timePoints=where(mask==True)[0]
        values=sample[mask]
        sample3=array(range(len(sample)))
        try:           
            #sp2 = Rbf(timePoints,int16(values),smooth=0.0001)
            sp3 = Rbf(timePoints,int16(values),smooth=1,function='thin_plate')
            y4=sp3(sample3)
            try:
                front,back = self.epspAnaliser(y4)
                self.epsp=append(self.epsp,[sample3+tmpObject.responseStart,y4],axis=1)
                if self.debug==1:
                    print((front,back,"front,back"))
                return front,back
            except:
                print "Unexpected error in epspAnaliser:", sys.exc_info()
                ###raise
                return 0,0            
        except:
            print "Unexpected error in reconstruction:", sys.exc_info()
            ###raise
            return 0,0
        
    
    def epspAnaliser(self,y):
        maxvalue=max(y)
        maxPoint=where(y==maxvalue)[0]
        minvalue1=min(y[:maxPoint])
        minvalue2=min(y[maxPoint:])
        ampl1=maxvalue-minvalue1
        ampl2=maxvalue-minvalue2
        firstPoint=where(y[:maxPoint]>minvalue1+ampl1*0.2)[0][0]
        firstValue=y[firstPoint]
        secondPoint=where(y[:maxPoint]>minvalue1+ampl1*0.8)[0][0]
        secondValue=y[secondPoint]
        thirdPoint=maxPoint+where(y[maxPoint:]<minvalue2+ampl2*0.8)[0][0]
        thirdValue=y[thirdPoint]
        fourthPoint=maxPoint+where(y[maxPoint:]<minvalue2+ampl2*0.2)[0][0]
        fourthValue=y[fourthPoint]
        front=(secondValue-firstValue)/(secondPoint-firstPoint)
        back=(fourthValue-thirdValue)/(fourthPoint-thirdPoint)
        return float32(front), float32(back[0])
                        
 
    
    def responsAnalysis(self):
        rMatrix=self.responsMatrix
        self.epsp=array([[],[]])
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
            try:
                tmpObject.epspFront,tmpObject.epspBack = self.epspReconstructor(tmpObject)
                
            except:
                print "Unexpected error wile response %s reconstruction:" % i, sys.exc_info()
                ###raise
        print(self.responseDict)
        

    def plotData(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data,'y')
        ax.plot(self.result,'b')
        ax.plot(self.epsp[0],self.epsp[1],'r')
        ax.grid(color='k', linestyle='-', linewidth=0.4)
        try:
            for i in self.responseDict.values():
                tmpObject=getattr(self,i)
                rect = Rectangle((tmpObject.responseStart, tmpObject.response_bottom), tmpObject.responseEnd-tmpObject.responseStart, tmpObject.response_top-tmpObject.response_bottom, facecolor="#aaaaaa", alpha=0.3)
                ax.add_patch(rect)
                ax.text(tmpObject.responseStart,tmpObject.response_top+20, "VPSP="+str(tmpObject.vpsp), fontsize=12, va='bottom')
                try:
                    for j in tmpObject.spikes:
                        tmpObject2=getattr(self,j)
                        tex = str((tmpObject2.responsNumber,tmpObject2.spikeNumber,tmpObject2.spikeAmpl))
                        ax.plot(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin],'or')
                        ax.plot(tmpObject2.spikeMax1,self.result[tmpObject2.spikeMax1],'og')
                        ax.plot(tmpObject2.spikeMax2,self.result[tmpObject2.spikeMax2],'og')
                        ax.vlines(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin], self.result[tmpObject2.spikeMin]+tmpObject2.spikeAmpl, color='k', linestyles='dashed')
                        ax.text(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin]-15, tex, fontsize=12, va='bottom')
                except:
                    print "Unexpected error wile spike ploating:", sys.exc_info()
                    ###raise
        except:
            print "Unexpected error wile ploating:", sys.exc_info()
            self.errorState=1
            ##raise
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
            
            
            
            
            
            
            
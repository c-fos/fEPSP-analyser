#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''
#library for filtering
import sys
from mimetypes import guess_type
from numpy import zeros, asmatrix, append, math, empty, sqrt, histogram, ones,ptp, diff, loadtxt, fromfile, array, int16, unique, where,float16,float32
import pywt
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from externalFunctions import iswt,extrema,smooth
from objects import Spike,Response
from scipy.interpolate import Rbf
from scipy import polyval, polyfit
from clussterization import clusterization, clusterAnalyser
from checkResult import resultAnalysis

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
        self.hardError=0
        self.softError=0
    
    def dataProcessing(self):       
        try:
            self.argReading()
        except:
            print "argReading() error:", sys.exc_info()
            self.hardError=1
        if self.debug==1:
            print(self.fileName)
        try:
            self.dataLoading()
        except:
            print "dataLoading() error:", sys.exc_info()
            self.hardError=1
        try:
            self.tresholdCreating() 
        except:
            print "tresholdCreating() error:", sys.exc_info()
            self.hardError=1
        try:
            self.cleanData=self.cutStimuli(self.data)
        except:
            print "cutStimuli() error:", sys.exc_info()
            self.hardError=1
        try: 
            self.snr=self.snrFinding(self.cleanData[self.deltaLen+self.stimulyDuration:],self.defaultFrame)
        except:
            print "snrFinding() error:", sys.exc_info()
            self.hardError=1
        try:
            self.mainLevelFinding()
        except:
            print "mainLevelFinding() error:", sys.exc_info()
            self.hardError=1
        try:
            self.filtering()
        except:
            print "filtering() error:", sys.exc_info()
            self.hardError=1
        try:
            self.spikeFinding()
        except:
            print "spikeFinding() error:", sys.exc_info()
            self.hardError=1
        try:
            self.clusters=clusterization(self,self.spikeDict,self.stimuli,self.debug,self.isClusterOn)
            #clusterization(fromObject,spikeDict,stimuli,debug,isClusterOn)
        except:
            print "clusterization() error:", sys.exc_info()
            self.hardError=1
        try:
            clusterAnalyser(self,self.spikeDict,self.clusters)
            #clusterAnalyser(fromObject,spikeDict,clusters)
        except:
            print "clusterAnalyser() error:", sys.exc_info()
            self.hardError=1
        try:
            self.responsMatrix=self.responsLength()
        except:
            print "responsLength() error:", sys.exc_info()
            self.hardError=1
        try:
            self.responsAnalysis()
        except:
            print "responsAnalysis() error:", sys.exc_info()
            self.hardError=1
        try:
            self.plotData()
        except:
            print "plotData() error:", sys.exc_info()
            self.hardError=1
        if self.mysql_writer!="pass":
            try:
                self.writeData()
                
            except:
                print "writeData() error:", sys.exc_info()
                self.hardError=1
            if self.hardError!=0 or self.softError!=0:
                try:
                    self.mysql_writer.dbWriteError(self.softError,self.hardError)
                except:
                    print "dbWriteError error:", sys.exc_info()
                    self.hardError=1
        try:
            resultAnalysis(self)
        except:
            print "resultAnalysis error:", sys.exc_info()
            self.hardError=1        


    def argReading(self):
        self.wavelet='sym5'#'bior3.5
        self.frequency = int(self.arguments[2])
        self.destination = str(self.arguments[3])
        self.coeffTreshold = int(self.arguments[4])
        self.debug = int(self.arguments[6])
        self.write = int(self.arguments[7])
        self.isClusterOn = int(self.arguments[8])

    def dataLoading(self):
        if guess_type(self.fileName)[0] == 'text/plain' or guess_type(self.fileName)[0] == 'chemical/x-mopac-input':
            data = loadtxt(self.fileName,dtype='float32')
            self.data = self.dataFitting(data,self.frequency)
        else:
            data=fromfile(self.fileName,int16)
            self.data=self.dataFitting(data,self.frequency)
            
    def histMean(self,sample):
        unique1=unique(sample)
        if len(unique1)>1:
            dataHist=histogram(sample,bins=unique1)
            histMax=dataHist[0].max()
            meanTmp = dataHist[1][where(dataHist[0]==histMax)[0]]
            return meanTmp[0]
        else:
            return unique1
    
    def dataFitting(self,data,frequency):
        tmp=2**(math.ceil(math.log(len(data),2)))
        self.deltaLen=tmp-len(data)
        meanTmp=self.histMean(data)
        newData=empty(tmp, dtype='float32')
        newData.fill(meanTmp)
        newData[self.deltaLen:]=data
        return newData

    def tresholdCreating(self):
        msec=self.frequency/1000 # 1 msec =  frequency/1000 points        
        self.defaultFrame=5*msec #frame size for mean() and std() finding must depend on frequency. assume it equeal to 5msec
        self.stimulyDuration=int(1*msec) # treshold for stimuli filtering ~20points==2msec==2*self.msec
        self.level=int((math.log(100.0/self.frequency,0.5)-1)) #wavelet decomposition level. level 6 to 10kHz signal.
        self.baseFrequency=300 #we must separate the levels of wavelet decomposition wich contains most part of the signal
        highNoiseFrequency=14000.0
        self.highNoiseLevel=int((math.log(highNoiseFrequency/self.frequency,0.5)))
        self.localDelay=1*msec #time delay before spike(filtering of local responses)
        
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
            
    def snrFinding(self,data,frameSize):
        minSD=self.stdFinder(data,frameSize)
        maxSD=ptp(data)
        snr=float16(maxSD/minSD)
        self.signalPtp=maxSD
        if self.debug==1:
            print((minSD,maxSD,snr,"minSD,maxSD,snr in snrFinding function"))
        return snr


      


    def findStimuli(self,data):
        pwr=pywt.swt(data, 'haar', 2)
        pwr2=abs(array(pwr[0][1]))
        treshold=max(pwr2)/2
        pwr2[pwr2<treshold]=0
        pwr2[pwr2>0]=treshold 
        dpwr=where((diff(pwr2)>treshold/2)==True)[0]
        dpwrMask=ones(len(dpwr),dtype='bool')
        for i in range(len(dpwr)-1):
            if dpwr[i+1]-dpwr[i]<self.stimulyDuration:
                dpwrMask[i+1]=0
        dpwr=dpwr[dpwrMask]
        if self.debug==1:
                print("number of finded stimuls - %s" % len(dpwr))          
        stimList=[[],[]]
        for i in range(len(dpwr)):
            start=dpwr[i]-10
            length=self.stimulyDuration/4
            baseline=self.histMean(data[start-30:start])
            baseStd=data[start-30:start].std()
            tmpStop=start+length
            if tmpStop+self.defaultFrame+5>len(data):
                tmpStop=len(data)-self.defaultFrame-5
            stimMean=data[start+10:tmpStop].mean()
            if self.debug==1:
                print((baseline,baseStd,stimMean,"baseline,baseStd,stimMean"))
            sample=smooth(data[tmpStop-5:tmpStop+self.defaultFrame+5],10)
            firstArray=abs(diff(sample))<=0.1
            secondArray=abs(sample[1:]-baseline)<baseStd/2
            thirdArray=abs(sample[1:]-baseline)<abs(baseline-stimMean)/20
            try:
                shift=where((firstArray+secondArray)*thirdArray==True)[0]
                if len(shift)>0:
                    realStop=tmpStop+shift[0]
                else:
                    realStop=tmpStop
                    print("can`t find stimulum end =(")
                    self.softError=1
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
        if len(self.stimuli[0])!=0:
            processedData=zeros(len(data),dtype='int')
            processedData+=data
            if self.debug==1:
                print("cut stimuli")
            for i in range(len(self.stimuli[0])):
                try:
                    if self.stimuli[0][i]<self.stimulyDuration:
                        patchValue=self.histMean(data[:self.stimuli[0][i]+self.stimulyDuration])
                    else:
                        patchValue=self.histMean(data[self.stimuli[0][i]-(self.stimuli[1][i]-self.stimuli[0][i]):self.stimuli[0][i]])
                    processedData[self.stimuli[0][i]:self.stimuli[1][i]]=patchValue
                except:
                    self.hardError=1   
            return processedData
        else:
            return data


    def mainLevelFinding(self):
        self.mainLevel=int((math.log((self.baseFrequency*(1/2+sqrt(sqrt(self.snr))/2))/self.frequency,0.5))-2)#
        if self.debug==1:
            print("self.mainLevel,self.snr",self.mainLevel,self.snr)
    

    def filtering(self):        
        self.coeffs=pywt.swt(self.cleanData, self.wavelet, level=self.mainLevel+1)
        self.coefsBeforeF=asmatrix([self.coeffs[i][1] for i in range(len(self.coeffs))])
        for i in range(len(self.coeffs)):
            cA, cD = self.coeffs[i]
            cDbk = cD
            if i>=(len(self.coeffs)-self.highNoiseLevel):
                cD=zeros(len(cA),dtype='float32')
                if self.debug==1:
                    print(("noisLevel",i))
            else:
                minSD=self.stdFinder(cD[self.deltaLen:],self.defaultFrame)
                cD=pywt.thresholding.soft(cD,minSD*(self.coeffTreshold+i**4))
            self.coeffs[i]=cA, cD
        self.coefsAfterF=asmatrix([self.coeffs[i][1] for i in range(len(self.coeffs))])
        self.result=iswt(self.coeffs,self.wavelet)
   

    def spikeFinding(self):
        resultData=self.result
        start=self.defaultFrame/4
        stop=-self.defaultFrame/4
        minimum,minimumValue = extrema(resultData[start:stop],_max = False, _min = True, strict = False, withend = True)
        maximum,maximumValue = extrema(resultData[start:stop],_max = True, _min = False, strict = False, withend = True)
        std=self.stdFinder(self.cleanData[self.deltaLen:],self.defaultFrame)
        self.signalStd=std
        SD=float16(std+std*sqrt(sqrt(self.snr))/2)#-5.0*self.coeffTreshold/self.snr))#? maybe we must add the snr check?
        if self.debug==1:
            print ((self.snr,std,SD,"self.snr,std,SD"))
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
                tmpObject.spikeLength=self.getSpikeLength(spikePoints[i][0],spikePoints[i][1],spikePoints[i][2])
                tmpObject.calculate()
        if self.debug==1:
            print((self.fileName,self.spikeDict))
            
    def getSpikeLength(self,max1,min1,max2):
        h1=self.result[max1]-self.result[min1]
        h1Part=self.result[max1]-h1*(1-0.5)#50% of first spike front
        firstPoint=where(self.result[max1:min1]<h1Part)[0][0]
        h2=self.result[max2]-self.result[min1]
        h2Part=self.result[max2]-h2*(1-0.5)#50% of first spike front
        secondPoint=where(self.result[min1:max2]<h2Part)[0][0]
        length=(secondPoint-firstPoint)*2
        return length
     
            
    def responsLength(self):
        responsMatrix=zeros((max(self.clusters),2),dtype=int)#[[start1,stop1],[start2,stop2]]
        length=len(self.result)
        smallFrame=self.defaultFrame/10
        if self.debug==1:
            print(len(unique(self.clusters)),len(self.clusters),"len(unique(self.clusters)),len(self.clusters)")
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
                    while((abs(self.result[k:k+smallFrame*4].mean()-baseLevel)>std2/4 or self.result[k:k+smallFrame*4].std()>std2/4) and (k<length-smallFrame*4 and k<self.stimuli[0][i])):
                        k+=smallFrame
                except:
                    if self.debug==1:
                        print "finding end of last response:", sys.exc_info()
                    k=lastMax
                    if k>length-smallFrame*4:
                        k=length-1
                    else:
                        while(abs(self.result[k:k+smallFrame*4].mean()-baseLevel)>std2/6 or self.result[k:k+smallFrame*4].std()>std2/6):
                            k+=smallFrame
                            if k>length-smallFrame*5:
                                k=length-1
                                break
                        
                stop=k
                if self.debug==1:
                    print((lastMax,stop,length,"lastMax,stop,length"))
                try:
                    sampleLen=stop-lastMax
                    if self.debug==1:
                        print(("respons end sample length:",sampleLen))
                    if sampleLen>0:
                        ar=polyfit(array(range(sampleLen)),self.result[lastMax:stop],2)
                        sample=polyval(ar,array(range(sampleLen)))
                        extrem=where(diff(sample)==0)[0]
                        if self.debug==1:
                            print(("number of extremums in respons length curve:", len(extrem)))
                        if len(extrem)>0:
                            lastExtremum=extrem[-1]
                            if sample[lastExtremum]>sample[lastExtremum+1]:
                                realStop=stop
                            else:
                                realStop=lastExtremum
                        else:
                            if sample[0]>sample[1]:
                                shift=where(diff(sample)>-0.000001)[0]
                                if len(shift)>0:
                                    realStop=lastMax+shift[0]
                                else:
                                    realStop=stop
                            else:
                                realStop=stop
                        responsMatrix[i-1]=start,realStop
                    else:
                        responsMatrix[i-1]=start,stop
                except:
                    responsMatrix[i-1]=start,stop
                    print "Unexpected error in response length finding:", sys.exc_info()
        return responsMatrix
    
    def epspReconstructor(self,tmpObject):
        tmpObject2=getattr(self,tmpObject.spikes[-1])
        sample1=self.result[tmpObject.responseStart:tmpObject2.spikeMax2]     
        mask=zeros(tmpObject.responseEnd-tmpObject.responseStart,dtype='bool')
        tmpObject2=getattr(self,tmpObject.spikes[0])
        mask[0]=1
        mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
        mask[tmpObject2.spikeMax1-tmpObject.responseStart]=1
        mask[(tmpObject2.spikeMax1-tmpObject.responseStart)/2]=1
        if len(tmpObject.spikes)>=2:
            try:
                tmpObject2=getattr(self,tmpObject.spikes[1])
                mask[tmpObject2.spikeMax1-tmpObject.responseStart]=1
                mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
            except:
                print "there is no second spike?:", sys.exc_info()
        tmpObject2=getattr(self,tmpObject.spikes[-1])
        sample2=self.result[tmpObject2.spikeMax2:tmpObject.responseEnd]
        sample2Points=array(range(len(sample2)))
        if len(sample2Points)>301:
            ar1=Rbf(sample2Points[sample2Points%100==0],sample2[sample2Points%100==0],smooth=0.01)#,function='gaussian')
        else:
            ar1=Rbf(sample2Points[sample2Points%10==0],sample2[sample2Points%10==0],smooth=0.01)
        xr1=ar1(sample2Points)
        step=len(sample2Points)/8
        for i in range(8):
            mask[tmpObject2.spikeMax2-tmpObject.responseStart+step*i]=1
        mask[-1]=1
        sample=append(sample1,xr1)
        if self.debug==1:
            print((len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd]),"len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd])"))
        timePoints=where(mask==True)[0]
        values=sample[mask]
        sample3=array(range(len(sample)))
        try:           
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
                return 0,0      
        except:
            print "Unexpected error in reconstruction:", sys.exc_info()
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
            tmpObject.epspFront=0
            tmpObject.epspBack=0
            try:
                tmpObject.epspFront,tmpObject.epspBack = self.epspReconstructor(tmpObject)
            except:
                print "Unexpected error wile response %s reconstruction:" % i, sys.exc_info()
                self.hardError=1
        print(self.fileName.split('/')[-1],self.responseDict)
        

    def plotData(self):
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        if self.debug==1:
            fig, axes_list = plt.subplots(3, 1, sharex=True)
            ax = axes_list[0]
        else:
            fig, ax = plt.subplots(1, 1)
        ax.plot(self.data,'y')
        try:
            ax.plot(self.result,'b')
            try:
                ax.plot(self.epsp[0],self.epsp[1],'r')
            except:
                self.softError=1
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
            except:
                print "Unexpected error wile ploating:", sys.exc_info()
                self.hardError=1
            for i in range(len(self.stimuli[0])):
                ax.axvline(x=self.stimuli[0][i],color='g')
        except:
            print "Unexpected error wile ploating computedData:", sys.exc_info()
            self.hardError=1
         
        if self.debug==1:   
            #    ax.locator_params(nbins=3)
            bx = axes_list[1]
            bx.set(xlabel="x-label", ylabel="y-label", title="before filtering")
            #print self.coefsBeforeF.shape
            normMatrixBefore=self.coefsBeforeF[:3]/sqrt(self.coefsBeforeF[:3].var())
            bx.imshow(normMatrixBefore, aspect='auto')
            #print(self.coefsBeforeF.size())
            cx = axes_list[2]
            cx.set(xlabel="x-label", ylabel="y-label", title="after filtering")
            normMatrixAfter=self.coefsAfterF[:3]/sqrt(self.coefsAfterF[:3].var())
            cx.imshow(normMatrixAfter, aspect='auto')
            plt.tight_layout()
            #
            fig2, ax2 = plt.subplots(1, 1)
            ax2.plot(self.data,'y')
            #
            plt.show()
        else:
            plt.savefig(self.fileName+"_graph.png")
            #plt.show()
            plt.close()# very important to stop memory leak
        
        del fig
        
    def writeData(self):
        if self.write:
            self.mysql_writer.dbWriteSignalProperties(self.signalPtp,self.snr,self.signalStd,self.mainLevel)
            for i in self.responseDict.values():
                tmpObject=getattr(self,i)
                self.mysql_writer.dbWriteResponse(tmpObject)
                for j in tmpObject.spikes:
                    tmpObject2=getattr(self,j)
                    self.mysql_writer.dbWriteSpike(tmpObject2)
                #    del tmpObject2
                #del tmpObject
            
            
            
            
            
            
            
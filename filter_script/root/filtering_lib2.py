#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''
#library for filtering
import sys
from mimetypes import guess_type
from numpy import zeros, log, asmatrix, append, arange, math, empty, sqrt, histogram, ones, diff, loadtxt, fromfile, array, int16, unique, where,float16,float32,argmin,median
import pywt
import os    
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('agg') 
from matplotlib.patches import Rectangle
from matplotlib.widgets import Cursor
import matplotlib.pyplot as plt
from root.externalFunctions import iswt,extrema
from root.objects import Spike,Response
from scipy.interpolate import Rbf
from scipy import polyval, polyfit, signal
from root.clussterization import clusterization, clusterAnalyser
from root.checkResult import resultAnalysis
import root.rInterface as rInterface 

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
            self.frequency=self.freqLoad(self.fileName,self.defFrequency)
        except:
            print "self.freqLoad() error:", sys.exc_info()
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
            if self.debug==1:        
                fig, ax = plt.subplots(1, 1)
                ax.grid(color='k', linestyle='-', linewidth=0.4)
                ax.plot(self.data)
                #ax.plot(self.cleanData,'g')       
                plt.show()   
                plt.close()# very important to stop memory leak
                del fig,ax
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
            self.resultRough=self.filtering(self.coeffTreshold)
            self.result=self.filtering(self.coeffTreshold-5)
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
            resultAnalysis(self,self.debug)
        except:
            print "resultAnalysis error:", sys.exc_info()
            self.hardError=1        


    def argReading(self):
        self.wavelet='sym3'#'bior3.5
        self.defFrequency = int(self.arguments[2])
        self.destination = str(self.arguments[3])
        self.coeffTreshold = int(self.arguments[4])+7#4
        self.debug = int(self.arguments[6])
        self.write = int(self.arguments[7])
        self.isClusterOn = int(self.arguments[8])
        try:
            self.manual = int(self.arguments[9])
        except:
            self.manual = 1

    def dataLoading(self):
        #print("mimetype = ",guess_type(self.fileName)[0])
        #if guess_type(self.fileName)[0] == 'text/plain' or guess_type(self.fileName)[0] == 'chemical/x-mopac-input':
        tmpData = loadtxt(self.fileName,dtype='float32')
        data=self.amplLoad(self.fileName,tmpData)
        self.data = self.dataFitting(data,self.frequency)
        #else:
        #    tmpData=fromfile(self.fileName,int16)
        #    data=self.amplLoad(self.fileName,tmpData)
        #    self.data=self.dataFitting(data,self.frequency)
    
    def freqLoad(self,filename,defFrequ):
        try:
            dictionary=dict()
            iniName=filename.strip('.dat')+".ins"
            print(filename)
            #print(iniName)
            fd=file(iniName,'r')
            for line in fd.readlines():
                pairs = line.split("=")
                for variable, value in zip(pairs[::2],pairs[1::2]):
                    dictionary[variable] = value.strip()
            frequency=int(dictionary['CodingFreq'].strip())*1000
        except:
            print "freqLoad() error:", sys.exc_info()
            frequency=defFrequ
        return frequency
    
    def amplLoad(self,filename,data):
            if ("5мв" in filename) or ("5mv" in filename):
                print("# 5mv amplifier #")
                return(data*5.0/2)
            else:
                return data        
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
        dataLen=len(data)
        tmp=2**(math.ceil(math.log(dataLen,2)))
        self.deltaLen=tmp-dataLen
        meanTmp=self.histMean(data[:dataLen/4])
        newData=empty(tmp, dtype='float32')
        newData.fill(meanTmp)
        newData[self.deltaLen:]=data
        return newData

    def tresholdCreating(self):
        msec=self.frequency/1000 # 1 msec =  frequency/1000 points        
        self.defaultFrame=4*msec #frame size for mean() and std() finding must depend on frequency. assume it equal to 4 msec
        self.stimulyDuration=int(0.8*msec) #1msec# treshold for stimuli filtering ~20points==2msec==2*self.msec
        self.level=int((math.log(100.0/self.frequency,0.5)-1)) #wavelet decomposition level. level 6 to 10kHz signal.
        self.baseFrequency=300 #we must separate the levels of wavelet decomposition that contains most part of the signal
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
        if self.debug==1:
            print((data,len(data),frameSize,"data,len(data),framesize"))
        minSD=self.stdFinder(data,frameSize)
        maxSD=self.getLocalPtp(data,frameSize*0.8)
        snr=float16(maxSD/minSD)
        self.signalPtp=maxSD
        self.signalStd=minSD
        if self.debug==1:
            print((minSD,maxSD,snr,"minSD,maxSD,snr in snrFinding function"))
        return snr 
    
    def getLocalPtp(self,data,framesize):
        ptpList=[]
        smoothedData=signal.medfilt(data, 5)
        smoothedDataLen=len(smoothedData)
        if self.debug==1:
            print((data,len(data),framesize,"data,len(data),framesize"))
        try:
            ptpList=[i.ptp() for i in [smoothedData[j:j+framesize] for j in arange(0,smoothedDataLen-framesize,framesize/3)]]
        except:
            print "Unexpected error in finding local ptp:", sys.exc_info()
        if len(ptpList)!=0:
            return max(ptpList)
        else:
            return 0

    def findStimuli(self,data):
        wavelet='haar'
        filterSize=pywt.Wavelet(wavelet).dec_len
        pwr=pywt.swt(data, wavelet, 2)
        pwr2=array(pwr[0][1])
        pwr4Std=self.stdFinder(pwr2[self.deltaLen:], self.defaultFrame)
        pwr4ptp=pwr2.ptp()
        treshold=pwr4Std*1.5*log(pwr4ptp/pwr4Std)#11 - empirical finding coef         
        #treshold=pwr4Std*1.5*log(pwr4ptp/pwr4Std)/1.1
        self.stimTreshold=treshold
        pwr5=zeros(len(pwr2))
        pwr5[1:]+=abs(diff(pwr2))
        pwr5[0]=pwr5[1]
        pwr5[pwr5<treshold]=0
        pwr5[pwr5>0]=treshold
        #for i in pwr5:
        #    try:
        #        if(pwr5[i]==treshold and pwr5[i+2]==treshold):
        #            pwr5[i+1]=treshold
        #    except:
        #        pass
        self.HiFrequNoise2=pwr5
        self.HiFrequNoise1=pwr2
        dpwr=where((diff(pwr5)==treshold)==True)[0]
        dpwrMask=ones(len(dpwr),dtype='bool')
        for i in range(len(dpwr)-1):
            if dpwr[i+1]-dpwr[i]<self.stimulyDuration/3.0:
                dpwrMask[i+1]=0
        dpwr=dpwr[dpwrMask]
        dpwrMask=ones(len(dpwr),dtype='bool')
        for i in range(len(dpwr)):
            length1=len(where((abs(diff(pwr2[dpwr[i]+filterSize:dpwr[i]+filterSize+self.stimulyDuration/2]))>=treshold)==True)[0])
            sample1=data[dpwr[i]+filterSize*2:dpwr[i]+filterSize*2+self.stimulyDuration/2]
            sample2=data[dpwr[i]+filterSize*2-self.stimulyDuration/2:dpwr[i]+filterSize*2]
            sample3=pwr2[dpwr[i]+filterSize:dpwr[i]+filterSize+self.stimulyDuration/2]
            sample4=pwr2[dpwr[i]+filterSize-self.stimulyDuration/2:dpwr[i]+filterSize]
            sample5=data[dpwr[i]+filterSize*2-self.stimulyDuration/2:dpwr[i]+filterSize*2+self.stimulyDuration/2]
            #sampleSumDiff=abs(sample1).sum()/abs(sample2).sum()
            sampleSumDiff=abs(median(sample1)-median(sample2))
            sampleStdDiff=sample3.std()/sample4.std()
            samplePtpDiff=abs(sample3).mean()/abs(sample4).mean()
            sampleGlobalStd=median(sample5)/median(sample2)
            if sampleSumDiff>0:
                #neuroTestResult=rInterface.stimNeuroCheck(length1*100000.0/self.frequency,sample1.max()*0.1,sample1.min()*0.1,sample1.std(),median(sample1)*0.1,sample1.mean()*0.1,sample2.mean()*0.1,median(sample2)*0.1,sample2.std(),median(diff(sample1*1.0)),median(diff(sample2*1.0)),sampleStdDiff,samplePtpDiff)>=0.5
                neuroTestResult=rInterface.stimNeuroCheck(length1*100000.0/self.frequency,sample1.ptp()*0.1,sample2.ptp()*0.1,sample1.std(),median(sample1)*0.1,sample1.mean()*0.1,sample2.mean()*0.1,median(sample2)*0.1,sample2.std(),median(diff(sample1*1.0)),median(diff(sample2*1.0)),sampleStdDiff,samplePtpDiff)>=0.5
                #print(neuroTestResult)
                if neuroTestResult:
                    print("Accepted","start,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd",dpwr[i]+filterSize,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd)
                    #self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"1",self.frequency,sampleStdDiff,samplePtpDiff)
                    pass
                else:
                    dpwrMask[i]=0
                    #self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"0",self.frequency,sampleStdDiff,samplePtpDiff)
                    print("Dropped","start,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd",dpwr[i]+filterSize,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd)
            else:
                dpwrMask[i]=0
                #self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"0",self.frequency,sampleStdDiff,samplePtpDiff)
                print("Dropped","start,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd",dpwr[i]+filterSize,samplePtpDiff,sampleStdDiff,sampleSumDiff,sampleGlobalStd)
            #
        
            #===================================================================
            # if (sampleStdDiff>1.1 and samplePtpDiff>1.2 and (sampleGlobalStd>1.005 or sampleGlobalStd<0.995)):
            #    print(("len of stimulum",i,dpwr[i]+filterSize,length1,length1/(self.stimulyDuration/(self.frequency/16000.0)),self.stimulyDuration/(self.frequency/16000.0),sampleSumDiff,sampleStdDiff,samplePtpDiff))
            #    if length1>=self.stimulyDuration/12 or (sampleSumDiff>10 and sampleStdDiff>1.2 and samplePtpDiff>1.5) or (sampleSumDiff>10 and samplePtpDiff>1.5 and sampleStdDiff>1.2) or (sampleSumDiff>50) or (sampleStdDiff>1.4 and samplePtpDiff>1.6)or (sampleStdDiff>1.6 and samplePtpDiff>1.4) or (sampleStdDiff>1.4 and samplePtpDiff>1.5 and sampleGlobalStd>0.98) or (sampleStdDiff>1.5 and samplePtpDiff>1.4 and sampleGlobalStd>0.98):
            #        if self.write:
            #            try:
            #                self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"1",self.frequency,sampleStdDiff,samplePtpDiff) 
            #            except:
            #                print "Unexpected error wile dbWriteStim:", sys.exc_info()
            #        else:
            #            pass
            #    #elif sampleSumDiff<0.5 and length1>=self.stimulyDuration/20:
            #    #    pass
            #    else:
            #        if self.write:
            #            try:
            #                self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"0",self.frequency,sampleStdDiff,samplePtpDiff) 
            #            except:
            #                print "Unexpected error wile dbWriteStim:", sys.exc_info()
            #        else:
            #            pass
            #        dpwrMask[i]=0
            #        print("stimulum dropped")
            # else:
            #    if self.write:
            #        try:
            #            self.mysql_writer.dbWriteStim(sample2,sample1,length1,i,"0",self.frequency,sampleStdDiff,samplePtpDiff) 
            #        except:
            #            print "Unexpected error wile dbWriteStim:", sys.exc_info()
            #    else:
            #        pass
            #    dpwrMask[i]=0
            #===================================================================
        dpwr=dpwr[dpwrMask]
        if self.debug==1:
                print("number of finded stimuls - %s" % len(dpwr))          
        stimList=[[],[]]
        for i in range(len(dpwr)):
            start=dpwr[i]+filterSize*2
            length=self.stimulyDuration/4
            if self.stimulyDuration/7>5:
                baseline=self.histMean(data[start-int(self.stimulyDuration/7):start])
                baseStd=data[start-int(self.stimulyDuration/7):start].std()
            else:
                baseline=data[start-1]
                baseStd=data[start-6:start-1].std()
            tmpStop=start+length
            if tmpStop+self.defaultFrame+5>len(data):
                tmpStop=len(data)-self.defaultFrame-5
            stimMean=data[start:tmpStop].mean()
            if self.debug==1:
                print((baseline,baseStd,start,tmpStop,stimMean,"baseline,baseStd,start,tmpStop,stimMean"))
            try:
                sample=signal.medfilt(data[tmpStop-5:tmpStop+self.defaultFrame+5],int(self.stimulyDuration/20))
            except:
                sample=data[tmpStop-int(self.stimulyDuration/40):tmpStop+self.defaultFrame+int(self.stimulyDuration/40)]
            firstArray=abs(diff(sample))<=baseStd
            secondArray=abs(sample[1:]-baseline)<baseStd/2
            fourthArray=self.stdArray(sample[1:],5)<=baseStd*3
            fivthArray=abs(sample[1:]-baseline)<baseStd*4
            if self.debug==1:
                print((any(firstArray==True),any(secondArray==True),any(fourthArray==True),any(fivthArray==True)))
            try:
                shift=where((firstArray*fivthArray+secondArray)*fourthArray==True)[0]
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

    def stdArray(self,sample,frame):
        test=empty(len(sample)+frame*2)
        result=empty(len(sample))
        test[frame:-frame]=sample
        test[:frame]=sample[:frame]
        test[-frame:]=sample[-frame:]
        for i in range(len(sample)):
            result[i]=test[i:i+frame*2].std()
        return result

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
        if self.debug==1:
            print("self.snr",self.snr)
        self.mainLevel=int(math.log((self.baseFrequency*(1/2+sqrt(sqrt(self.snr))/2))/self.frequency,0.5)-1.4)#
        if self.debug==1:
            print("self.mainLevel:",self.mainLevel)    

    def filtering(self,coeffTreshold):        
        self.coeffs=pywt.swt(self.cleanData, self.wavelet, level=self.mainLevel+1)
        self.coefsBeforeF=asmatrix([self.coeffs[i][1] for i in range(len(self.coeffs))])
        for i in range(len(self.coeffs)):
            cA, cD = self.coeffs[i]
            if i>=(len(self.coeffs)-self.highNoiseLevel):
                cD=zeros(len(cA),dtype='float32')
                if self.debug==1:
                    print(("noisLevel",i))
            else:
                #if self.debug==1:        
                #    fig, ax = plt.subplots(1, 1)
                #    ax.grid(color='k', linestyle='-', linewidth=0.4)
                #    ax.plot(cD)       
                #    plt.show()   
                #    plt.close()# very important to stop memory leak
                #    del fig,ax
                minSD=self.stdFinder(cD[self.deltaLen:],self.defaultFrame)
                maxSD=self.getLocalPtp(cD[self.deltaLen:],self.defaultFrame*0.8)
                snr=maxSD/minSD
                if self.debug==1:
                    print(minSD,maxSD,snr,i,minSD*(coeffTreshold*(snr**(i**(0.7)/(i**(1.5)+i+1)))+i*2))
                cD=pywt.thresholding.soft(cD,minSD*(coeffTreshold*(snr**(i**(0.7)/(i**(1.5)+i+1)))+i*2))
            self.coeffs[i]=cA, cD
        self.coefsAfterF=asmatrix([self.coeffs[i][1] for i in range(len(self.coeffs))])
        return iswt(self.coeffs,self.wavelet)   

    def spikeFinding(self):
        resultDataForSearch=self.resultRough
        resultData=self.result
        start=self.defaultFrame/4
        stop=-self.defaultFrame/4
        minimum,minimumValue = extrema(resultDataForSearch[start:stop],_max = False, _min = True, strict = False, withend = True)
        maximum,maximumValue = extrema(resultDataForSearch[start:stop],_max = True, _min = False, strict = False, withend = True)
        minimumValue=resultData[minimum+start]
        tmpMaximum=maximum.tolist()
        tmpMaximum.extend(array(self.stimuli[1])-start)
        maximum=array(tmpMaximum)
        maximum.sort()
        maximumValue=resultData[maximum+start]
        #if self.debug==1:        
            #fig, ax = plt.subplots(1, 1)
            #ax.grid(color='k', linestyle='-', linewidth=0.4)
            #ax.plot(resultDataForSearch[start:stop],'b')
            #for i in range(len(minimum)):
            #    ax.axvline(x=minimum[i],color='g')
            #for i in range(len(maximum)):
            #    ax.axvline(x=maximum[i],color='r')
            #plt.show()   
            #plt.close()# very important to stop memory leak
            #del fig,ax
        #=======================================================================
        # if self.debug==1:        
        #    fig, ax = plt.subplots(1, 1)
        #    ax.grid(color='k', linestyle='-', linewidth=0.4)
        #    ax.plot(resultDataForSearch,'b')
        #    ax.plot(resultData,'g')
        #    plt.show()   
        #    plt.close()# very important to stop memory leak
        #    del fig,ax
        #=======================================================================
        std=self.stdFinder(self.cleanData[self.deltaLen:],self.defaultFrame)
        SD=float16(std+std*self.snr**(1/4)/4)#-5.0*self.coeffTreshold/self.snr))#? maybe we must add the snr check?
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
            if maximumValue[tmpMaximum1]-minimumValue[i]>SD and maximumValue[tmpMaximum2]-minimumValue[i]>SD\
             and (maximum[tmpMaximum2]-maximum[tmpMaximum1])>self.stimulyDuration/2:
                spikePoints.append([start+maximum[tmpMaximum1],start+minimum[i],start+maximum[tmpMaximum2]])
        for i in range(len(spikePoints)):
            try:
                ampl=round(resultData[spikePoints[i][0]]-resultData[spikePoints[i][1]]+\
                       (resultData[spikePoints[i][2]]-resultData[spikePoints[i][0]])/\
                       (spikePoints[i][2]-spikePoints[i][0])*(spikePoints[i][1]-spikePoints[i][0]),1)
                if self.debug==1:
                    print(("Len of SpikePoints=%s,Spike %s, ampl=%s,SD=%s" % (len(spikePoints),i,ampl,SD)))
                if ampl>SD:
                    index=len(self.spikeDict)
                    self.spikeDict[index]="n"+str(i)
                    setattr(self,self.spikeDict[index],Spike(self.frequency))
                    tmpObject=getattr(self,self.spikeDict[index])                
                    tmpObject.responseStart=start
                    tmpObject.responseEnd=stop
                    tmpObject.spikeMax1=spikePoints[i][0]
                    tmpObject.spikeMax1Val=self.result[spikePoints[i][0]]
                    tmpObject.spikeMin=spikePoints[i][1]
                    tmpObject.spikeMinVal=self.result[spikePoints[i][1]]
                    tmpObject.spikeMax2=spikePoints[i][2]
                    tmpObject.spikeMax2Val=self.result[spikePoints[i][2]]
                    tmpObject.spikeMaxToMin=tmpObject.spikeMax2Val-tmpObject.spikeMinVal
                    tmpObject.spikeAmpl=ampl
                    tmpObject.manual=self.manual
                    tmpObject.spikeLength=self.getSpikeLength(spikePoints[i][0],spikePoints[i][1],spikePoints[i][2])*1000.0/self.frequency
                    tmpObject.spikeFront,tmpObject.spikeBack=self.getSpikeAngles(resultData[spikePoints[i][0]:spikePoints[i][2]])
            except:
                print "Unexpected error in finding of stimuli end:", sys.exc_info()
        if self.debug==1:
            print((self.fileName,self.spikeDict))
        
    def checkForFibrePotential(self,spikeList,respDictValue):
        try:
            #print(self.responseDict.values())
            tmpObject=getattr(self,respDictValue)
            tmpObject2=getattr(self,tmpObject.spikes[0])
            if rInterface.neuroCheck(tmpObject2.spikeMax2Val-tmpObject2.spikeMinVal,tmpObject2.spikeDelay,tmpObject2.spikeLength,tmpObject2.spikeFront,tmpObject2.spikeBack,self.frequency/1000)>=0.5:
                tmpObject2.fibre=1
                print("There are AP at zero position")
        except:
            print "Unexpected error wile fibre potential checking:", sys.exc_info()
            
                     
    def interactiveFibreSearch(self,spikeList,respDictValue):
        try:
            print(self.responseDict.values())
            tmpObject=getattr(self,respDictValue)
            tmpObject2=getattr(self,tmpObject.spikes[0])
            start=tmpObject2.spikeMax1-(tmpObject2.spikeMin-tmpObject2.spikeMax1)
            tmpObject2=getattr(self,tmpObject.spikes[-1])
            stop=tmpObject2.spikeMax2+(tmpObject2.spikeMax2-tmpObject2.spikeMin)
            fig, ax = plt.subplots(1, 1)
            #self.overview=plt.subplot2grid((2,2),(0,0),rowspan=2,colspan=2)
            ax.grid(color='k', linestyle='-', linewidth=0.4)
            ax.plot(self.cleanData,'y')
            ax.plot(self.result,'b')
            ax.set_xlim((start,stop))
            tmpObject=getattr(self,respDictValue)
            try:
                self.minimumsForFibre=[]
                print(tmpObject.spikes)
                for j in tmpObject.spikes:
                    print(j)
                    tmpObject2=getattr(self,j)
                    ax.plot(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin],'or')
                    ax.text(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin]-15, str(j), fontsize=12, va='bottom')
                    self.minimumsForFibre.append(tmpObject2.spikeMin)
                          
            except:
                print "Unexpected error wile spike ploating:", sys.exc_info()
            self.fibreIndex=''
            cursor=Cursor(ax, useblit=True, color='black', linewidth=2 )
            _widgets=[cursor]
            fig.canvas.mpl_connect('button_press_event',self.click)
            plt.show()   
            # very important to stop memory leak
            del fig,ax
            if self.fibreIndex!='':
                #print(tmpObject.spikes[self.fibreIndex])
                print("Fibre spike: %s" % tmpObject.spikes[self.fibreIndex])
            else:
                print "O_o"
            fibre=tmpObject.spikes[self.fibreIndex]#input()
            if fibre=='':
                pass
            else:
                tmpObject2=getattr(self,fibre)
                tmpObject2.fibre=1
        except:
            print "Unexpected error wile ploating:", sys.exc_info()


    def click(self,event):
        """
        What to do, if a click on the figure happens:
            1. Check which axis
            2. Get data coord's.
            3. Plot resulting data.
            4. Update Figure
        """
        #if event.inaxes==ax:
        #Get nearest data
        xpos=argmin(abs(event.xdata-self.minimumsForFibre))
        #Check which mouse button:
        if event.button==1:
            #Plot it                
            self.fibreIndex=xpos

        elif event.button==3:
            self.fibreIndex=''
        #plt.draw()
        plt.close()
        
    def getSpikeAngles(self,sample):
        sampleMin = min(sample)
        try:
            minPoint = where(sample==sampleMin)[0][0]
        except:
            minPoint = where(sample==sampleMin)[0]
        max1 = sample[0]
        max2 = sample[-1]
        h1 = max1-sampleMin
        point1 = where(sample[:minPoint]<(max1-h1*0.2))[0][0]
        point2 = where(sample[:minPoint]<(max1-h1*0.8))[0][0]
        h2 = max2-sampleMin
        point3 = where(sample[minPoint:]>(sampleMin+h2*0.2))[0][0]
        point4 = where(sample[minPoint:]>(sampleMin+h2*0.8))[0][0]
        
        if point2-point1>0:
            angle1 = (sample[point1]-sample[point2])/(point2-point1)
        else:
            angle1 = sample[point1]-sample[point2]
        if point4-point3>0:
            angle2 = (sample[point4]-sample[point3])/(point4-point3)
        else:
            angle2 = sample[point4]-sample[point3]
        return angle1*1000.0/self.frequency,-angle2*1000.0/self.frequency
            
    def getSpikeLength(self,max1,min1,max2):
        h1=self.result[max1]-self.result[min1]
        h1Part=self.result[max1]-h1*(1-0.5)#50% of first spike front
        h2=self.result[max2]-self.result[min1]
        h2Part=self.result[max2]-h2*(1-0.5)#50% of first spike front
        try:
            firstPoint=where(self.result[max1:min1]<h1Part)[0][0]        
            secondPoint=where(self.result[min1:max2]>h2Part)[0][0]
        except:
            try:
                firstPoint=where(self.result[max1:min1]<h1Part)[0]        
                secondPoint=where(self.result[min1:max2]>h2Part)[0]
            except:
                firstPoint=(min1-max1)/2
                secondPoint=(max2-min1)/2
        length=(secondPoint+(min1-(max1+firstPoint)))*2
        if self.debug==1:
            print "spike lendth finding"
            print((max1,min1,max2,h1,h1Part,h2,h2Part,firstPoint,secondPoint,length,"\nmax1,min1,max2,h1,h1Part,h2,h2Part,firstPoint,secondPoint,length"))
        return length
     
            
    def responsLength(self):
        """
        The function for calculation the response length - the distance between the stimuli and the moment of potential returning to baseline 
        """
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
                """
                The case when distinct responses are separated by clusterization method rather then stimuli
                This mean that response start search starts from the first max of the first spike
                """
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
                """
                The case when distinct responses are separated by electrical stimuli search
                This mean that the response starts from electrical artifact
                """
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
                """
                Shift the response end to left for long tail preventing
                """
                difference=self.result[lastMax]-self.result[stop]
                threshold=self.result[lastMax]-0.9*difference
                if difference>0:                    
                    while self.result[stop]<threshold:
                        stop-=1
                if difference<0:
                    while self.result[stop]>threshold:
                        stop-=1
                if self.debug==1:
                    print((lastMax,stop,length,"lastMax,stop,length"))
                try:
                    """
                    Calculate response end precisely using polynomial smoothing 
                    """
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
        """
        The function for reconstructing shape of the inverted fEPSP
        Where:
            tmpObject - fEPSP event object
        """
        tmpObject2=getattr(self,tmpObject.spikes[-1])
        sample1=self.result[tmpObject.responseStart:tmpObject2.spikeMax2]     
        mask=zeros(tmpObject.responseEnd-tmpObject.responseStart,dtype='bool')
        tmpObject2=getattr(self,tmpObject.spikes[0])
        print("Is the first spike a fibre spike? : %s" % str(tmpObject2.fibre))
        if(tmpObject2.fibre==1):
            #mask[0]=1
            mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
            mask[(tmpObject2.spikeMin+tmpObject2.spikeMax2)/2-tmpObject.responseStart]=1
            epspStart=tmpObject2.spikeMin-tmpObject.responseStart
            mask[epspStart]=1
        else:
            intend=self.getFirstMaxIndent(tmpObject2.spikeMax1,self.result)
            #print(intend)
            mask[tmpObject2.spikeMax1-tmpObject.responseStart]=1
            mask[tmpObject2.spikeMax2-tmpObject.responseStart]=1
            epspStart=tmpObject2.spikeMax1-intend-tmpObject.responseStart
            mask[epspStart]=1
            mask[epspStart+intend/2]=1
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
        if len(sample2Points)>30:
            if len(sample2Points)>301:
                ar1=Rbf(sample2Points[sample2Points%100==0],sample2[sample2Points%100==0],smooth=0.01)#,function='gaussian')
            else:
                ar1=Rbf(sample2Points[sample2Points%10==0],sample2[sample2Points%10==0],smooth=0.01)
            xr1=ar1(sample2Points)
            self.epileptStd=self.calculateEpilept(sample2,xr1)
            step=len(sample2Points)/8
            for i in range(8):
                mask[tmpObject2.spikeMax2-tmpObject.responseStart+step*i]=1
            
            sample=append(sample1,xr1)
            if self.debug==1:
                print((len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd]),"len(mask),len(sample),len(self.result[tmpObject.responseStart:tmpObject.responseEnd])"))
        else:
            sample=sample1
            self.epileptStd=0
        mask[-1]=1    
        timePoints=where(mask==True)[0]
        values=sample[mask]
        sample3=array(range(len(sample)))
        try:           
            sp3 = Rbf(timePoints,int16(values),smooth=1,function='thin_plate')
            y4=sp3(sample3)
            try:
                front,back = self.epspAnaliser(y4[epspStart:])
                self.epsp=append(self.epsp,[sample3+tmpObject.responseStart,y4],axis=1)
                if self.debug==1:
                    print((front,back,"front,back"))
                return epspStart,front,back
            except:
                print "Unexpected error in epspAnaliser:", sys.exc_info()
                return epspStart,0,0      
        except:
            print "Unexpected error in reconstruction:", sys.exc_info()
            return epspStart,0,0
        
    def getFirstMaxIndent(self,stop,sample):
        """
        This function used for calculate the distance between first max of the first spike and prior minimal extremum
        in the case when fibre potential didn`t find 
        """
        i=0
        if sample[stop]<sample[stop-1]:
            while(sample[stop-i]<sample[stop-i-1]):
                i+=1
        else:       
            while(sample[stop-i]>=sample[stop-i-1]):
                i+=1
        return i
    
    def calculateEpilept(self,sample2,xr1):
        """
        Function for rouge estimating the "epileptiform activity". The estimation based on the oscillations after second maximum of the last founded spike in response.
        """
        diffSample=sample2-xr1
        return diffSample.std()
        
    def epspAnaliser(self,y):
        """
        Function for different fEPSP properties calculation.
        fEPSP front: "ADC codes"/"ADC points" at the 20%-80% of fEPSP amplitude interval
        fEPSP back: "ADC codes"/"ADC points" at the 80%-20% of fEPSP amplitude interval 
        """
        #preparing
        maxvalue=max(y)
        maxPoint=where(y==maxvalue)[0]
        minvalue1=min(y[:maxPoint])
        minvalue2=min(y[maxPoint:])
        ampl1=maxvalue-minvalue1
        ampl2=maxvalue-minvalue2
        #computation
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
        """
        Main function for response(fEPSP and all spikes induced by one stimulus)
        The following functions are called from there:
        
            self.checkForFibrePotential - is the first "spike" a real spike or fibre potential
            self.getResponsLength - get the length of response in ADC points (20% to 20% of fEPSP amplitude)
            self.setSpikeDelays - calculate the time delay between stimulus and spike`s minimum
            self.epspReconstructor - reconstruct the shape of fEPSP
            self.epileptStd - rouge estimate the "epileptiform activity"
            self.spikeArea - calculation of area between signal curve and fEPSP for one spike
            self.epspArea - calculation of area between signal and fEPSP for all spikes
             
        """
        rMatrix=self.responsMatrix
        self.epsp=array([[],[]])
        epspAreaStart=0
        for i in unique(self.clusters):
            index=len(self.responseDict)
            self.responseDict[index]="r"+str(i)
            setattr(self,self.responseDict[index],Response())
            tmpObject=getattr(self,self.responseDict[index])
            try:
                #tmpObject.fibre,tmpObject.spikes=self.checkForFibrePotential(array(self.spikeDict.values())[self.clusters==i])
                try:
                    tmpObject.spikes=array(self.spikeDict.values())[self.clusters==i]
                    tmpObject.responseStart=rMatrix[i-1][0]
                    tmpObject.responseEnd=rMatrix[i-1][1]
                    tmpObject.length=self.getResponsLength(self.result[rMatrix[i-1][0]:rMatrix[i-1][1]])*1000.0/self.frequency
                    tmpObject.response_top=max(self.result[rMatrix[i-1][0]:rMatrix[i-1][1]])
                    tmpObject.response_bottom=min(self.result[rMatrix[i-1][0]:rMatrix[i-1][1]])
                    tmpObject.baselevel=self.result[rMatrix[i-1][0]-self.defaultFrame/2:rMatrix[i-1][0]].mean()
                    tmpObject.responsNumber=i
                    tmpObject.vpsp=round(tmpObject.response_top-tmpObject.baselevel,1)
                except:
                    print "Unexpected error wile fill response properties:", sys.exc_info()
                try:
                    self.setSpikeDelays(tmpObject.spikes,tmpObject.responseStart)
                except:
                    print "Unexpected error wile setSpikeDelays:", sys.exc_info() 
                if self.manual==1:
                    self.interactiveFibreSearch(array(self.spikeDict.values())[self.clusters==i],self.responseDict[index])
                else:
                    self.checkForFibrePotential(array(self.spikeDict.values())[self.clusters==i],self.responseDict[index])
                self.spikeNumberShift(tmpObject)    
                    #tmpObject.fibre,tmpObject.spikes=self.checkForFibrePotential(array(self.spikeDict.values())[self.clusters==i])
            except:
                print "Unexpected error wile checkForFibrePotential:", sys.exc_info()
            if len(unique(self.clusters))<4:     
                try:
                    epspStartCorrection,tmpObject.epspFront,tmpObject.epspBack = self.epspReconstructor(tmpObject)
                    epspAreaStart,tmpObject.epspArea=self.Area(tmpObject,epspAreaStart,epspStartCorrection)
                    tmpObject.epspEpileptStd=self.epileptStd
                except:
                    print "Unexpected error wile response %s reconstruction:" % i, sys.exc_info()
                    self.hardError=1
        print(self.fileName.split('/')[-1],self.responseDict)
        
    def spikeNumberShift(self,tmpObject):
        tmpObject2=getattr(self,tmpObject.spikes[0])
        #print (tmpObject2.spikeNumber,tmpObject2.fibre)
        if tmpObject2.spikeNumber==0 and tmpObject2.fibre!=1:
            print("spike number shift")
            for i in tmpObject.spikes:
                #print(i)
                tmpObject2=getattr(self,i)
                #print(tmpObject2.spikeNumber)
                tmpObject2.spikeNumber+=1
                #print(tmpObject2.spikeNumber)
             
    def Area(self,tmpObject,epspAreaStart,epspStartCorrection):
        start=tmpObject.responseStart
        stop=tmpObject.responseEnd
        areaArray=self.epsp[1][epspAreaStart+epspStartCorrection:stop-start+epspAreaStart]-self.result[start+epspStartCorrection:stop]
        #fig, ax = plt.subplots(1, 1) 
        #ax.plot(areaArray)
        for j in tmpObject.spikes:
            tmpObject2=getattr(self,j)            
            spikeStart=tmpObject2.spikeMax1-start-epspStartCorrection
            #ax.vlines(spikeStart,100, 0, color='g', linestyles='dashed')
            spikeStop=tmpObject2.spikeMax2-start-epspStartCorrection
            #ax.vlines(spikeStop,100, 0, color='r', linestyles='dashed')
            tmpObject2.area=round(sum(areaArray[spikeStart:spikeStop]))
            #print("spikeArea",j,tmpObject2.area)
            #print(areaArray[spikeStart],areaArray[spikeStop])
        #plt.show()
        area=round(sum(areaArray))
        #print ("area",area,epspAreaStart,epspStartCorrection)
        epspAreaStart=epspAreaStart+stop-start
        return epspAreaStart,area
        
    def getResponsLength(self,sample):
        h1=sample[0]-max(sample)
        h1Part=h1*0.2#20% of first spike front
        h2=sample[-1]-max(sample)
        h2Part=h2*0.2
        try:
            firstPoint=where(sample>h1Part)[0][0]        
            secondPoint=where(sample.revers()>h2Part)[0][0]
        except:
            try:
                firstPoint=where(sample>h1Part)[0]        
                secondPoint=where(sample.revers()>h2Part)[0]
            except:
                firstPoint=0
                secondPoint=0
        length=len(sample)-firstPoint-secondPoint
        return length
        
        
    def setSpikeDelays(self,spikeDict,startPoint):
        for i in spikeDict:
            tmpObject=getattr(self,i)
            tmpObject.spikeDelay=(tmpObject.spikeMin-startPoint)*1000.0/self.frequency
               
    def plotData(self):
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        if self.debug==1:
            fig, axes_list = plt.subplots(3, 1, sharex=True)
            ax = axes_list[0]
        else:
            fig, ax = plt.subplots(1, 1)
        fig.canvas.set_window_title(self.fileName.split("/")[-1]) 
        ax.plot(self.cleanData,'y')
        try:
            ax.plot(self.result,'b')
            try:
                #pass
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
            ax.locator_params(nbins=3)
            bx = axes_list[1]
            bx.set(xlabel="x-label", ylabel="y-label", title="before filtering")
            #print self.coefsBeforeF.shape
            normMatrixBefore=self.coefsBeforeF[:3]/sqrt(self.coefsBeforeF[:3].var())
            bx.imshow(normMatrixBefore, aspect='auto')
            #print(self.coefsBeforeF.size())
            cx = axes_list[2]
            cx.set(xlabel="x-label", ylabel="y-label", title="hiFrequNoise")
            normMatrixAfter=self.coefsAfterF[:3]/sqrt(self.coefsAfterF[:3].var())
            #cx.imshow(normMatrixAfter, aspect='auto')
            cx.plot(abs(diff(self.HiFrequNoise1)),'g')
            
            cx.plot(abs(diff(self.HiFrequNoise2)),'b')
            cx.axhline(self.stimTreshold)
            plt.tight_layout()
            #
            #fig2, ax2 = plt.subplots(1, 1)
            #ax2.plot(self.data,'y')
            #
            plt.show()
        else:
            plt.savefig(self.fileName+"_graph.png")
            #plt.show()
            plt.close()# very important to stop memory leak
        
        del fig
        
        if self.debug==1:        
                fig, ax = plt.subplots(1, 1)
                ax.grid(color='k', linestyle='-', linewidth=0.4)
                ax.plot(self.result,'b')       
                plt.show()   
                plt.close()# very important to stop memory leak
                del fig,ax
        if self.debug==1:        
                fig, ax = plt.subplots(1, 1)
                ax.grid(color='k', linestyle='-', linewidth=0.4)
                ax.plot(self.cleanData,'y')
                ax.plot(self.result,'b')       
                plt.show()   
                plt.close()# very important to stop memory leak
                del fig,ax
        if self.debug==1:        
                fig, ax = plt.subplots(1, 1)
                ax.grid(color='k', linestyle='-', linewidth=0.4)
                ax.plot(self.result,'b')
                try:
                    for i in self.responseDict.values():
                        tmpObject=getattr(self,i)
                        try:
                            for j in tmpObject.spikes:
                                tmpObject2=getattr(self,j)
                                ax.plot(tmpObject2.spikeMin,self.result[tmpObject2.spikeMin],'or')
                                ax.plot(tmpObject2.spikeMax1,self.result[tmpObject2.spikeMax1],'og')
                                ax.plot(tmpObject2.spikeMax2,self.result[tmpObject2.spikeMax2],'og')
                        except:
                            print "Unexpected error wile spike ploating:", sys.exc_info()
                except:
                    print "Unexpected error wile ploating:", sys.exc_info()
                    self.hardError=1       
                plt.show()   
                plt.close()# very important to stop memory leak
                del fig,ax
        
    def writeData(self):
        if self.write:
            try:
                self.mysql_writer.dbWriteSignalProperties(self.signalPtp,self.snr,self.signalStd,self.mainLevel) 
            except:
                print "Unexpected error wile dbWriteSignalProperties:", sys.exc_info()
            try:
                self.mysql_writer.dbWriteNumberOfResponses(len(self.responseDict.values()))
            except:
                print "Unexpected error wile dbWriteNumberOfResponses:", sys.exc_info()
            for i in self.responseDict.values():
                tmpObject=getattr(self,i)
                self.mysql_writer.dbWriteResponse(tmpObject)
        
                try:
                    for j in tmpObject.spikes:
                        tmpObject2=getattr(self,j)
                        self.mysql_writer.dbWriteSpike(tmpObject2)
                except:
                    print "Unexpected error wile dbWriteSpike:", sys.exc_info()
            
            
            
    #===========================================================================
    # if self.debug==1:        
    #    fig, ax = plt.subplots(1, 1)
    #    ax.grid(color='k', linestyle='-', linewidth=0.4)
    #    ax.plot(self.cleanData,'y')       
    #    plt.show()   
    #    plt.close()# very important to stop memory leak
    #    del fig
    #===========================================================================

#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''
from time import strftime,localtime
from os import stat
import sys,MySQLdb,pickle

class Mysql_writer:
    
    def __init__(self,filePath,tagString):
        self.filePath=filePath
        self.variables_global()
        self.dbConnect()
        self.tagString=tagString
        self.rTagDict={"коф":"инкубация","КОФ":"инкубация","ФИЗ":"инкубация","физ":"инкубация","тета":"тетанизация","инк":"инкубация","teta":"тетанизация"}
        self.rTagMask=["до","перед","макс","отмыв"]
        
    def tagWriter(self):
        tagList= self.tagString.split(',')
        if len(tagList)>0:
            cursor = self.conn.cursor()
            for i in tagList:
                tagId = self.tagCheck(i,'tagTable','tagId')
                cursor.execute("INSERT INTO experimentTags(experiment_idexperiment,tagTable_tagId)\
                             VALUES (%s,%s);", (self.idExperiment,tagId))
            self.conn.commit()
            cursor.close()   
        

    def tagCheck(self,tag,table,idName):
        cursor = self.conn.cursor()
        cursor.execute("SELECT {0} FROM {1} WHERE tagName='{2}';".format(idName,table,tag))
        try:
            tagId = cursor.fetchone()[0]
        except:
            print "there are no '%s' tag" % tag, sys.exc_info()
            cursor.execute("INSERT INTO {0}(tagName)\
                             VALUES ('{1}');".format(table,tag))
            self.conn.commit()
            cursor.execute("SELECT {0} FROM {1} WHERE tagName='{2}';".format(idName,table,tag))
            tagId = cursor.fetchone()[0]
        return tagId
                
    def variables_global(self):
        try:
            fd=open("dbConfig",'r')
            dbAccessVars=pickle.load(fd)
            fd.close()
            print(dbAccessVars)
            self.dbServerIp=dbAccessVars[0]
            self.userName=dbAccessVars[2]
            self.userPassword=dbAccessVars[3]
            self.dbName=dbAccessVars[1]
        except:
            print("Load default database config")
            self.dbServerIp='localhost'
            self.userName='filteruser_local'
            self.userPassword='filter123'
            self.dbName='filterdb'
        self.date= strftime('%Y%m%d',localtime(stat(self.filePath).st_mtime))

    def variables_local(self,tmpObject):#
        self.responsNumber=tmpObject.responsNumber
        self.numberofspikes=tmpObject.allSpikes
        self.responslength=tmpObject.responseEnd-tmpObject.responseStart
        self.ampl=tmpObject.spikeAmpl
        self.number=tmpObject.spikeNumber
                
    def dbConnect(self):
        try:
            self.conn = MySQLdb.connect(self.dbServerIp,self.userName ,self.userPassword ,self.dbName );
        except:
            print "Db connect error"
            sys.exit(1)
            
    def dbWriteExperiment(self):
        experimentName=str(self.filePath.split('/')[-2])
        cursor = self.conn.cursor()
        cursor.execute("INSERT INTO experiment(date,experimentName)\
                             VALUES (%s,%s);", (self.date,experimentName))
        self.conn.commit()
        cursor.execute("SELECT idexperiment \
                             FROM experiment \
                             ORDER BY idexperiment\
                             DESC LIMIT 1;")
        self.idExperiment = cursor.fetchone()[0]
        cursor.close()
        
    def dbWriteRecord(self):
        fileName=str(self.filePath.split('/')[-1])
        cursor = self.conn.cursor()
        cursor.execute("INSERT INTO record(filename,time,\
                                            experiment_idexperiment)\
                        VALUES(%s,%s,%s);", (fileName,self.time,str(self.idExperiment)))
        self.conn.commit()
        cursor.execute("SELECT idrecord \
                             FROM record \
                             ORDER BY idrecord\
                             DESC LIMIT 1;")
        self.idRecord = cursor.fetchone()[0]
        cursor.close()
         
    def findTags(self,tagString,tagDict,tagMask):
        tagList=[]
        for i in tagDict.keys():
            if (i in tagString) and (all([j not in tagString for j in tagMask])==True):
                tagList.append(tagDict[i])
        return tagList
        
    def dbWriteRecordTags(self,filename):
        tagList = self.findTags(filename,(self.rTagDict),(self.rTagMask))
        if len(tagList)==0:
            tagList = ["-"]
        #print(tagList)
        cursor = self.conn.cursor()
        for i in tagList:
            tagId = self.tagCheck(i,"recordTags","idrecordTags")
            #print((self.idRecord,tagId))
            cursor.execute("INSERT INTO recordToTags(record_idrecord,recordTags_idrecordTags)\
                             VALUES (%s,%s);", (self.idRecord,tagId))
            self.conn.commit()
        cursor.close()
    
    def dbWriteResponse(self,tmpObject):
        rNumber=tmpObject.responsNumber
        vpsp=tmpObject.vpsp
        nOfSpikes=len(tmpObject.spikes)
        rLength=tmpObject.length
        epspArea=tmpObject.epspArea     
        cursor = self.conn.cursor()
        cursor.execute("INSERT INTO responses(number,numberofspikes,length,record_idrecord,vpsp,epspFront,epspBack,epileptStd,epspArea)\
                             VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)", (str(rNumber),str(nOfSpikes),str(rLength),str(self.idRecord),str(vpsp),\
                             str(tmpObject.epspFront),str(tmpObject.epspBack),str(tmpObject.epspEpileptStd),str(epspArea)))
        self.conn.commit()
        cursor.close()
        
    def dbWriteNumberOfResponses(self,number):
        cursor = self.conn.cursor()
        cursor.execute("UPDATE record\
                        SET numberofresponses=%s\
                        WHERE record.idrecord=%s;", (str(number),str(self.idRecord)))
        self.conn.commit()
        cursor.close()
        
    def dbWriteSignalProperties(self,ptp,snr,std,mainLevel):
        cursor = self.conn.cursor()
        cursor.execute("SELECT idrecord \
                             FROM record \
                             ORDER BY idrecord\
                             DESC LIMIT 1;")
        self.idRecord = cursor.fetchone()[0]
        try:
            cursor.execute("INSERT INTO signalProperties(record_idrecord,ptp,snr,std,mainLevel)\
                             VALUES(%s,%s,%s,%s,%s)", (str(self.idRecord),str(ptp),str(snr),str(std),str(mainLevel)))
        except:
            pass
        self.conn.commit()
        cursor.close()
        
    def dbWriteError(self,soft,hard):
        cursor = self.conn.cursor()
        cursor.execute("UPDATE record\
                        SET softError=%s,hardError=%s\
                        WHERE record.idrecord=%s;", (str(soft),str(hard),str(self.idRecord)))
        self.conn.commit()
        cursor.close()   

    def dbWriteSpike(self,tmpObject):
        ampl=str(tmpObject.spikeAmpl)
        number=tmpObject.spikeNumber
        sLength=str(tmpObject.spikeLength)#must be changed to length at 80% or something like that
        maxdiff=str(tmpObject.spikeMax2Val-tmpObject.spikeMax1Val)
        angle1=str(tmpObject.spikeFront)
        angle2=str(tmpObject.spikeBack)
        delay=str(tmpObject.spikeDelay)
        maxToMin=str(tmpObject.spikeMaxToMin)
        area=str(tmpObject.area)
        fibre=str(tmpObject.fibre)
        cursor = self.conn.cursor()
        cursor.execute("SELECT idresponses\
                             FROM responses\
                             ORDER BY idresponses\
                             DESC LIMIT 1;")
        idResponse = cursor.fetchone()[0]
        cursor.execute("INSERT INTO spikes(ampl,number,responses_idresponses,length,maxDiff,angle1,angle2,delay,maxtomin,area,fibre)\
                             VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);", (ampl,str(number),str(idResponse),sLength,maxdiff,angle1,angle2,delay,maxToMin,area,fibre))
        self.conn.commit()
        cursor.close()
    def dbDisconnect(self):
        self.conn.close()
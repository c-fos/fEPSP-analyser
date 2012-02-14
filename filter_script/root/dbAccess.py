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
    
    def __init__(self,filePath,substance):
        self.filePath=filePath
        self.variables_global()
        self.dbConnect()
        self.substanceName=substance
        self.numberOfResponses=1#must be refactored!
        
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
        cursor = self.conn.cursor()
        cursor.execute("INSERT INTO experiment(date,substance)\
                             VALUES (%s,%s);", (self.date,self.substanceName))
        self.conn.commit()
        cursor.close()
    def dbWriteRecord(self):
        fileName="%s/%s" % (self.filePath.split('/')[-2],self.filePath.split('/')[-1])
        cursor = self.conn.cursor()
        cursor.execute("SELECT idexperiment \
                             FROM experiment \
                             ORDER BY idexperiment\
                             DESC LIMIT 1;")
        idExperiment = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO record(filename,time,numberofresponses,\
                                            experiment_idexperiment)\
                        VALUES(%s,%s,%s,%s);", (fileName,self.time,str(self.numberOfResponses),str(idExperiment)))
        self.conn.commit()
        cursor.close()

    def dbWriteResponse(self,tmpObject):
        #print("dbWriteResponse")
        rNumber=tmpObject.responsNumber
        vpsp=tmpObject.vpsp
        nOfSpikes=len(tmpObject.spikes)
        rLength=tmpObject.responseEnd-tmpObject.responseStart
        
        cursor = self.conn.cursor()
        cursor.execute("SELECT idrecord \
                             FROM record \
                             ORDER BY idrecord\
                             DESC LIMIT 1;")
        idRecord = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO responses(number,numberofspikes,length,record_idrecord,vpsp,epspFront, epspBack)\
                             VALUES(%s,%s,%s,%s,%s,%s,%s)", (str(rNumber),str(nOfSpikes),str(rLength),str(idRecord),str(vpsp),\
                             str(tmpObject.epspFront),str(tmpObject.epspBack)))
        self.conn.commit()
        cursor.close()

    def dbWriteSpike(self,tmpObject):
        #print("dbWriteSpike")
        ampl=tmpObject.spikeAmpl
        number=tmpObject.spikeNumber
        sLength=tmpObject.spikeMax2-tmpObject.spikeMax1#must be changed to length at 80% or something like that
        maxdiff=tmpObject.spikeMax2Val-tmpObject.spikeMax1Val
        cursor = self.conn.cursor()
        cursor.execute("SELECT idresponses\
                             FROM responses\
                             ORDER BY idresponses\
                             DESC LIMIT 1;")
        idResponse = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO spikes(ampl,number,responses_idresponses,length,maxDiff)\
                             VALUES(%s,%s,%s,%s,%s);", (str(ampl),str(number+1),str(idResponse),str(sLength),str(maxdiff)))
        self.conn.commit()
        cursor.close()
    def dbDisconnect(self):
        self.conn.close()
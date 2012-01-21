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
    
    def __init__(self,filename,substance):
        self.fileName=filename
        self.variables_global()
        self.dbConnect()
        self.substanceName=substance
        
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
            self.userName='filteruser'
            self.userPassword='filter123'
            self.dbName='filterdb'
        self.date= strftime('%Y%m%d',localtime(stat(self.fileName).st_mtime))
        self.numberofresponses=1

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
        cursor = self.conn.cursor()
        cursor.execute("SELECT idexperiment \
                             FROM experiment \
                             ORDER BY idexperiment\
                             DESC LIMIT 1;")
        self.idExperiment = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO record(filename,time,numberofresponses,\
                                            experiment_idexperiment)\
                        VALUES(%s,%s,%s,%s);", (self.fileName,self.time,str(self.numberofresponses),str(self.idExperiment)))
        self.conn.commit()
        cursor.close()

    def dbWriteResponse(self):
        cursor = self.conn.cursor()
        cursor.execute("SELECT idrecord \
                             FROM record \
                             ORDER BY idrecord\
                             DESC LIMIT 1;")
        self.idRecord = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO responses(number,numberofspikes,length,record_idrecord)\
                             VALUES(%s,%s,%s,%s)", (str(self.responsNumber+1),str(self.numberofspikes),str(self.responslength),str(self.idRecord)))
        self.conn.commit()
        cursor.close() 

    def dbWriteSpike(self):
        cursor = self.conn.cursor()
        cursor.execute("SELECT idresponses\
                             FROM responses\
                             ORDER BY idresponses\
                             DESC LIMIT 1;")
        self.idResponse = cursor.fetchall()[0][0]
        cursor.execute("INSERT INTO spikes(ampl,number,responses_idresponses)\
                             VALUES(%s,%s,%s);", (str(self.ampl),str(self.number+1),str(self.idResponse)))
        self.conn.commit()
        cursor.close()
    def dbDisconnect(self):
        self.conn.close()
#!/usr/bin/python2
# -*- coding: utf-8 -*-
#signal filtering and spike finding
"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
'''
Created on 23.11.2011
@author: Malakhin Ilya (pilat1988@gmail.com)

usage: "python2 fEPSP-analyser.py $arg1 $arg2 $arg3 $arg4 $arg5"
    where:
        arg1=path to dir with raw data
        arg2=record frequency in Hz
        arg3=path to dir for saved graphs images
        arg4=smoothing coefficient. 1 is a good choice here
        arg5=name of the experiment. This name will accompany calculated fEPSP amplitudes in mysql db.
    
    output: fEPSP amplitudes writen in mysqldb.
            pickled representations of filtered and smoothed records in $arg3/test/ for \
                using with graph.py
            saved graph images of every finished input data file in ./results/$input_folder_name/ 
'''
from os import stat
from time import strftime,localtime
from glob import glob
import sys,pickle
from root.dbAccess import Mysql_writer
from root.filtering_lib2 import dataSample
import shutil
import codecs


class fepspAnalyser:
    
    def errorProcessing(self,filename,errorType):
        try:
            print "copy file with trouble to separate directory"
            if errorType=="soft":
                shutil.copy2(filename,"./softErrors/")
            else:
                shutil.copy2(filename,"./hardErrors/")
        except:
            print "Unexpected error during copy:", sys.exc_info()
            raise
        
    def __init__(self,arguments):
        print(arguments)
        UTF8Writer = codecs.getwriter('utf8')
        sys.stdout = UTF8Writer(sys.stdout)
	#sys.setdefaultencoding('utf-8')

        #print (sys.stdout.encoding)
	dirPath = unicode(arguments[1])
        print(dirPath)
        fileList=glob(dirPath+'/*.dat')
        print(fileList)
        if int(arguments[7]):
            mysql_writer=Mysql_writer(fileList[0],arguments[5])
            print("write to database enabled")    
            mysql_writer.dbWriteExperiment()
            print("experiment has been written")
            mysql_writer.tagWriter()
            print("tag has been written")
        else:
            mysql_writer="pass"
        for i in fileList:
            fileName=i.split('/')[-1]
            #print(unicode.encode(fileName,'utf-8'))
            try:
                creatingTime=strftime('%H%M%S',localtime(stat(unicode.encode(i,'utf-8')).st_mtime))
                print(creatingTime)
                if int(arguments[7]):
                    print("mysql in main body start")
                    mysql_writer.filePath=i
                    mysql_writer.time = creatingTime
                    mysql_writer.dbWriteRecord()
                    print("recodr has been written")
                    mysql_writer.dbWriteRecordTags(fileName)
                    print("mysql in main body end")
                dataSample1=dataSample(i,mysql_writer,arguments)
                print("dataSample object created")
                dataSample1.dataProcessing()
                print("dataSample processing finished")
                if dataSample1.hardError==1:
                    self.errorProcessing(i, "hard")
                elif dataSample1.softError==1:
                    self.errorProcessing(i, "soft")
                else:
                    try:
                        name=dirPath+u"/"+unicode(creatingTime)+u"_"+fileName+u".pic"
                        with open(name,"w") as fd:
                            pickle.dump(dataSample1.result[dataSample1.stimuli[0][0]:],fd)
                    except: 
                        print "Unexpected error in pickle:", sys.exc_info()
                        raise
                del dataSample1
            except:
                print "Unexpected error in main body:", sys.exc_info()
                try:
                    del dataSample1
                except:
                    pass
        if int(arguments[7]):
            try:
                mysql_writer.dbDisconnect()
            except:
                print "Unexpected error in dbDisconect:", sys.exc_info()
        
if __name__ == "__main__":
    analyserObject=fepspAnalyser(sys.argv)

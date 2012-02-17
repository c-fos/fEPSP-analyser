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
@author: pilat (pilat1988@gmail.com)

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
from dbAccess import Mysql_writer
from filtering_lib2 import dataSample
import shutil



class fepspAnalyser:
    def __init__(self,arguments):
        print(arguments)
        dirPath = str(arguments[1])
        fileList=glob(dirPath+"/*.dat")
        #print(fileList)
        if int(arguments[7]):
            mysql_writer=Mysql_writer(fileList[0],arguments[5])
            print("write to database enabled")    
            mysql_writer.dbWriteExperiment()
            mysql_writer.tagWriter()
        else:
            mysql_writer="pass"
        for i in fileList:
            fileName=i.split('/')[-1]
            try:
                creatingTime=strftime('%H%M%S',localtime(stat(i).st_mtime))
                if int(arguments[7]):
                    mysql_writer.filePath=i
                    mysql_writer.time = creatingTime
                    mysql_writer.dbWriteRecord()
                dataSample1=dataSample(i,mysql_writer,arguments)
                dataSample1.dataProcessing()
                if dataSample1.errorState==1:
                    try:
                        print "copy file with trouble to separate directory"
                        shutil.copy2(i,"./mustChecked/")
                    except:
                        print "Unexpected error during copy:", sys.exc_info()
                        raise
                else:
                    try:
                        name=dirPath+"/"+str(creatingTime)+"_"+fileName+".pic"
                        with open(name,"w") as fd:
                            pickle.dump(dataSample1.result[dataSample1.stimuli[0][0]:],fd)
                    except: 
                        print "Unexpected error in pickle:", sys.exc_info()
                        raise
                del dataSample1
            except:
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
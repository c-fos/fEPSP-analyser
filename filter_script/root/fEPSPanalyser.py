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

class fepspAnalyser:
    def __init__(self,arguments):
        print(arguments)
        dirName = str(arguments[1])
        fileList=glob("./"+dirName+"/*.dat")
        mysql_writer=Mysql_writer(str(fileList[0]),arguments[5])
        if int(arguments[7]):
            print("write to database enabled")    
            mysql_writer.dbWriteExperiment()
        for i in range(len(fileList)):
            try:
                mysql_writer.fileName=str(fileList[i])
                mysql_writer.time= strftime('%H%M%S',localtime(stat(fileList[i]).st_mtime))
                if int(arguments[7]):
                    mysql_writer.dbWriteRecord()
                dataSample1=dataSample(str(fileList[i]),mysql_writer,arguments)
                try:
                    name="./"+dirName+"/"+str(mysql_writer.time)+".pic"
                    with open(name,"w") as fd:
                        pickle.dump(dataSample1.result[dataSample1.stimuli[0][0]:],fd)
                except: 
                    print("pickle with error")
                del dataSample1
            except:
                try:
                    del dataSample1
                except:
                    pass
        try:
            mysql_writer.dbDisconnect()
        except:
            pass
        
if __name__ == "__main__":
    analyserObject=fepspAnalyser(sys.argv)
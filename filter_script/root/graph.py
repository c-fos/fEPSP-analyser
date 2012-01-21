#!/usr/bin/python2
# -*- coding: utf-8 -*-
#this is fEPSP-analyser addon which plotting 3d reconstruction of experiment
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
Created on 29.11.2011

@author: pilat
'''

from glob import glob
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pickle,os,sys,pylab,numpy

class graphReconstruction:
    def __init__(self, fileList, start=0, stop=10000, rstride=2, cstride=20, debug="0"):
        list1=[]
        for i in fileList:
            if debug=="1":
                print(i) 
            if os.stat(i).st_size>0:
                fd=open(i,'r')
                data=list(pickle.load(fd))
                if debug=="1":
                    print (len(data))
                if len(data)>=stop:
                    list1+=[data[start:stop]]
                fd.close()
        x = numpy.arange(start,stop,1)
        y = numpy.arange(0,len(list1),1)
        xgrid, ygrid = numpy.meshgrid(x, y)
        zgrid = numpy.array(list1)

        fig = pylab.figure()
        axes = Axes3D(fig)
        if debug=="1":
            print(numpy.size(xgrid),numpy.size(ygrid),numpy.size(zgrid))
        axes.plot_surface(xgrid, ygrid, zgrid,rstride=rstride, cstride=cstride, cmap = cm.jet)
        pylab.show()
        
if __name__ == "__main__":       
    dirName = str(sys.argv[1])
    fileList=glob(dirName+"/*.pic")
    fileList.sort()
    start=int(sys.argv[2])
    stop=int(sys.argv[3])
    rstride=int(sys.argv[4])
    cstride=int(sys.argv[5])
    debug=sys.argv[6]
    reconst=graphReconstruction(fileList,start=start,stop=stop,rstride=rstride,cstride=cstride,debug=debug)


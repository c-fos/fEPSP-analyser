#!/usr/bin/python2
# -*- coding: utf-8 -*-
'''
Created on 05.12.2011

@author: pilat
'''

class Spike:
    def __init__(self,frequency):
        #primary variables
        self.responseNumber=0 #int
        self.responseStart=0  #point
        self.responseEnd=0    #point
        self.spikeNumber=0    #int
        self.allSpikes=0      #int
        self.spikeMax1=0      #point
        self.spikeMin=0       #point
        self.spikeMax2=0      #point
        self.spikeAmpl=0      #depend on input
        self.frequency=frequency
        #secondary variables
        self.spikeLength=0    #msec        
        self.spikeFrequency=0 #hZ
    
    def calculate(self):
        #self.spikeLength=1000.0*(self.spikeMax2-self.spikeMax1)/self.frequency#msec
        self.spikeFrequency=(1.0/self.spikeLength)*1000


class Response:
    def __init__(self):
        pass
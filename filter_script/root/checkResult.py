'''
Created on 15.02.2012

@author: pilat

This function is just simple placeholder.
'''
<<<<<<< HEAD
import matplotlib.pyplot as plt
import pywt
import numpy
from externalFunctions import iswt

def resultAnalysis(parentObj,debug):
    if len(parentObj.responseDict)>=1:# and len(parentObj.epsp[0])==len(parentObj.responseDict):
        print "Test passed succsessful"
    else:
        print "Test error"
        parentObj.errorState=1
    diffArray=parentObj.cleanData-parentObj.result
    coeffs=pywt.swt(diffArray, 'haar', level=5)
    for i in range(5):
        cA, cD = coeffs[i]
        cD=numpy.zeros(len(cA),dtype='float32')
        coeffs[i]=cA, cD
    lowFrquDiffArray=iswt(coeffs,'haar')
    if debug==1:    
        fig, axes_list = plt.subplots(2, 1, sharex=True)
        ax = axes_list[0]
        ax.plot(parentObj.cleanData, 'y')
        ax.plot(parentObj.result,'b')    
        ax.grid(color='k', linestyle='-', linewidth=0.4)
        bx = axes_list[1]
        bx.plot(diffArray)
        bx.plot(lowFrquDiffArray,'k')
        bx.grid(color='k', linestyle='-', linewidth=0.4)
        plt.tight_layout()
        tex = str((lowFrquDiffArray.std(),lowFrquDiffArray.ptp(),"diffArray.std(),diffArray.ptp()"))
        bx.text(50,50, tex, fontsize=12, va='bottom')
        plt.show()
        del fig
=======


def resultAnalysis(parentObj):
    if len(parentObj.responseDict)>=1:# and len(parentObj.epsp[0])==len(parentObj.responseDict):
        print "Test passed succsessful"
    else:
        print "Test error"
        parentObj.errorState=1
    
>>>>>>> refs/heads/Stable

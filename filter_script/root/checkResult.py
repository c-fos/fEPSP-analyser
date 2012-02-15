'''
Created on 15.02.2012

@author: pilat
'''
def resultAnalysis(parentObj):
    if len(parentObj.responseDict)>=1 and len(parentObj.epsp[0])==len(parentObj.responseDict):
        print "Test passed succsessful"
    else:
        print "Test error"
        parentObj.errorState=1
    
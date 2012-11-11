'''
Created on 26.09.2012

@author: pilat
'''
import unittest
import filtering_lib2


class Test(unittest.TestCase):


    def testgetLocalPtp(self):
        data=[-5,-5,-5,-5,-5,-5,5,5,5,5,5,5]
        data2=[0,0,0,0,0,0]
        dataSample=filtering_lib2.dataSample("test","pass",[0,0,0,0,0,0,0,0,0])
        dataSample.argReading()
        print(dataSample.getLocalPtp(data,3))
        self.assertTrue(dataSample.getLocalPtp(data,3)==10)
        print(dataSample.getLocalPtp(data2,3))
        self.assertEqual(dataSample.getLocalPtp(data2,3),0)
        print(dataSample.getLocalPtp(data,15))
        self.assertFalse(dataSample.getLocalPtp(data,15)==10)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
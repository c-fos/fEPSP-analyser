'''
Created on 18.02.2012

@author: pilat
'''
import numpy as np
import hcluster as hc
y=[]
x=np.asmatrix([[0,0,1,1,0,0,0,0,0,1,0,0,1,1,2,1,1,0,0,0],[0,0,0,1,1,0,0,0,0,0,0,0,1,2,3,2,1,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,1,1,2,1,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]],dtype='int')
for i in range(len(x)):
    for j in range(x[0].size):
        y.append([x[i,j]]) 
print y
z=hc.fclusterdata(y,1.153)
print(z)
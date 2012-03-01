'''
Created on 18.02.2012

@author: pilat
'''
import numpy as np
import scipy.cluster.hierarchy as hc
import scipy.spatial.distance as dist

y=[]
x=np.asmatrix([[0,0,1,1,0,0,0,0,0,1,0,0,1,1,2,1,1,0,0,0],[0,0,0,1,1,0,0,0,0,0,0,0,1,2,3,2,1,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,1,1,2,1,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]],dtype='int')
shape1=x.shape
for i in range(len(x)):
    for j in range(x[0].size):
        y.append([i*x[i,j],j*x[i,j]]) 
print y
dist1=dist.pdist(np.array([np.array(y[i]) for i in range(len(y))]))
z1=hc.complete(dist1)
print(dist1)
z=hc.fcluster(z1,1.1)
print(x,np.reshape(z,shape1))
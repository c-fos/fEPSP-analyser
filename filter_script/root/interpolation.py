'''
Created on 09.02.2012

@author: pilat

import numpy as np
from scipy.interpolate import interp1d
import pylab as pl

measured_time = np.array([0,10,55,70,78,100])#np.linspace(0, 1, 10)
measures = np.array([1,1.1,2.3,1.7,1.6,1.4])

computed_time = np.linspace(0, 100, 1000)
cubic_interp = interp1d(measured_time, measures, kind='cubic')
cubic_results = cubic_interp(computed_time)
print(cubic_results)
pl.plot(measured_time, measures, 'o', ms=6, label='measures')
pl.plot(computed_time, cubic_results, label='cubic interp')
pl.legend()
pl.show()
'''
import numpy as np
import pylab as pl
import scipy
import scipy.interpolate
x = np.array([0,1,2,3,4,10,55,70,78,100,101,102,103])#np.linspace(0, 1, 10)
y = np.array([1,1,1.1,1.1,1.1,1.1,2.3,1.7,1.6,1.43,1.41,1.43,1.41])
sp = scipy.interpolate.spline(x,y, np.array(range(103)),kind='smoothest')
print(sp) # should be a good approximation of cos(0.5)=0.8775
pl.plot(x,y,'o',ms=6)
pl.plot(np.array(range(103)),sp)
pl.show()
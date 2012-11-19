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
#import pylab as pl
#from scipy import polyval, polyfit
#import scipy.interpolate as inter
x = np.array([1,2,3,4,10,55,70,78,100,101,102,103])#np.linspace(0, 1, 10)
y = np.array([1,1.05,1.1,1.1,1.1,2.3,1.7,1.6,1.43,1.41,1.43,1.41])
#sp = inter.spline(x,y, np.array(range(103)),kind='smoothest')
x2=np.array(range(103))
#f=inter.Rbf(x,y,smooth=0.0005)
#f2=inter.Rbf(x,y,smooth=1,function='thin_plate')
#ar=polyfit(x[5:],y[5:],2)
#xr=polyval(ar,x[5:])
'''

        'multiquadric': sqrt((r/self.epsilon)**2 + 1)
        'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
        'gaussian': exp(-(r/self.epsilon)**2)
        'linear': r
        'cubic': r**3
        'quintic': r**5
        'thin_plate': r**2 * log(r)
'''

#y2=f(x2)
#y3=f2(x2)
#print(sp) # should be a good approximation of cos(0.5)=0.8775
#pl.plot(x,y,'o',ms=6)
#pl.plot(np.array(range(103)),sp,'r')
#pl.plot(x2,y2,'g')
#pl.plot(x2,y3,'y')
#pl.plot(x[5:],xr)
#pl.show()
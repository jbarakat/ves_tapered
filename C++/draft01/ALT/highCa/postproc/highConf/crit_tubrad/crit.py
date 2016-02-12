from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from pylab import *
import numpy as np
from scipy import interpolate, interpolate
from math import log10, floor, pi
#from cmath import sqrt
import cmath

n = 10000
vspan = np.linspace(0.6,1.0,n)
aspan = np.linspace(0.6,1.0,n)

# initialize plots
fs = 10
plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

cf = (cmath.sqrt(3)/2.0)*1j
for i in range(len(vspan)) :
	v = vspan[i]

	a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
	a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
	a = -0.5*(a1 + a2) + cf*(a1 - a2)
	a = a.real

	if abs(v - 1.0) < 1e-12 :
		a = 1.0

	aspan[i] = a

# red blood cell
v = 0.61
a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
a = -0.5*(a1 + a2) + cf*(a1 - a2)
a = a.real

plt.plot(vspan,aspan,'k')
plt.plot(v,a,'k.')
plt.text(v+0.01, a-0.02, 'RBC',fontsize=fs)
plt.xlim((min(vspan),max(vspan)))
plt.xlabel('$v$',fontsize=fs)
plt.ylabel('$a_c$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show()

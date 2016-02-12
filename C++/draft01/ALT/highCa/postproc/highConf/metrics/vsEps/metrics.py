def load_src(name, fpath):
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from pylab import *
import numpy as np
from scipy import interpolate, interpolate
from math import log10, floor, pi
import cmath
from readMetrics import *

# initialize plots
fs = 10
plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

# array of reduced volume
vmin = 65
vmax = 100
dv   = 5
N = (vmax - vmin)/dv + 1
vspan = linspace(vmin,vmax,N)
V = []
for i in range(len(vspan)) :
	V.append(str(int(vspan[i])))
V = np.column_stack(V)
V = V[0]
print V

# array of confinement
cmin = 991
cmax = 999
dc   = 1
N = (cmax - cmin)/dc + 1
cspan = linspace(cmin,cmax,N)
C = ['98','99']
for i in range(len(cspan)) :
	C.append(str(int(cspan[i])))
C = np.column_stack(C)
C = C[0]
print C

for strv in V :
	PDROP  = []
	FLUX   = []
	VELRAT = []
	VOLRAT = []
	FILM   = []
	for strc in C :
		# read data
		(redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau) = readMetrics(strv,strc)
		
		v = redvol
		a = tubrad
	
		# film thickness (epsilon)
		n  = len(r)
		nh = n/2
		h  = 1.0 - r[nh]
		
		# scale pressure and surface tension
		cf = pow(h,1.5)
		for i in range(len(p)) :
			p  [i] = cf*p  [i]
			gam[i] = cf*gam[i]
		
		# scale pressure drop and flux
		cf = pow(h,1.0)
		pdrop = cf*pdrop
		flux  = flux/cf
	
		# calculate critical tube radius
		a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
		a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
		cf = (cmath.sqrt(3)/2.0)*1j
		
		ac = -0.5*(a1 + a2) + cf*(a1 - a2)
		ac = float(ac.real)
		ac = round(ac,3)
		
		frac = round(a/ac,2)
	#	h = round(h,2)
	
		PDROP .append(pdrop)
		FLUX  .append(flux)
		VOLRAT.append(redvol)
		FILM  .append(h)
	
	PDROP  = np.column_stack(PDROP )
	FLUX   = np.column_stack(FLUX  )
	VOLRAT = np.column_stack(VOLRAT)
	FILM   = np.column_stack(FILM  )
	PDROP  = PDROP [0]
	FLUX   = FLUX  [0]
	VOLRAT = VOLRAT[0]
	FILM   = FILM  [0]

	# calculate critical tube radius
	v = VOLRAT[-1]
	
	a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
	a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
	
	ac = -0.5*(a1 + a2) + cf*(a1 - a2)
	ac = float(ac.real)
	ac = round(ac,3)
	
	# calculate length of cylindrical section (for a --> ac)
	L = 2.0*(pow(ac*ac,-1.0) - 1.0)
	
	N = 10000
	espan = linspace(0.0,0.05,N)
	pspan = linspace(0.0,0.05,N)
	cf = (cmath.sqrt(3)/2.0)*1j
	for i in range(len(espan)) :
		eps = espan[i]
	
		dp = 2.0*L + 17.8*pow(eps,0.5)
		pspan[i] = dp
	
	plt.plot(FILM,PDROP,'k.')
	plt.plot(espan,pspan,'k--')

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlim((0.0,0.01))
plt.xlabel(r'$\epsilon$')
plt.ylabel(r'$\epsilon a {\it \Delta} p$')

plt.show()

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
vmax = 99
dv   = 1
N = (vmax - vmin)/dv + 1
vspan = linspace(vmin,vmax,N)
V = []
for i in range(len(vspan)) :
	V.append(str(int(vspan[i])))
V = np.column_stack(V)
V = V[0]

PDROP  = []
PDROP2 = []
FLUX   = []
VELRAT = []
VOLRAT = []
FILM   = []
FRONTTENS = []
for strv in V :
	# read data
	(redvol, tubrad, flux, pdrop, gamf, s, x, r, p, gam, tau) = readMetrics(strv)
	
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
	gamf = cf*gamf
	
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
	PDROP2.append((p[0]-p[-1])*pow(h,-0.5))
	FLUX  .append(flux)
	VOLRAT.append(redvol)
	FILM  .append(h)
	FRONTTENS.append(gamf)

PDROP  = np.column_stack(PDROP )
PDROP2 = np.column_stack(PDROP2 )
FLUX   = np.column_stack(FLUX  )
VOLRAT = np.column_stack(VOLRAT)
FILM   = np.column_stack(FILM  )
FRONTTENS = np.column_stack(FRONTTENS)
PDROP  = PDROP [0]
PDROP2 = PDROP2[0]
FLUX   = FLUX  [0]
VOLRAT = VOLRAT[0]
FILM   = FILM  [0]
FRONTTENS = FRONTTENS[0]

plt.plot(VOLRAT,PDROP2,'wo')
plt.plot(VOLRAT,PDROP,'k.')

N = 10000
vspan = linspace(0.65,0.99999,N)
pspan = linspace(0.65,0.99999,N)
qspan = linspace(0.65,0.99999,N)
cf = (cmath.sqrt(3)/2.0)*1j
eps = np.mean(FILM)
for i in range(len(vspan)) :
	v = vspan[i]
	# calculate critical tube radius
	a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
	a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
	
	ac = -0.5*(a1 + a2) + cf*(a1 - a2)
	ac = float(ac.real)
	ac = round(ac,3)

	L = 2.0*(pow(ac*ac,-1.0) - 1.0)

	dp = 2.0*L #+ 17.8*pow(eps,0.5)
	pspan[i] = dp

	q  = 1.0 - eps/3.0 - eps*eps/12.0
	qspan[i] = q

#plt.plot(vspan,pspan,'k-')

for i in range(len(vspan)) :
	pspan[i] = pspan[i] + 17.8*pow(eps,0.5)
	#pspan[i] = pspan[i] + 20*pow(eps,0.5)

plt.plot(vspan,pspan,'k-')

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel(r'$v$')
plt.ylabel(r'$\epsilon a {\it \Delta} p$')
leg = plt.legend(('numerical solution','asymptotic solution'),
loc='best',fontsize=fs)
leg.get_frame().set_linewidth(0.0)
print FILM

#plt.plot(VOLRAT,FLUX,'k.')
#plt.plot(vspan,qspan,'k--')
#plt.xticks(fontsize=fs)
#plt.yticks(fontsize=fs)
#plt.xlabel(r'$v$')
#plt.ylabel(r'$(\epsilon a)^{-1} Q$')

#plt.plot(VOLRAT,FLUX,'k.')
#plt.xticks(fontsize=fs)
#plt.yticks(fontsize=fs)
#plt.xlabel(r'$v$')
#plt.ylabel(r'$\epsilon$')

#plt.plot(VOLRAT,FRONTTENS,'k.')
#plt.xticks(fontsize=fs)
#plt.yticks(fontsize=fs)
#plt.xlabel(r'$v$')
#plt.ylabel(r'$\epsilon^{\frac{3}{2}}\tau^f$')
#print FILM
#print FRONTTENS

plt.show()

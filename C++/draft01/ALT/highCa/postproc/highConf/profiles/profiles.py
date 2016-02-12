def load_src(name, fpath):
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from pylab import *
import numpy as np
from scipy import interpolate, integrate
from math import log10, floor, pi
import cmath
from readProfiles import readProfiles

# initialize plots
fs = 10
#plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(6,6))
#ax1 = plt.subplot(2,1,1)
#ax2 = ax1.twinx()
#ax3 = plt.subplot(2,1,2, sharex=ax1)
ax1 = plt.subplot(3,1,1)
ax2 = plt.subplot(3,1,2, sharex=ax1)
ax3 = plt.subplot(3,1,3, sharex=ax2)

fig.subplots_adjust(hspace=0)


# array of reduced volume
V = ['90','85','80','75','70','65']
#V = ['99', '98','97','96','95','94','93','92','91','90']

for strv in V :
	
	# read data
	(redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau, tauw) = readProfiles(strv)
	
	v = redvol
	a = tubrad
	
	# scale pressure and surface tension
	n  = len(r)
	nh = n/2
	h  = 1.0 - r[nh]

	h = 0.0
	for i in [nh-2,nh-1, nh, nh+1, nh+2] :
		hh = 1.0 - r[i]
		h = h + hh
	h = h/5.0

	cf = pow(h,1.5)
	for i in range(len(p)) :
		p  [i] = cf*p  [i]
		gam[i] = cf*gam[i]
	#h = round(h,2)

	
	# integrate 1/h over domain
	kernel = np.linspace(0,1,len(s))
	for i in range(len(s)) :
		hh = a*(1.0 - r[i])
		if s[i] > 0.5*pi and s[i] < s[-1] - 0.5*pi :
			kernel[i] = (4.0*hh - 3.0*h*a)/(hh*hh)
		else :
			kernel[i] = 0.0
		kernel[i] = (4.0*hh - 3.0*h*a)/(hh*hh)
		kernel[i] = -tauw[i]
		kernel[i] = -tau[i]
	
	#	hh = 1.0 - r[i]
	#	kernel[i] = h/hh
	#	if s[i] > 0.5*pi and s[i] < s[-1] - 0.5*pi :
	#		hh = 1.0 - r[i]
	#		kernel[i] = h/hh
	#	else :
	#		kernel[i] = 0.0
	length = 2.0*(pow(a, -2.0) - 1.0)
	
	cf = pow(h,0.5)
	integral = integrate.trapz(kernel,s,0)
	integral = integral - length/h
	integral = integral*cf
	print integral

	#print pow(h,0.5)*length
	#print gam[-1] - gam[0]
	#print gam[-1] - gam[0] - (pow(h,0.5)*length)

	
	# scale radial and axial coordinates
	minx = max(x)
	for i in range(len(r)) :
		x[i] = x[i] - minx
	#	r[i] = a*r[i]
	#	x[i] = a*x[i]
	
	# reflect shape across symmetry axis
	xx  = np.concatenate((x, np.flipud( x)))
	rr  = np.concatenate((r, np.flipud(-r)))
	
	# calculate critical tube radius
	a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
	a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
	cf = (cmath.sqrt(3)/2.0)*1j
	
	ac = -0.5*(a1 + a2) + cf*(a1 - a2)
	ac = float(ac.real)
	ac = round(ac,3)
	
	frac = round(a/ac,2)
	
	# wall
	xw = [min(x)-0.1,max(x)+0.1]
	uw = [a,a]
	lw = [-a,-a]
	uw = [1,1]
	lw = [-1,-1]
	
	#fig = plt.subplots(figsize=(12,12))
	ax1.plot(xx,rr, 'k-')
	ax2.plot(x,  p, 'k-')
	ax3.plot(x,gam, 'k-')
#	ax3.plot(xw,uw, 'k--')
#	ax3.plot(xw,lw, 'k--')
	
	# plot preferences
	if strv == V[-1] :
		ax1.text(xx[0]+0.1, rr[0], '$v$ = 0.' + str(int(v*100)),fontsize=fs)
	else :
		ax1.text(xx[0]+0.1, rr[0], '0.' + str(int(v*100)),fontsize=fs)
	
	ax1.set_ylabel(r'$r/a$',fontsize=fs)
	ax2.set_ylabel(r'$\epsilon^{\frac{3}{2}}a p$', color='k',fontsize=fs)
	ax3.set_ylabel(r'$\epsilon^{\frac{3}{2}}\tau$', color='k',fontsize=fs)
	ax3.set_xlabel(r'$x/a$',fontsize=fs)
	
	setp(ax1.get_xticklabels(), visible=False)
	setp(ax1.get_yticklabels(), fontsize=fs)
	setp(ax2.get_yticklabels(), fontsize=fs)
	setp(ax3.get_yticklabels(), fontsize=fs)
	setp(ax3.get_xticklabels(), fontsize=fs)
	#ax1.set_title('$v$ =' + str(v) + ', $a$ = ' + str(frac) + '$a_c$' + ', $\epsilon$ = ' + str(h),fontsize=fs)
	
	#plt.savefig(homedir + '/../plots/profiles/shape/v' + redvol + '/rmax' + maxrad + '/shape_v' + redvol + '_rmax' + maxrad + '_Ca' + capnum + '.pdf')
	#plt.close()
ax1.set_ylim((min(rr )-0.2,max(rr )+0.2))
ax2.set_ylim((min(p  )-1.5,max(p  )+0.5))
ax3.set_ylim((min(gam)-0.12,max(gam)+0.12))
ax3.set_xlim((min(x)-0.1,max(x)+0.1))
ax1.set_yticklabels([str(abs(i)) for i in ax1.get_yticks()])
ax2.set_yticklabels([str(int(i)) for i in ax2.get_yticks()])
ax3.set_yticklabels([str(i) for i in ax3.get_yticks()])
ax3.set_xticklabels([str(int(i)) for i in ax3.get_xticks()])
#ax3.set_ylim(-a-0.02,a+0.02)
plt.show()

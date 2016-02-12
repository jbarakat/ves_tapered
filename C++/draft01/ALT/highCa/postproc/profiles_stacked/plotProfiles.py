def load_src(name, fpath):
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def plotProfiles(redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau) :
	from matplotlib import pyplot as plt
	from matplotlib.pyplot import cm
	from pylab import *
	import numpy as np
	from scipy import interpolate, interpolate
	from math import log10, floor, pi
	import cmath

	v = redvol
	a = tubrad

	# scale pressure and surface tension
	n  = len(r)
	nh = n/2
	h  = 1.0 - r[nh]
	cf = pow(h,1.5)
	for i in range(len(p)) :
		p  [i] = cf*p  [i]/a
		gam[i] = cf*gam[i]
	h = round(h,2)

	# scale radial and axial coordinates
	for i in range(len(r)) :
		r[i] = a*r[i]
		x[i] = a*x[i]

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
	#wall = []
	#for i in range(len(x)) :
	#	wall.append(1.0)
	#wall = np.column_stack(wall)
	#wall = wall[0]
	
	#fig = plt.subplots(figsize=(12,12))
	fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True, figsize=(12,12))
	ax1 = plt.subplot(2,1,1)
	ax2 = ax1.twinx()
	ax3 = plt.subplot(2,1,2, sharex=ax1)
	fig.subplots_adjust(hspace=0)
	
	ax1.plot(x,  p, 'b-')
	ax2.plot(x,gam, 'g-')
	ax3.plot(xx,rr, 'k-')
	ax3.plot(xw,uw, 'k--')
	ax3.plot(xw,lw, 'k--')
	
	fs = 24
	ax1.set_ylabel(r'$\epsilon^{\frac{3}{2}}p$', color='b',fontsize=fs)
	ax2.set_ylabel(r'$\epsilon^{\frac{3}{2}}\tau$', color='g',fontsize=fs)
	ax3.set_ylabel(r'$r$',fontsize=fs)
	ax3.set_xlabel(r'$x$',fontsize=fs)
	
	setp(ax1.get_xticklabels(), visible=False)
	setp(ax1.get_yticklabels(), color='b',fontsize=fs-4)
	setp(ax2.get_yticklabels(), color='g',fontsize=fs-4)
	setp(ax3.get_yticklabels(), fontsize=fs-4)
	setp(ax3.get_xticklabels(), fontsize=fs-4)
	ax3.set_xlim((min(x)-0.1,max(x)+0.1))
	ax3.set_ylim(-a-0.02,a+0.02)
	ax3.set_yticklabels([str(abs(i)) for i in ax3.get_yticks()])
	#ax1.set_title('$v$ =' + str(v) + ', $a_c$ = ' + str(ac) + ', $a$ = ' + str(frac) + '$a_c$' + ', $\epsilon$ = ' + str(h),fontsize=fs)
	ax1.set_title('$v$ =' + str(v) + ', $a$ = ' + str(frac) + '$a_c$' + ', $\epsilon$ = ' + str(h),fontsize=fs)
	
	#ax1.set_ylim((-120,20))
	#ax2.set_ylim((-5,40))
	
	#plt.savefig(homedir + '/../plots/profiles/shape/v' + redvol + '/rmax' + maxrad + '/shape_v' + redvol + '_rmax' + maxrad + '_Ca' + capnum + '.pdf')
	#plt.close()
	plt.show()

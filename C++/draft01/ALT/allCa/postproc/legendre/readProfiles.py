def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readProfiles() :
	import numpy as np
	from scipy import interpolate, integrate
	from math import log10, floor, pi, sqrt, atan, acos, cos, sin
	import os.path
	
	altdir = os.path.abspath("../")
	
	load_src("readLT", altdir + "/tools/readLT.py")
	from readLT import *
	
	###################################################
	#                     SETUP                       #
	###################################################
	
	# read parameters
	ifile = open('./params.in', 'r')
	print ifile
	for line in ifile :
		params = line.split()		# space delimited text
		if params[0] == 'PROFILE' :
			STR = params[2]
		if params[0] == 'REDVOL' :
			redvol = params[2]
		if params[0] == 'CONFIN' :
			confin = params[2]
		if params[0] == 'CAPNUM' :
			capnum = params[2]
	ifile.close()
		
	###################################################
	#               LUBRICATION THEORY                #
	###################################################
		
	xL = []
	yL = []

	vL   = 0
	epsL = 0
	CaL  = 0
	
	# load data file
	(Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLT(altdir, redvol, confin, capnum)

	# reflect shape across symmetry axis
	xx  = np.concatenate((x, np.flipud( x)))
	rr  = np.concatenate((r, np.flipud(-r)))

	# determine polar angle and spherical radius
	thet = []
	rho  = []
	for i in range(len(x)) :
		xi = x[i]
		ri = r[i]
		rhoi = sqrt(xi*xi + ri*ri)
		#theti = atan(ri/xi)
		theti = acos(xi/rhoi)
	
		thet.append(theti)
		rho.append(rhoi)
	thet = np.column_stack(thet)
	rho  = np.column_stack(rho)
	thet = thet[0]
	rho  = rho[0] 
	
	# calculate coefficients for Legendre representation
	N = 7
	cf = []
	for n in range(N) :
		fint = []
		for i in range(len(x)) : 
			rhoi  = rho[i]
			theti = thet[i]
			costh = cos(theti)
			costh2 = costh*costh
			costh4 = costh2*costh2
			costh6 = costh4*costh2
			if n == 0 :
				Pni = 1.0
			if n == 1 :
				Pni = costh
			if n == 2 :
				Pni = 0.5*(3.0*costh2 - 1.0)
			if n == 3 :
				Pni = 0.5*(5.0*costh2 - 3.0)*costh
			if n == 4 :
				Pni = 0.125*(35.0*costh4 - 30.0*costh2 + 3.0)
			if n == 5 :
				Pni = 0.125*(63.0*costh4 - 70.0*costh2 + 15.0)*costh
			if n == 6 :
				Pni = 0.0625*(231.0*costh6 - 315.0*costh4 + 105.0*costh2 - 5.0)

			finti = rhoi*Pni*sin(theti)
			fint.append(finti)
		
		fint = np.column_stack(fint)
		fint = fint[0]

		cfn = -0.5*(2*n + 1.0)*integrate.trapz(fint,thet)
		cf.append(cfn)
	
	cf = np.column_stack(cf)
	cf = cf[0]

	rhoproj = []
	xproj = []
	rproj = []
	thetproj = np.linspace(0,pi,1000)
	for i in range(len(thetproj)) :
		theti = thetproj[i]
		costh = cos(theti)
		costh2 = costh*costh
		costh4 = costh2*costh2
		costh6 = costh4*costh2

		rhoi = 0.0
		N = 5
		for n in range(N) :
			if n == 0 :
				Pni = 1.0
			if n == 1 :
				Pni = costh
			if n == 2 :
				Pni = 0.5*(3.0*costh2 - 1.0)
			if n == 3 :
				Pni = 0.5*(5.0*costh2 - 3.0)*costh
			if n == 4 :
				Pni = 0.125*(35.0*costh4 - 30.0*costh2 + 3.0)
			if n == 5 :
				Pni = 0.125*(63.0*costh4 - 70.0*costh2 + 15.0)*costh
			if n == 6 :
				Pni = 0.0625*(231.0*costh6 - 315.0*costh4 + 105.0*costh2 - 5.0)

			rhoi = rhoi + cf[n]*Pni

		xi   = rhoi*cos(theti)
		ri   = rhoi*sin(theti)
		
		xproj.append(xi)
		rproj.append(ri)
	
	xproj = np.column_stack(xproj)
	rproj = np.column_stack(rproj)
	xproj = xproj[0]
	rproj = rproj[0]

	print cf


	from matplotlib import pyplot as plt
	plt.figure()
	plt.plot(x,r,xproj,rproj)
	plt.show()
	
#	plt.figure()
##	plt.plot(thet, rho,'.')
#	plt.plot(thet, rho)
#	plt.show()
	
	## interpolate to xi using splines
	#tck = interpolate.splrep(x, y, s=0)
	#yi = interpolate.splev(xi, tck, der=0)

	if   (STR == 'SHAPE') :
		xL = xx
		yL = rr
	elif (STR == 'PRESS') :
		xL = x
		yL = p
	elif (STR == 'SHEAR') :
		xL = x
		yL = tau
	elif (STR == 'TENS' ) :
		xL = x
		yL = gam
	elif (STR == 'TRANS' ) :
		xL = x
		yL = qs
	
	# add a decimal point to the Ca string if need be
	if (capnum[0] == '0') :
		capnum = capnum[0] + '.' + capnum[1:]
	
	a  = sqrt(area/(4.0*pi))
	a3 = a*a*a
	v  = vlme/((4.0/3.0)*pi*a3)
	v  = round(v, 2)

	vL   = v
	epsL = eps
	CaL  = Ca

	return (STR, vL, epsL, CaL, xL, yL)

def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readProfiles() :
	import numpy as np
	#from scipy import interpolate, interpolate
	from math import log10, floor, pi, sqrt, log
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

	tauw = []
	for i in range(len(s)) :
		rr = r[i]
		r2 = rr*rr
		if (rr > 0.0001) :
			from math import log
			logr = log(rr)
		
			g =   (8.0/(1.0 - r2))
			g = g*(eps - 1.0 - (1.0 - r2)/(2.0*logr))
			g = g/(1.0 + r2 + (1.0 - r2)/logr)

			e = 0.25*g*(2.0 + (1.0 - r2)/logr) + 1.0/logr
		else : 
			g  = 8.0*(eps - 1.0) # limiting value as r --> 0
			e = 0.5*g
		#tauw.append(-e*tx[i])
		tauw.append(e)
	tauw = np.column_stack(tauw)
	tauw = tauw[0]

	lengthx = len(x)
	mid = lengthx/2
#	print mid
#	print (1.0 - r[mid])/eps
#	print 1.0 - ((1.0 - r[mid])/eps - 1.0)*6.0/eps

	rmid = r[mid]
	Q2 = eps
#	print eps - (1.0 - rmid)
	rmid = 0.990147
	Q2 = 1.0 - rmid - (1.0/3.0)*pow((1.0 - rmid),2.0) + (1.0/4.0)*pow((1.0 - rmid),3.0)
#	Q2 = 0.00983
#	print 0.009821 - (1.0 - rmid + (1.0/4.0)*pow((1.0 - rmid),3.0))
	print Q2
	r2 = rmid*rmid
	g = (8.0/(1.0 - r2))
	g = g*(Q2 - 1.0 - (1.0 - r2)/(2.0*log(rmid)))
	g = g/(1.0 + r2 + (1.0 - r2)/log(rmid))

	print rmid
	print g

	# reflect shape across symmetry axis
	xx  = np.concatenate((x, np.flipud( x)))
	rr  = np.concatenate((r, np.flipud(-r)))
	
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
	elif (STR == 'WALLSHEAR') :
		xL = x
		yL = tauw
	elif (STR == 'TENS' ) :
		xL = x
		yL = gam
	#	for i in range(len(yL)):
	#		yL[i] = -p[i]*r[i]/gam[i]
	elif (STR == 'CURV' ) :
		xL = x
		yL = cs
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
	
	###################################################
	#    HIGH CAPILLARY NUMBER LUBRICATION THEORY     #
	###################################################
		
	xH = []
	yH = []

	# load data file
	(eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLTHighCa(altdir, redvol, confin)

	tauw = []
	for i in range(len(s)) :
		rr = r[i]
		r2 = rr*rr
		if (rr > 0.0001) :
			from math import log
			logr = log(rr)
		
			g =   (8.0/(1.0 - r2))
			g = g*(eps - 1.0 - (1.0 - r2)/(2.0*logr))
			g = g/(1.0 + r2 + (1.0 - r2)/logr)

			e = 0.25*g*(2.0 + (1.0 - r2)/logr) + 1.0/logr
		else : 
			g  = 8.0*(eps - 1.0) # limiting value as r --> 0
			e = 0.5*g
		#tauw.append(-e*tx[i])
		tauw.append(e)
	tauw = np.column_stack(tauw)
	tauw = tauw[0]

	# reflect shape across symmetry axis
	xx  = np.concatenate((x, np.flipud( x)))
	rr  = np.concatenate((r, np.flipud(-r)))
	
	## interpolate to xi using splines
	#tck = interpolate.splrep(x, y, s=0)
	#yi = interpolate.splev(xi, tck, der=0)

	if   (STR == 'SHAPE') :
		xH = xx
		yH = rr
	elif (STR == 'PRESS') :
		xH = x
		yH = p
	elif (STR == 'SHEAR') :
		xH = x
		yH = tau
	#	for i in range(len(tau)) :
	#		yH[i] = -tau[i]#/(Ca)
	elif (STR == 'WALLSHEAR') :
		xH = x
		yH = tauw
	elif (STR == 'TENS' ) :
		xH = x
		yH = gam
	#	for i in range(len(yH)):
	#		yH[i] = -p[i]*r[i]/gam[i]
	elif (STR == 'CURV' ) :
		xH = x
		yH = cs
	elif (STR == 'TRANS' ) :
		xH = x
		yH = qs
	
	return (STR, vL, epsL, CaL, xL, yL, xH, yH)

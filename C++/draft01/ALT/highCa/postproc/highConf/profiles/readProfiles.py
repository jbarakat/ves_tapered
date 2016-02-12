def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readProfiles(strv) :
	import numpy as np
	#from scipy import interpolate, interpolate
	from math import log10, floor, pi, sqrt, log
	import os.path
	
	altdir = os.path.abspath("../../")
	
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
	redvol = strv
		
	###################################################
	#    HIGH CAPILLARY NUMBER LUBRICATION THEORY     #
	###################################################
		
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
	
	a  = sqrt(area/(4.0*pi))
	a3 = a*a*a
	v  = vlme/((4.0/3.0)*pi*a3)

	tubrad = round(1.0/a,3)
	redvol = round(v,2)
	flux   = round(eps,3)
	pdrop  = p[0] - p[-1]
	
	return (redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau, tauw)

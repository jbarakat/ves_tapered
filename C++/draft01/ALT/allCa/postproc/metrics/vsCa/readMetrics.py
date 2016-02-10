def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readMetrics(redvol,confin) :
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import cm
	import numpy as np
	from scipy import interpolate, interpolate
	from math import log10, floor, pi, sqrt
	import os.path
	
	homedir = os.path.abspath(os.path.join("../", os.pardir))

	load_src("readLT", homedir + "/tools/readLT.py")
	from readLT import *
	
	###################################################
	#                     SETUP                       #
	###################################################
	
	# read parameters
	ifile = open('./params.in', 'r')
	print ifile
	for line in ifile :
		params = line.split()		# space delimited text
		if params[0] == 'METRIC' :
			STR = params[2]
#		if params[0] == 'REDVOL' :
#			redvol = params[2]
#		if params[0] == 'CONFIN' :
#			confin = params[2]
#		if params[0] == 'BENMOD' :
#			benmod = params[2]
#		if params[0] == 'CAPNUM' :
#			capnum = params[2]
	ifile.close()
	
	###################################################
	#               LUBRICATION THEORY                #
	###################################################
	
	# load data file
	(CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA, CURV) = readLT2(homedir, redvol, confin)

	REDVOL = []
	LAMBDA = []
	for i in range(len(AREA)) :
		a  = sqrt(AREA[i]/(4.0*pi))
		a3 = a*a*a
		REDVOL.append(VLME[i]/((4.0/3.0)*pi*a3))
		LAMBDA.append(a)
	REDVOL = np.column_stack(REDVOL)
	REDVOL = REDVOL[0]
	LAMBDA = np.column_stack(LAMBDA)
	LAMBDA = LAMBDA[0]
		
	if   (STR == 'DP') :
		xL = CA
		yL = DP
	elif (STR == 'EPS') :
		xL = CA
		yL = EPS
	elif (STR == 'AREA') :
		xL = CA
		yL = AREA
	elif (STR == 'VLME') :
		xL = CA
		yL = VLME
	elif (STR == 'REDVOL') :
		xL = CA
		yL = REDVOL
	elif (STR == 'GAMMAX') :
		xL = CA
		yL = GAMMAX
	elif (STR == 'TAUMAX') :
		xL = CA
		yL = TAUMAX
	
	# get reduced volume
	v     = float(redvol)/100.0

	return (STR, CA, EPS, REDVOL, LAMBDA, AREA, VLME, DP, ETA, CURV)

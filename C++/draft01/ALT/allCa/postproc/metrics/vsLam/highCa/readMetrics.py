def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readMetrics(redvol, maxrad, capnum) :
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import cm
	import numpy as np
	from scipy import interpolate, interpolate
	from math import log10, floor, pi, sqrt
	import os.path
	
	homedir = os.path.abspath(os.path.join("../../", os.pardir))

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
#		if params[0] == 'MAXRAD' :
#			maxrad = params[2]
#		if params[0] == 'BENMOD' :
#			benmod = params[2]
#		if params[0] == 'CAPNUM' :
#			capnum = params[2]
	ifile.close()
	
	###################################################
	#               LUBRICATION THEORY                #
	###################################################
	
	# load data file
#	(CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA) = readLT2(homedir, redvol, maxrad)
	(Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLT(homedir, redvol, maxrad, capnum)

	Dp = p[0] - p[-1]
	length = max(x) - min(x)
	eta = (pi*length/vlme)*(Dp/(8.0*(1.0 - eps)*length) - 1.0)
	maxgam = max(gam)

	return (STR, Ca, eps, area, vlme, Dp, eta, length, maxgam)
	#return (STR, CA, EPS, REDVOL, LAMBDA, AREA, VLME, DP, ETA)

def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readMetrics(redvol,capnum) :
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import cm
	import numpy as np
	from scipy import interpolate, integrate
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
	
	REDVOL = ['90','91','92','93','94','95','96','97','98','99','100']
	REDVOL = ['80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100']
	REDVOL = ['70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100']
	REDVOL = ['65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100']

	CONFIN = ['80','85','90','95']
	CONFIN = ['95','96','97','98','99']
	CONFIN = ['80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
	#CONFIN = ['70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99',]
#	CONFIN = ['80','82','84','86','88','90','92','94','96','98']
#	CONFIN = ['80','85','90','95','98']
#	CONFIN = ['90','91','92','93','94','95','96','97','98','99']
	#CONFIN = ['99']
	#CONFIN = ['99']

	DP = []
	EPS = []
	AREA = []
	VLME = []
	LAMBDA = []
	VOLRAT = []
	for confin in CONFIN :
		
		###################################################
		#               LUBRICATION THEORY                #
		###################################################
		
		# load data file
		#(CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA) = readLT2(homedir, redvol, confin)
		(Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLT(homedir, redvol, confin, capnum)

		# alternative way of calculating the pressure drop
		# - integrate the shear stress on the wall
		shear = []
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
			shear.append(-e*tx[i])
		shear = np.column_stack(shear)
		shear = shear[0]
		
		dp = 2.0*integrate.trapz(shear, s) 

		# reduced volume
		a = sqrt(area/(4.0*pi))
		a3 = a*a*a
		lam = a
		v = vlme/((4.0/3.0)*pi*a3)

		DP .append(dp)
		EPS.append(eps)
		AREA.append(area)
		VLME.append(vlme)
		LAMBDA.append(a)
		VOLRAT.append(v)
	
	DP     = np.column_stack(DP)
	EPS    = np.column_stack(EPS)
	AREA   = np.column_stack(AREA)
	VLME   = np.column_stack(VLME)
	LAMBDA = np.column_stack(LAMBDA)
	VOLRAT = np.column_stack(VOLRAT)

	DP     = DP[0]
	EPS    = EPS[0]
	AREA   = AREA[0]
	VLME   = VLME[0]
	LAMBDA = LAMBDA[0]
	VOLRAT = VOLRAT[0]
	
	return (STR, EPS, VOLRAT, LAMBDA, AREA, VLME, DP)

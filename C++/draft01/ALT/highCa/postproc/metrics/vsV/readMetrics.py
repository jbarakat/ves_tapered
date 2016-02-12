def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readMetrics(confin,capnum) :
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
	REDVOL = ['65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
	DP = []
	EPS = []
	AREA = []
	VLME = []
	LAMBDA = []
	for redvol in REDVOL :
		
		###################################################
		#               LUBRICATION THEORY                #
		###################################################
		
		# load data file
		#(CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA) = readLT2(homedir, redvol, confin)
		#(Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLT(homedir, redvol, confin, capnum)
		(eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLTHighCa(homedir, redvol, confin)
		
		v = float(redvol)/100.0
	
		# calculate lamc from the cubic equuation: 2 v x^3 - 3 x + 1 = 0
		lamc = 1.0
		if (abs(v - 1.00) < 1e-12) :
			lamc = 1.000000000000000000000000
	
		if (abs(v - 0.99) < 1e-12) :
			lamc = 1.090275107548953613318768
		if (abs(v - 0.98) < 1e-12) :
			lamc = 1.133537870150389282270677
		if (abs(v - 0.97) < 1e-12) :
			lamc = 1.169545928031930431602916
		if (abs(v - 0.96) < 1e-12) :
			lamc = 1.202032009320693326333173
		if (abs(v - 0.95) < 1e-12) :
			lamc = 1.232435708492447912665580
		if (abs(v - 0.94) < 1e-12) :
			lamc = 1.261494864992458404813790 
		if (abs(v - 0.93) < 1e-12) :
			lamc = 1.289649308437067534243365 
		if (abs(v - 0.92) < 1e-12) :
			lamc = 1.317187962872965134660298  
		if (abs(v - 0.91) < 1e-12) :
			lamc = 1.344314250171293006672972 
		if (abs(v - 0.90) < 1e-12) :
			lamc = 1.371179203625788627288043
		
		if (abs(v - 0.89) < 1e-12) :
			lamc = 1.397899861949529379050635
		if (abs(v - 0.88) < 1e-12) :
			lamc = 1.424570228180523274316834
		if (abs(v - 0.87) < 1e-12) :
			lamc = 1.451268169201854679431355
		if (abs(v - 0.86) < 1e-12) :
			lamc = 1.478059959100682705404164
		if (abs(v - 0.85) < 1e-12) :
			lamc = 1.505003385688966825138370
		if (abs(v - 0.84) < 1e-12) :
			lamc = 1.532149944429383541688818
		if (abs(v - 0.83) < 1e-12) :
			lamc = 1.559546432719025655752993
		if (abs(v - 0.82) < 1e-12) :
			lamc = 1.587236138747886422020625
		if (abs(v - 0.81) < 1e-12) :
			lamc = 1.615259749564796639821514
		if (abs(v - 0.80) < 1e-12) :
			lamc = 1.643656060707450592693672 

		if (abs(v - 0.79) < 1e-12) :
			lamc = 1.672462543253015249713266
		if (abs(v - 0.78) < 1e-12) :
			lamc = 1.701715807075084193076668
		if (abs(v - 0.77) < 1e-12) :
			lamc = 1.731451987830022394460647
		if (abs(v - 0.76) < 1e-12) :
			lamc = 1.761707077607608484461251
		if (abs(v - 0.75) < 1e-12) :
			lamc = 1.792517213974340291730173 
		if (abs(v - 0.74) < 1e-12) :
			lamc = 1.823918938509152593286319
		if (abs(v - 0.73) < 1e-12) :
			lamc = 1.855949433369632771788232
		if (abs(v - 0.72) < 1e-12) :
			lamc = 1.888646742600595040664463
		if (abs(v - 0.71) < 1e-12) :
			lamc = 1.922049983587091034764999
		if (abs(v - 0.70) < 1e-12) :
			lamc = 1.956199553113622973856379  
		
		if (abs(v - 0.69) < 1e-12) :
			lamc = 1.991137331820519094597204
		if (abs(v - 0.68) < 1e-12) :
			lamc = 2.026906890378410524999735
		if (abs(v - 0.67) < 1e-12) :
			lamc = 2.063553700384946879653673
		if (abs(v - 0.66) < 1e-12) :
			lamc = 2.101125352791323200747707
		if (abs(v - 0.65) < 1e-12) :
			lamc = 2.139671786567143914510814  
		if (abs(v - 0.64) < 1e-12) :
			lamc = 2.179245530295288117177149
		if (abs(v - 0.63) < 1e-12) :
			lamc = 2.219901959443900932366550
		if (abs(v - 0.62) < 1e-12) :
			lamc = 2.261699572184745456134028
		if (abs(v - 0.61) < 1e-12) :
			lamc = 2.304700286813593219417380
		if (abs(v - 0.60) < 1e-12) :
			lamc = 2.348969764079628399293969




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
				g = 8.0*(eps - 1.0) # limiting value as r --> 0
				e = 0.5*g

		Lcyl = 2.0*(lamc*lamc - 1.0)
		Lh   = 0.5*Lcyl
		for i in range(len(s)) :
			rr = r[i]
			xx = x[i]
			r2 = rr*rr
			if ((xx > -Lh or xx < Lh) and rr > 0.00001) :
				from math import log
				logr = log(rr)
			
				g =   (8.0/(1.0 - r2))
				g = g*(eps - 1.0 - (1.0 - r2)/(2.0*logr))
				g = g/(1.0 + r2 + (1.0 - r2)/logr)

				e = 0.25*g*(2.0 + (1.0 - r2)/logr) + 1.0/logr
			else : 
				e = 0.0
			#shear.append(-e*tx[i])
			shear.append(-e)
		shear = np.column_stack(shear)
		shear = shear[0]
		
		dp = 2.0*integrate.trapz(shear, x) 

		eps = 1.0 - r[len(r)/2]

		DP .append(dp)
		EPS.append(eps)
		AREA.append(area)
		VLME.append(vlme)
	
		for i in range(len(AREA)) :
			a  = sqrt(AREA[i]/(4.0*pi))
			a3 = a*a*a
			LAMBDA.append(a)
	DP     = np.column_stack(DP)
	EPS    = np.column_stack(EPS)
	AREA   = np.column_stack(AREA)
	VLME   = np.column_stack(VLME)
	LAMBDA = np.column_stack(LAMBDA)

	DP     = DP[0]
	EPS    = EPS[0]
	AREA   = AREA[0]
	VLME   = VLME[0]
	LAMBDA = LAMBDA[0]
			
	if   (STR == 'DP') :
		xL = REDVOL
		yL = DP
	elif (STR == 'EPS') :
		xL = REDVOL
		yL = EPS
	elif (STR == 'AREA') :
		xL = REDVOL
		yL = AREA
#	elif (STR == 'VLME') :
#		xL = REDVOL
#		yL = VLME
#	elif (STR == 'GAMMAX') :
#		xL = CA
#		yL = GAMMAX
#	elif (STR == 'TAUMAX') :
#		xL = CA
#		yL = TAUMAX
	
#	# get reduced volume
#	v     = float(redvol)/100.0

	return (STR, EPS, REDVOL, LAMBDA, AREA, VLME, DP)

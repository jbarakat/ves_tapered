def load_src(name, fpath):
	import os, imp 
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

# read file (lubrication theory)
def readFile(path) :
	import numpy as np
	from numpy import sin, cos
	from scipy import integrate
	import math
	from math import pi

	data = np.loadtxt(path, unpack=True,skiprows=1)
	#print 'Reading ' + path + ' ...'
		
	# initialize arrays
	s	  = []
	r	  = []
	psi = []
	cs	= []
	qs	= []
	p	  = []
	gam = []
	A	  = []
	V	  = []
	
	# grab data
	Q2  = data[ 9,0]
	S	  = data[10,0]
	s	  .append(data[0,:]*S)
	r	  .append(data[1,:])
	psi .append(data[2,:])
	cs	.append(data[3,:])
	qs	.append(data[4,:])
	p	  .append(data[5,:])
	gam .append(data[6,:])
	A	  .append(data[7,:])
	V	  .append(data[8,:])
	
	s	  = np.column_stack(s	)
	r	  = np.column_stack(r	)
	psi = np.column_stack(psi)
	cs	= np.column_stack(cs )
	qs	= np.column_stack(qs )
	p	  = np.column_stack(p	)
	gam = np.column_stack(gam)
	A	  = np.column_stack(A	)
	V	  = np.column_stack(V	)
	
	s   = s  [:,0]
	r	  = r  [:,0]
	psi = psi[:,0]
	cs	= cs [:,0]
	qs	= qs [:,0]
	p	  = p  [:,0]
	gam = gam[:,0]
	A	  = A  [:,0]
	V	  = V  [:,0]

	area = A[-1]
	vlme = V[-1]
	
	cospsi = cos(psi)
	sinpsi = sin(psi)
		
	# get axial coordinate
	x	 = integrate.cumtrapz(cospsi, s)
	x	 = np.concatenate(([0.0], x))
	x	 = x - 0.5*(max(x) + min(x))

	# get film parameter
	eps = Q2

	# get normal and tangent
	nx = []
	nr = []
	tx = []
	tr = []
	for i in range(len(x)) :
		nx.append( sinpsi[i])
		nr.append( cospsi[i])
		tx.append( cospsi[i])
		tr.append(-sinpsi[i])
	
	# get azimuthal curvature
	cphi = []
	for i in range(len(x)) :
		if (r[i] > 0.02) :
			ccphi = -cospsi[i]/r[i]
		else :
			ccphi = cs[i]
		cphi.append(ccphi)
	cphi = np.column_stack(cphi)
	cphi = cphi[0]

	# get shear stress
	Qpi = -eps

	tau = []
	for i in range(len(r)) :
	#	if (r[i] > 232) :
	#		from math import log
	#		r2 = r[i]*r[i]
	#		logr = log(r[i])
	#		rlogr = r[i]*logr
	#		irlogr = 1.0/rlogr

	#		g = -(16.0/(1.0 - r2))
	#		g = g*(0.5 + (1.0 - r2)/(4.0*logr) + 0.5*Qpi)
	#		g = g/(1.0 + r2 + (1.0 - r2)/logr)
	#		
	#		e = 0.25*g*((1 - r2)*irlogr + 2.0*r[i]) + irlogr
	#	else :
	#	g  = -8.0*(1.0 + Qpi); # limiting value as r --> 0
	
		Dr   = 1.0 - r[i]
		Dr2  = Dr  *Dr     
		Dr3  = Dr  *Dr2 
		Dr4  = Dr2 *Dr2 
		Dr5  = Dr2 *Dr3 
		Dr6  = Dr3 *Dr3 
		Dr7  = Dr3 *Dr4 
		Dr8  = Dr4 *Dr4 
		Dr9  = Dr4 *Dr5 
		Dr10 = Dr5 *Dr5 
		Dr11 = Dr5 *Dr6 
		Dr12 = Dr6 *Dr6 
		Dr13 = Dr6 *Dr7 
		Dr14 = Dr7 *Dr7 
		Dr15 = Dr7 *Dr8 
		Dr16 = Dr8 *Dr8 
		Dr17 = Dr8 *Dr9 
		Dr18 = Dr9 *Dr9 
		Dr19 = Dr9 *Dr10
		Dr20 = Dr10*Dr10
	
		e = (3.0/Dr2)*Qpi                  
		e = e + (2.0/Dr )*(1.0 + Qpi)          
		e = e + 1.0 + 1.45*Qpi                 
		e = e + Dr  *(0.7 + 1.15*Qpi)          
		e = e + Dr2 *(0.55 + 0.981786*Qpi)     
		e = e + Dr3 *(0.469286 + 0.882857*Qpi) 
		e = e + Dr4 *(0.423214 + 0.820839*Qpi) 
	#	e = e + Dr5 *(0.395036 + 0.778982*Qpi) 
	#	e = e + Dr6 *(0.376375 + 0.748554*Qpi) 
	#	e = e + Dr7 *(0.362969 + 0.724934*Qpi) 
	#	e = e + Dr8 *(0.352610 + 0.705622*Qpi) 
	#	e = e + Dr9 *(0.344131 + 0.689218*Qpi) 
	#	e = e + Dr10*(0.336895 + 0.674906*Qpi) 
	#	e = e + Dr11*(0.330541 + 0.662184*Qpi) 
	#	e = e + Dr12*(0.324853 + 0.650726*Qpi) 
	#	e = e + Dr13*(0.319695 + 0.640305*Qpi) 
	#	e = e + Dr14*(0.314973 + 0.630758*Qpi) 
	#	e = e + Dr15*(0.310620 + 0.621958*Qpi) 
	#	e = e + Dr16*(0.306587 + 0.613808*Qpi) 
	#	e = e + Dr17*(0.302833 + 0.606228*Qpi) 
	#	e = e + Dr18*(0.299326 + 0.599152*Qpi) 
	#	e = e + Dr19*(0.296039 + 0.592525*Qpi) 
	#	e = e + Dr20*(0.292949 + 0.586299*Qpi) 

		tau.append(e)
	
#	from matplotlib import pyplot as plt
#	plt.plot(x, cs)
#	plt.show()

	Q = 0.5*Q2

	return (eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S)




# get profiles at fixed v, conf, Ca
def readLT(homedir, redvol, confin, capnum) :	
	import numpy as np
	from numpy import nan
	import glob
	from math import pi
	
	# load data files
	outdir	= homedir + '/../output/v' + redvol + '/conf' + confin + '/*Ca' + capnum + '.dat'
	substring = homedir + '/../output/v' + redvol + '/conf' + confin + '/sln_v' + redvol + '_conf' + confin + '_Ca'
	files	= sorted(glob.glob(outdir))
	nfiles = len(files)
	lss = len(substring)

	Ca   =  nan 
	eps  =  nan 
	area =  nan 
	vlme =  nan 
	s    = [nan]
	x    = [nan]
	r    = [nan]
	nx   = [nan]
	nr   = [nan]
	tx   = [nan]
	tr   = [nan]
	cs   = [nan]
	cphi = [nan]
	psi  = [nan]
	qs   = [nan]
	p    = [nan]
	tau  = [nan]
	gam  = [nan]

	# import data
	if nfiles > 0 :   #### NEED TO REVISE THIS .. WHAT TO DO WHEN THERE IS NO FILE...
		f = files[0]
		#print 'Reading ' + f + ' ...'
		(eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S) = readFile(f)
	
		# get capillary number
		lf = len(f)
		Ca = f[lss:lf-4]
		if (Ca[0] == '0') :
			Ca = Ca[0] + '.' + Ca[1:]
		Ca = float(Ca)
	
		# scale pressure and surface tensions
		p   = p  /Ca
		gam = gam/Ca
		qs  = qs /Ca
	
	# return data
	return (Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S)
		




# get mobility parameters over range of Ca at fixed v, conf
def readLT2(homedir, redvol, confin) :	
	import numpy as np
	import glob
	#from math import pi, sin, cos, log
	from scipy import integrate
	
	# load data files
	outdir	= homedir + '/../output/v' + redvol + '/conf' + confin + '/*.dat'
	substring = homedir + '/../output/v' + redvol + '/conf' + confin + '/sln_v' + redvol + '_conf' + confin + '_Ca'
	files	= sorted(glob.glob(outdir))
	nfiles = len(files)
	lss = len(substring)

	# import data
	CA     = []	# capillary number
	EPS    = []	# film or leakback parameter = 1 - V/U
	LENGTH = [] # axial length
	AREA   = [] # volume
	VLME   = [] # volume
	FX     = [] # axial force per unit circumference on the particle
	FR     = [] # radial force per unit circumference on the particle
	DP     = [] # pressure difference across the particle
	GAMMAX = [] # maximum tension
	TAUMAX = [] # maximum shear stress
	ETA    = [] # effective viscosity
	CURV   = [] # rear mean curvature
	for f in files :
		# read file
		(eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S) = readFile(f)
		
		EPS   .append(eps)
		AREA  .append(area)
		VLME  .append(vlme)
	
		# get capillary number
		lf = len(f)
		Ca = f[lss:lf-4]
		if (Ca[0] == '0') :
			Ca = Ca[0] + '.' + Ca[1:]
		Ca = float(Ca)
		CA .append(Ca)
	
		# scale pressure and surface tension
		p   = p  /Ca
		gam = gam/Ca

		# get axial length
		L = max(x) - min(x)
		LENGTH.append(L)
		
		# get force components (per unit circumference) on the particle
		fnr = []
		ftr = []
		for i in range(len(s)) :
			fnr.append( -p[i]*r[i])
			ftr.append(tau[i]*r[i])
		fnr = np.column_stack(fnr)
		ftr = np.column_stack(ftr)
		fnr = fnr[0]
		ftr = ftr[0]
		
		ffx = []
		ffr = []
		for i in range(len(s)) :
			ffx.append(fnr[i]*nx[i] + ftr*tx[i])
			ffr.append(fnr[i]*nr[i] + ftr*tr[i])
		ffx = np.column_stack(ffx)
		ffr = np.column_stack(ffr)
		ffx = ffx[0]
		ffr = ffr[0]
		fx = integrate.trapz(ffx, s)
		fr = integrate.trapz(ffr, s)

		FX.append(fx)
		FR.append(fr)
		
		# get pressure drop
		p1 = p[-1]
		for i in range(len(s)) :
			p0 = p[i]
			if (cs[i] <= 0) :
				break
		dp = p0 - p1

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
			#shear.append(-e*tx[i])
			shear.append(-e)
		shear = np.column_stack(shear)
		shear = shear[0]
		
		dp = 2.0*integrate.trapz(shear, x) 

		DP .append(dp)

		# get maximum tension
		gammax = max(gam)
		GAMMAX.append(gammax)

		# get maximum shear stress
		taumax = max(tau)
		TAUMAX.append(taumax)
	
		# get effective viscosity
		V   = 1 - eps
		#eta = dp/(8*Ca*V*L)
		eta = dp/(8*V*L)
		ETA .append(eta)

		# get value of rear curvature
		curv = cs[0]
		CURV.append(curv)


	DATA = np.column_stack([CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA, CURV])
	DATA = sorted(DATA, key=lambda i:i[0])
	DATA = np.column_stack(DATA)

	CA     = DATA[0]
	EPS    = DATA[1]
	LENGTH = DATA[2]
	AREA   = DATA[3]
	VLME   = DATA[4]
	FX     = DATA[5]
	FR     = DATA[6]
	DP     = DATA[7]
	GAMMAX = DATA[8]
	TAUMAX = DATA[9]
	ETA    = DATA[10]
	CURV   = DATA[10]

	return (CA, EPS, LENGTH, AREA, VLME, FX, FR, DP, GAMMAX, TAUMAX, ETA, CURV)



# get profiles at fixed v, conf, Ca
def readInput(homedir, redvol, confin, capnum) :
	import numpy as np
	from numpy import nan
	import glob
	from scipy import integrate, interpolate
	from math import log10, floor, pi
	import math
	import cmath
	
	# load data files
	outdir	= homedir + '/sln_*Ca' + capnum + '.dat'
	substring = homedir + '/sln_v' + redvol + '_conf' + confin + '_Ca'
	files	= sorted(glob.glob(outdir))
	nfiles = len(files)
	lss = len(substring)

	Ca   =  nan 
	eps  =  nan 
	area =  nan 
	vlme =  nan 
	s    = [nan]
	x    = [nan]
	r    = [nan]
	nx   = [nan]
	nr   = [nan]
	tx   = [nan]
	tr   = [nan]
	cs   = [nan]
	cphi = [nan]
	qs   = [nan]
	p    = [nan]
	tau  = [nan]
	gam  = [nan]

	# import data
	if nfiles > 0 :   #### NEED TO REVISE THIS .. WHAT TO DO WHEN THERE IS NO FILE...
		f = files[0]
		#print 'Reading ' + f + ' ...'
		(eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S) = readFile(f)
	
		# get capillary number
		lf = len(f)
		Ca = f[lss:lf-4]
		if (Ca[0] == '0') :
			Ca = Ca[0] + '.' + Ca[1:]
		Ca = float(Ca)
	
		# scale pressure and surface tensions
		p   = p  /Ca
		gam = gam/Ca
		qs  = qs /Ca

	################################################################ now modify for tapered problem...
	
	t = np.linspace(0.0,1.0,len(s));
	for i in range(len(s)) :
		t[i] = s[i]/s[-1]

	# calculate center of mass
	kernel = np.linspace(0.0,1.0,len(s))
	for i in range(len(s)) :
		kernel[i] = x[i]*r[i]*r[i]*math.cos(psi[i])
		
	xcom = (pi/vlme)*integrate.trapz(kernel, s)
	
	# center vesicle at center of mass
	for i in range(len(x)) :
		x[i] = x[i] - xcom
	
	# get the nominal radius of the vesicle
	Rves  = math.sqrt(area / (4.0*pi))
	Rv3   = Rves*Rves*Rves

	# define the original characteristic scales
	# (vesicle speed, tube radius, and viscosity)
	Rtube = 1.0
	U     = 1.0

	# set the center-of-mass to be the origin
	XCOM  = 0.0
	
	# get reduced volume and bending modulus
	v     = vlme/(4.0*pi*Rv3/3.0)
	kb    = 1.0/Ca
	
	# get critical tube radius
	cf = (cmath.sqrt(3)/2.0)*1j
	a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
	a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)
	Rc = -0.5*(a1 + a2) + cf*(a1 - a2)
	Rc = Rc.real
	if abs(v - 1.0) < 1e-12 :
		Rc = 1.0
	
	# get confinement
	conf  = Rc*Rves

	# rescale all variables wrt the pressure drop, 
	# the suspending fluid viscosity (mu = 1), 
	# and the nominal radius of the vesicle
	dp    = p[0] - p[-1]
	kb    = kb/(dp*Rv3)
	Rtube = Rtube / Rves
	U     = U /(dp*Rves)
	Q     = Q /(dp*Rves*Rves)
	Uflow = U - 2*Q/Rtube
	S     = S / Rves
	area = area/(Rves*Rves)
	vlme = vlme/Rv3
	for i in range(len(s)) :
		p  [i] = p  [i] / dp
		tau[i] = tau[i] / dp
		gam[i] = gam[i] / (dp*Rves)
		qs [i] = qs [i] / (dp*Rves)

		r   [i] = r   [i] / Rves
		x   [i] = x   [i] / Rves
		s   [i] = s   [i] / Rves
		cs  [i] = cs  [i] * Rves
		cphi[i] = cphi[i] * Rves
		A   [i] = A   [i] / (Rves*Rves)
		V   [i] = V   [i] / Rv3

	return (v, conf, kb, area, vlme, t, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S, Rtube, U, XCOM)

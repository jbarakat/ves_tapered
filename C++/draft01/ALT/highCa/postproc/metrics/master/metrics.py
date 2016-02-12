from readMetrics import readMetrics
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import splrep, splev, interp1d
from math import log10, floor

MAXRAD = ['90', '92', '94', '96']
MAXRAD = ['90', '99']
#MAXRAD = ['80','82','84','86','88','90', '92', '94', '96', '98']
MAXRAD = ['90','91','92','93','94','95','96','97','98','99']
#MAXRAD = ['99']

color = iter(cm.brg(np.linspace(0.0,0.8, len(MAXRAD))))
fs = 24
plt.figure(figsize=(14,12))
for maxrad in MAXRAD :
	c = next(color)
	(STR, CA, EPS, REDVOL, LAMBDA, AREA, VLME, LENGTH, DP, ETA) = readMetrics(maxrad)

	v = round(REDVOL[0],2)
	lam = LAMBDA[0]

	# calculate lamc from the cubic equuation: 2 v x^3 - 3 x + 1 = 0
	if (abs(v - 0.95) < 1e-12) :
		lamc = 1.2324
	if (abs(v - 0.90) < 1e-12) :
		lamc = 1.3712
	if (abs(v - 0.85) < 1e-12) :
		lamc = 1.5050
	if (abs(v - 0.80) < 1e-12) :
		lamc = 1.6437
	if (abs(v - 0.75) < 1e-12) :
		lamc = 1.7925
	if (abs(v - 0.70) < 1e-12) :
		lamc = 1.9562
	if (abs(v - 0.65) <= 0.01) :
		lamc = 2.1397
	
	lam = lam/lamc
	delta = 1 - lam

#	fc = 7.5112
#	n  = 0.8776
#	n  = 0.86
#	#fc = 19.041
#	#n  = 0.9931
#	delta = pow(delta,n)
#	delta = 1
#	fc = 1
	
	lam = round(lam,2)
	v = round(v,2)

#	CE = CA
#	DPE = DP
#	for i in range(len(CE)) :
#		CE[i] = CA[i]/(EPS[i]*EPS[i])
#		DPE[i] = DP[i]*EPS[i]
#	xL = CE
#	yL = DPE

	xL = CA/(delta*delta)
	if (STR == 'DP') :
		yL = DP*delta
	#	yL = DP
	#	for i in range(len(yL)) :
	#		yL[i] = yL[i]/(8.0*(1.0 - EPS[i])*LENGTH[i])
	if (STR == 'EPS') :
		yL = (EPS/delta)
	
	
	## interpolate to xi using splines
	#tck = interpolate.splrep(x, y, s=0)
	#yi = interpolate.splev(xi, tck, der=0)

	xi  = np.logspace(floor(log10(min(xL))), floor(log10(max(xL))), 10000)
	tck = splrep(xL, yL, s=0)
	yi  = splev (xi, tck, der=0)

#	f = interp1d(xL, yL)
#	yi = f(xi)

	plt.plot(xi,yi,'-',c=c)
	plt.plot(xL,yL,'.',c=c)
	
	if (maxrad == MAXRAD[0] or maxrad == MAXRAD[-1]) :
		plt.text(1000, yi[-1], str(lam),fontsize = fs-8)
	
	if (maxrad == MAXRAD[0]) :
		lam0 = lam
	
	if (maxrad == MAXRAD[-1]) :
		lam1 = lam

# set plot preferences
plt.title('$v$ = 0.' + str(int(v*100)) + ', $\lambda/\lambda_c$ = ' + str(lam0) + ' to ' + str(lam1),fontsize=fs)
plt.xscale('log')
#plt.yscale('symlog')
plt.xticks(fontsize=fs-4)
plt.yticks(fontsize=fs-4)
plt.xlim((0.01,1000))

plt.xlabel('$\mathrm{Ca} = \mu U R^2 / (\epsilon^2 \kappa)$',fontsize=fs)
if (STR == 'DP') :
	plt.ylabel('$\epsilon \Delta p R/ (\mu U)$',fontsize=fs)
if (STR == 'EPS') :
	plt.ylabel('$Q/(\pi \epsilon R^2 U) = 1 - V/U$',fontsize=fs)
if (STR == 'AREA') :
	plt.ylabel('$A/R^2$',fontsize=fs)
if (STR == 'VLME') :
	plt.ylabel('$\Omega/R^3$',fontsize=fs)
if (STR == 'REDVOL') :
	plt.ylabel('$v$',fontsize=fs)
if (STR == 'GAMMAX') :
	plt.ylabel('$\gamma_{\mathrm{max}}/ (\mu U)$',fontsize=fs)
if (STR == 'TAUMAX') :
	plt.ylabel(r'$\tau_{\mathrm{max}}R/ (\mu U)$',fontsize=fs)


plt.show()

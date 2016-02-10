from readMetrics import readMetrics
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import splrep, splev
from scipy.stats import linregress
from math import log10, floor, sqrt, pi

MAXRAD = ['90', '92', '94']
#MAXRAD = ['80','82','84','86','88','90', '92', '94', '96', '98']
MAXRAD = ['90','91','92','93','94','95','96','97','98','99']
#MAXRAD = ['90', '92', '94', '96']

REDVOL = ['95', '90', '85','80','75','70','65']
#REDVOL = ['95', '90', '85','80','75','70']
#REDVOL = ['80','75','70']

capnum = '0000010'

color = iter(cm.brg(np.linspace(0.0,0.8, len(REDVOL))))
fs = 24
plt.figure(figsize=(14,12))

maxL = 0
nspline = 1000

vecRedVol = []
vecExponent = []
vecCoeff = []
vecRSquared = []

for redvol in REDVOL :
	c = next(color)
	x = []
	y = []
	L = []
	for maxrad in MAXRAD :
		(STR, Ca, eps, area, vlme, Dp, l) = readMetrics(redvol, maxrad, capnum)

		lam = sqrt(area/(4.0*pi))
		v = round(float(redvol)/100,2)

		# calculate lamc from the cubic equuation: 2 v x^3 - 3 x + 1 = 0
		if (abs(v - 0.95) <= 0.01) :
			lamc = 1.2324
		if (abs(v - 0.90) <= 0.01) :
			lamc = 1.3712
		if (abs(v - 0.85) <= 0.01) :
			lamc = 1.5050
		if (abs(v - 0.80) <= 0.01) :
			lamc = 1.6437
		if (abs(v - 0.75) <= 0.01) :
			lamc = 1.7925
		if (abs(v - 0.70) <= 0.01) :
			lamc = 1.9562
		if (abs(v - 0.65) <= 0.01) :
			lamc = 2.1397
		
		lam = lam/lamc
		delta = 1 - lam
	
		x.append(delta)
		if (STR == 'DP') :
			y.append(Dp )
		elif (STR == 'EPS') :
			y.append(eps)
		L.append(l)
	
	x = np.column_stack(x)
	y = np.column_stack(y)
	L = np.column_stack(L)
	x = x[0]
	y = y[0]
	L = L[0]
	
	DATA = np.column_stack([x,y,L])
	DATA = sorted(DATA, key=lambda i:i[0])
	DATA = np.column_stack(DATA)

	x = DATA[0]
	y = DATA[1]
	L = DATA[2]

	# interpolate to xi using splines
	xi  = np.linspace(min(x),max(x),nspline)
	tck = splrep(x, y, s=0)
	yi  = splev (xi, tck, der=0)
	
	# fit linear curve to log-log data
	logx  = []
	logy  = []
	for i in range(len(x)) :
	  logx.append(log10(x[i]))
	  logy.append(log10(y[i]))
	logxp = logx
	
	fit   = np.polyfit(logx, logy, 1)
	poly  = np.poly1d (fit)
	dpoly = np.polyder(poly)
	logyp = poly(logxp)

	slope, intercept, r_value, p_value, std_err = linregress(logxp,logyp)

	coeff = pow(10,intercept)

	slope = round(slope,4)
	coeff = round(coeff,4)
	rsqr  = round(r_value**2,4)

	xi = np.linspace(min(x),max(x),nspline)
	yi = []
	for i in range(len(xi)) :
		yi.append(coeff*pow(xi[i],slope))

	vecRedVol  .append(v)
	vecExponent.append(slope)
	vecCoeff   .append(coeff)
	vecRSquared.append(rsqr)

#	plt.plot(xi,yi,'-',c=c, label = 'vesicle ($v$ = 0.' + redvol + '), $n$ = ' + str(slope))
	plt.plot(xi,yi,'-',c=c, label = '$v$ = 0.' + redvol)
	plt.plot(x,y,'.',c=c)

if (STR == 'DP') :
	# plot result for bullets
	l = round(maxL,1)
	l = 10
	fc  = 2.0
	yi = []
	x = np.linspace(0.06, 0.22, nspline)
	for m in range(len(x)) :
	  eps = pow(x[m],-1.0)
	  yi.append(fc*l*eps)
	
	#plt.plot(x, yi, 'k--', label = 'rigid bullet ($L$ = ' + str(l) + '$R$), $n$ = -1')
	plt.plot(x, yi, 'k--', label = 'bullet ($L$ = ' + str(l) + '$R$)')
	
	# plot results for spheres
	x = np.linspace(0.02, 0.2, nspline)
	fc  = 4.0*sqrt(2.0)*pi
	yi = []
	for m in range(len(x)) :
	  eps = pow(x[m],-0.5)
	  yi.append(fc*eps)
	
	#plt.plot(x, yi, 'k-.', label = 'rigid sphere, $n$ = -0.5')
	plt.plot(x, yi, 'k-.', label = 'sphere')

elif (STR == 'EPS') :
	# plot results for bullets
	x = np.linspace(0.02, 0.18, nspline)
	fc = 1.0
	yi = []
	for m in range(len(x)) :
		yi.append(fc*x[m])
	
	#plt.plot(x, yi, 'k--', label = 'rigid bullet ($L$ = ' + str(l) + '$R$), $n$ = -1')
	plt.plot(x, yi, 'k--', label = 'bullet')

	# plot results for spheres
	x = np.linspace(0.02, 0.16, nspline)
	fc = 4.0/3.0
	yi = []
	for m in range(len(x)) :
		yi.append(fc*x[m])
	
	#plt.plot(x, yi, 'k-.', label = 'rigid sphere, $n$ = -0.5')
	plt.plot(x, yi, 'k-.', label = 'sphere')


print vecRedVol
print vecExponent
print vecCoeff
print vecRSquared

vmin = REDVOL[-1]
vmax = REDVOL[0]
	
# set plot preferences
plt.title('$v$ = 0.' + vmin + ' to 0.' + vmax + ', $\mathrm{Ca}$ = ' + capnum,fontsize=fs)
plt.xticks(fontsize=fs-4)
plt.yticks(fontsize=fs-4)

plt.xlabel('$\epsilon = 1 - \lambda/\lambda_c$',fontsize=fs)

if (STR == 'DP') :
	plt.ylabel('$\Delta p R/ (\mu U)$',fontsize=fs)
elif (STR == 'EPS') :
	plt.ylabel('$Q/(\pi R^2 U) = 1 - V/U$',fontsize=fs)
#	plt.xlim((0,0.25))
#	plt.ylim((0,0.2))
	
leg = plt.legend(fontsize=fs-4,loc='best')
leg.get_frame().set_linewidth(0.0)

plt.xscale('log')
plt.yscale('log')
plt.show()

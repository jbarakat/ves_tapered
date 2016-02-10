from readMetrics import readMetrics
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import splrep, splev, interp1d
from math import log10, floor

CONFIN = ['80','85','90','95']
CONFIN = ['95','96','97','98','99']
CONFIN = ['80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98']
CONFIN = ['80','82','84','86','88','90','92','94','96','98']
CONFIN = ['80','85','90','95','98']
CONFIN = ['90','91','92','93','94','95','96','97','98','99']
#CONFIN = ['99']
#CONFIN = ['99']

CAPNUM = ['1']

#fs = 24
#plt.figure(figsize=(14,12))
fs = 10
plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
#for redvol in REDVOL :
for capnum in CAPNUM :
	color = iter(cm.brg(np.linspace(0.0,0.8, len(CONFIN))))
	for confin in CONFIN :
		c = next(color)
		#(STR, CA, EPS, REDVOL, LAMBDA, AREA, VLME, DP, ETA) = readMetrics(confin)
		(STR, EPS, REDVOL, LAMBDA, AREA, VLME, DP) = readMetrics(confin,capnum)
	
	#	CE = CA
	#	DPE = DP
	#	for i in range(len(CE)) :
	#		CE[i] = CA[i]/(EPS[i]*EPS[i])
	#		#CE[i] = CA[i]/(delta*delta)
	#		DPE[i] = DP[i]*EPS[i]
	#		#DPE[i] = DP[i]*delta
	#	CA = CE
	#	DP = DPE
		REDVOL2 = REDVOL
		for i in range(len(REDVOL)) :
			v = REDVOL[i]
			v = float(v)/100
			REDVOL2[i] = v
		REDVOL = REDVOL2
	
		xL = REDVOL
		if (STR == 'DP') :
			yL = DP
		if (STR == 'EPS') :
			yL = EPS
	
		## interpolate to xi using splines
		#tck = interpolate.splrep(x, y, s=0)
		#yi = interpolate.splev(xi, tck, der=0)
	
#		xi  = np.logspace(floor(log10(min(xL))), floor(log10(max(xL))), 10000)
#		xi  = np.logspace(-4, floor(log10(max(xL))), 10000)
#		tck = splrep(xL, yL, s=0)
#		yi  = splev (xi, tck, der=0)
	#	xi = xL
	#	yi = yL
	
	#	f = interp1d(xL, yL)
	#	yi = f(xi)
	
		#plt.plot(xi,yi,'-',c=c, label = '$\lambda/\lambda_c$ = ' + str(lam))
#		plt.plot(xi,yi,'-',c=c)
		plt.plot(xL,yL,'.',c=c)
			
		if (confin == CONFIN[0] or confin == CONFIN[-1]) :
			plt.text(xL[-1], yL[-1], '0.' + str(confin),fontsize = fs)
		
#		if (confin == CONFIN[0]) :
#			lam0 = lam
#		
#		if (confin == CONFIN[-1]) :
#			lam1 = lam

eps = EPS[0]
vvec = []
press = []
for i in range(65, 101, 1) :
	v = float(i)/100.0

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
	
	eps = 0.01
	dp = 4.0*(lamc*lamc - 1.0)/eps

	# for some reason, this is only accurate to a shift of 170...
#	dp = dp + 172.5
	dp = dp + 178
#	dp = dp + 63.2/pow(eps,0.5)

	vvec.append(v)
	press.append(dp)

vvec = np.column_stack(vvec)
press = np.column_stack(press)
vvec = vvec[0]
press = press[0]

print vvec
print press

plt.plot(vvec,press,'k--')



	
	#	v = round(REDVOL[0],2)
	#	lam = LAMBDA[0]
	#
	#	
	#	
	#	lam = lam/lamc
	#	delta = 1.0 - lam
	#	
	#	lam = round(lam,2)
	#	v = round(v,2)

# set plot preferences
#plt.title('$v$ = 0.' + str(int(v*100)) + ', $\lambda/\lambda_c$ = ' + str(lam0) + ' to ' + str(lam1),fontsize=fs)
#plt.xscale('log')
#plt.yscale('log')
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
#plt.xlim((0.000001,1000))
#plt.ylim((0,20))

plt.xlabel('$v$',fontsize=fs)
if (STR == 'DP') :
	plt.ylabel('$\Delta p R/ (\mu U)$',fontsize=fs)
if (STR == 'EPS') :
	#plt.ylabel('$Q/(\pi R^2 U) = 1 - V/U$',fontsize=fs)
	plt.ylabel('$\epsilon = 2Q/(RU)$',fontsize=fs)
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
	
#leg = plt.legend(fontsize=fs-4,loc='best')
#leg.get_frame().set_linewidth(0.0)


plt.show()

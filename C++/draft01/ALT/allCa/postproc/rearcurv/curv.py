from read import readMetrics
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from scipy.interpolate import splrep, splev, interp1d
from math import log10, floor

#CONFIN = ['80','85','90','95']
#CONFIN = ['95','96','97','98','99']
CONFIN = ['70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
#CONFIN = ['70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
#CONFIN = ['80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
#CONFIN = ['90','91','92','93','94','95','96','97','98','99']
#CONFIN = ['90']
#CONFIN = ['80','82','84','86','88','90','92','94','96','98']
#CONFIN = ['80','85','90','95','98']
#CONFIN = ['98']
#CONFIN = ['80','99']
#CONFIN = ['99']

REDVOL = ['65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80','81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99','100']
#REDVOL = ['100']
#REDVOL = ['99']
#REDVOL = ['95']
#REDVOL = ['90']
#REDVOL = ['85']
#REDVOL = ['80']
#REDVOL = ['75']
#REDVOL = ['70']
#REDVOL = ['65']

CAPNUM = ['50']
CAPNUM = ['1']
CAPNUM = ['2']
CAPNUM = ['5']
CAPNUM = ['10']
CAPNUM = ['100']
CAPNUM = ['500']
#CAPNUM = ['0010']

fs = 24
plt.figure(figsize=(12,6))
#fs = 10
#plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

lspan = np.arange(0.7,0.851,0.01)
lspan = np.arange(0.7,0.991,0.01)
#lspan = np.arange(0.8,0.991,0.01)

vspan = np.arange(0.65,0.951,0.01)
vspan = np.arange(0.65,1.001,0.01)

x, y = np.meshgrid(vspan,lspan)
xx, yy = np.meshgrid(vspan,lspan)
z = xx
#print y
#print x[0][0]
#print y[0]

capnum = CAPNUM[0]
for i in range(len(REDVOL)) :
	redvol = REDVOL[i]
	for j in range(len(CONFIN)) :
		confin = CONFIN[j]
		(v, lam, csrear) = readMetrics(redvol, confin, capnum)
		
		#print csrear
		#z[j][i] = csrear
		if csrear > 0 :
			z[j][i] = 1.0
		else :
			z[j][i] = 0.0
#print z
#print x

#plt.contourf(x,y,z,1000)
plt.contourf(x,y,z)
cbar = plt.colorbar()

## set plot preferences
##plt.title('$v$ = 0.' + str(int(v*100)) + ', $\lambda/\lambda_c$ = ' + str(lam0) + ' to ' + str(lam1),fontsize=fs)
#plt.xscale('log')
#plt.yscale('log')
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
#plt.xlim((0.000001,1000))
##plt.ylim((0,20))
#
plt.xlabel('$v$',fontsize=fs)
plt.ylabel('$\lambda/\lambda_c$',fontsize=fs)
##leg = plt.legend(fontsize=fs-4,loc='best')
##leg.get_frame().set_linewidth(0.0)
#
#
plt.show()

from readProfiles import readProfiles
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from pylab import *
import numpy as np
from scipy import interpolate, interpolate
from math import log10, floor, pi

(Ca, eps, v, lam, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readProfiles()

# reflect shape across symmetry axis
xx  = np.concatenate((x, np.flipud( x)))
rr  = np.concatenate((r, np.flipud(-r)))
		
v = round(v,2)

# calculate lamc from the cubic equuation: 2 v x^3 - 3 x^2 + 1 = 0
lamc = 1
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
if (abs(v - 0.61) <= 0.01) :
	lamc = 2.3047

lam = lam/lamc
lam = round(lam,2)
delta = 1 - lam

if (Ca >= 1) :
	Ca = int(Ca)

# wall
xw = [-6,6]
uw = [1,1]
lw = [-1,-1]
#wall = []
#for i in range(len(x)) :
#	wall.append(1.0)
#wall = np.column_stack(wall)
#wall = wall[0]

#fig = plt.subplots(figsize=(12,12))
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True, figsize=(12,12))
ax1 = plt.subplot(2,1,1)
ax2 = ax1.twinx()
ax3 = plt.subplot(2,1,2, sharex=ax1)
fig.subplots_adjust(hspace=0)

ax1.plot(x,  p, 'b-')
ax2.plot(x,gam, 'g-')
ax3.plot(xx,rr, 'k-')
#ax3.plot(xw,uw, 'k--')
##ax1.plot(xw,lw, 'k--')

fs = 24
ax1.set_ylabel('$p R/(\mu U)$', color='b',fontsize=fs)
ax2.set_ylabel('$\gamma/(\mu U)$', color='g',fontsize=fs)
ax3.set_ylabel('$r/R$',fontsize=fs)
ax3.set_xlabel('$x/R$',fontsize=fs)

setp(ax1.get_xticklabels(), visible=False)
setp(ax1.get_yticklabels(), color='b',fontsize=fs-4)
setp(ax2.get_yticklabels(), color='g',fontsize=fs-4)
setp(ax3.get_yticklabels(), fontsize=fs-4)
setp(ax3.get_xticklabels(), fontsize=fs-4)
ax3.set_xlim((min(x)-0.5,max(x)+0.5))
ax3.set_yticklabels([1,0.5,0])
ax1.set_title('$v$ = 0.' + str(int(v*100)) + ', $\lambda/\lambda_c$ = ' + str(lam) + ', $\mathrm{Ca}$ = ' + str(Ca),fontsize=fs)

#ax1.set_ylim((-120,20))
#ax2.set_ylim((-5,40))

#plt.savefig(homedir + '/../plots/profiles/shape/v' + redvol + '/rmax' + maxrad + '/shape_v' + redvol + '_rmax' + maxrad + '_Ca' + capnum + '.pdf')
#plt.close()
plt.show()

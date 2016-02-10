import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from scipy import interpolate, interpolate
from math import log10, floor, pi

vmin = 0.6
vmax = 1
vrange = np.linspace(vmin,vmax,10000)
lamc = []
for v in vrange :
	a = 2.0*v
	b = -3.0
	c = 0
	d = 1.0
	coeff = [a,b,c,d]
	roots = np.roots(coeff)
	lc = max(roots)
	lamc.append(lc)
lamc = np.column_stack(lamc)
lamc = lamc[0]

# plot
fs = 10
plt.figure(figsize=(6,4))
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
plt.plot(vrange, lamc, 'k-')
plt.xlabel('$v$',fontsize=fs)
plt.ylabel('$\lambda_c$',fontsize=fs)
plt.xlim((vmin,vmax))
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show()


	


## set plot preferences
#fs = 24
#plt.xlabel('$x/R$', fontsize=fs+4)
#
#leg = plt.legend(('theory, ' + '$v$ = ' + strvL + ', $\epsilon$ = ' + streL + ', $\mathrm{Ca}$ = ' + strCL,),
#loc='best',fontsize=fs-4)
#leg.get_frame().set_linewidth(0.0)
#
##plt.savefig(homedir + '/../plots/profiles/shape/v' + redvol + '/rmax' + maxrad + '/shape_v' + redvol + '_rmax' + maxrad + '_Ca' + capnum + '.pdf')
##plt.close()
#plt.show()

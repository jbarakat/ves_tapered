def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

import numpy as np
from scipy import interpolate, integrate
from math import log10, floor, pi, sqrt
import os.path

altdir = os.path.abspath("../")
indir = os.path.abspath("../../input/raw")

load_src("readLT", altdir + "/tools/readLT.py")
from readLT import *

###################################################
#                     SETUP                       #
###################################################

# read parameters
ifile = open('./params.in', 'r')
print ifile
for line in ifile :
	params = line.split()		# space delimited text
	if params[0] == 'PROFILE' :
		STR = params[2]
	if params[0] == 'REDVOL' :
		redvol = params[2]
	if params[0] == 'CONFIN' :
		confin = params[2]
	if params[0] == 'CAPNUM' :
		capnum = params[2]
ifile.close()
	
###################################################
#               LUBRICATION THEORY                #
###################################################

# load data file
(v, conf, kb, area, vlme, t, s, x, r, nx, nr, tx, tr, cs, cphi, psi, qs, p, tau, gam, A, V, Q, S, Rtube, U, XCOM) = readInput(indir, redvol, confin, capnum)

#print 'PARAMETERS: v = ' + str(round(v,2)) + ', conf = ' + str(round(conf, 2)) + ', kb = ' + str(round(kb,5))
#print 'VARIABLES: t, r, x, psi, cs, qs, p, sig, A, V, Q, R, U, X, S'
#for i in range(len(s)) :
#	print '%6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f' % (t[i],r[i],x[i],psi[i],cs[i],qs[i],p[i],gam[i],A[i],V[i],Q,Rtube,U,XCOM,S)

f = open('input.dat', 'w')
f.write('PARAMETERS: v = ' + str(round(v,2)) + ', conf = ' + str(round(conf, 2)) + ', kb = ' + str(round(kb,5)) + '\n')
f.write('CHARACTERISTIC SCALES: pressure drop (dp = 1), suspending fluid viscosity (mu = 1), nominal vesicle radius (a = 1)\n')
f.write('VARIABLES: t, r, x, psi, cs, qs, p, sig, A, V, Q, R, U, X, S\n')
for i in range(len(s)) :
	f.write('%6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f\n' % (t[i],r[i],x[i],psi[i],cs[i],qs[i],p[i],gam[i],A[i],V[i],Q,Rtube,U,XCOM,S))
f.close()

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import glob
from scipy import interpolate, interpolate
from math import log10, floor, pi, sqrt
import os.path

homedir = os.path.abspath(os.path.join("./", os.pardir))
homedir = "../../output/v90/conf90/"

#outdir	= homedir + '/../output/v' + redvol + '/conf' + confin + '/*.dat'
#substring = homedir + '/../output/v' + redvol + '/conf' + confin + '/sln_v' + redvol + '_conf' + confin + '_Ca'
#files	= sorted(glob.glob(outdir))
#nfiles = len(files)
#lss = len(substring)

outdir = homedir + "*.dat"
substring = homedir + "sln_v90_conf90_Ca"
files	= sorted(glob.glob(outdir))
nfiles = len(files)
lss = len(substring)

print files[-1]

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from pylab import *
import numpy as np
from scipy import interpolate, interpolate
from math import log10, floor, pi
#from cmath import sqrt
import cmath

v = 0.8

#a1 = pow((-v - pow((v*v - 1.0),0.5)),2.0)
a1 = pow(-v - cmath.sqrt(v*v - 1.0),1.0/3.0)
a2 = pow(-v + cmath.sqrt(v*v - 1.0),1.0/3.0)

cf = (cmath.sqrt(3)/2.0)*1j

a = -0.5*(a1 + a2) + cf*(a1 - a2)
a = a.real

print a

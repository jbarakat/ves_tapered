from readProfiles import readProfiles
from plotProfiles import plotProfiles
from matplotlib import pyplot as plt

(STR,vL, epsL, CaL, xL, yL, vH, epsH, xH, yH) = readProfiles()
plotProfiles(STR, vL, epsL, CaL, xL, yL, vH, epsH, xH, yH)

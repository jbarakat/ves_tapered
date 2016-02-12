from readProfiles import readProfiles
from plotProfiles import plotProfiles
from matplotlib import pyplot as plt

(STR, v, eps, Ca, xL, yL, xH, yH) = readProfiles()
plotProfiles(STR, v, eps, Ca, xL, yL, xH, yH)

from readProfiles import readProfiles
from plotProfiles import plotProfiles

(redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau) = readProfiles()
plotProfiles(redvol, tubrad, flux, pdrop, s, x, r, p, gam, tau)

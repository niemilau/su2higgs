import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# from matplotlib.mlab import griddata
import scipy.interpolate
from scipy.interpolate import griddata
import scipy.optimize
import matplotlib
import sys
import string, math



# setup latex plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# make font size bigger
matplotlib.rcParams.update({'font.size': 14})
# but make legend smaller
matplotlib.rcParams.update({'legend.fontsize': 14})
matplotlib.rcParams['legend.numpoints'] = 1 # get rid of the annoying double-marker in legend
#plt.figure(figsize=(1,1))
#linestyles = ['_', '-', '--', ':']

from matplotlib import pyplot;
from pylab import genfromtxt;  


# read simulation results
dataphi = genfromtxt("phisq_collected.dat",names=True);
datasigma = genfromtxt("Sigmasq_collected.dat",names=True);

T = dataphi['T'];

# lattice field squared exp values of lattice fields, 
# defined as phiSq_lat = a * phiSq_cont + div
phisqlat = dataphi['phisqAvg']
sigmasqlat = datasigma['SigmasqAvg']
phiErr = dataphi['phisqErr']
sigmaErr = datasigma['SigmasqErr']

# remove lattice divergence: see Appendix A.2 in hep-lat/9504001 
# constants from lattice PT

Sigma = 3.1759115; 
delta = 1.942130;
zeta = 0.08849;

# read from data file (these are the same in both files):
g3sq = dataphi['gsq']  # units: GeV^1
RGscale = dataphi['RGscale']
a = dataphi['a']


# substract divergence. 
phisqCont = phisqlat / a  - Sigma/(2*math.pi*a) - 3*g3sq/(16*math.pi**2) * (np.log(6/(a*RGscale)) + zeta + 0.25 * Sigma**2 - delta );
sigmasqCont = sigmasqlat / a  - 3*Sigma/(4*math.pi*a) - 3*g3sq/(4*math.pi**2) * (np.log(6/(a*RGscale)) + zeta + 0.25 * Sigma**2 - delta );



phiErrCont = 0
sigmaErrCont = 0


# plot: error bars included. remember to change the errors to match what we plot
plt.errorbar(T, 2*phisqCont/T, yerr=2*phiErrCont, fmt='o', label = r"$2\phi^\dagger\phi / T$")
plt.errorbar(T, sigmasqCont/T, yerr=sigmaErrCont, fmt='o', label = r"$\Sigma^a\Sigma^a / T$")


# draw horizontal line at y = 0
plt.axhline(y=0, color='black', linestyle='-')

plt.title(r'$\beta_g = 8$. Note: fields are the 3d continuum fields',fontsize=12)
plt.xlabel(r'$T$')
#h = plt.ylabel(r'$2\langle\phi^\dagger\phi\rangle/T ')
#h.set_rotation(0)
#plt.xlim([0,0.2])
#plt.ylim([-0.5,2.0])
#plt.xticks(np.arange(min(T), max(T)+1, 1.0)) # set denser ticks
pyplot.legend();
pyplot.savefig('condensates.pdf')
plt.gcf().clear()


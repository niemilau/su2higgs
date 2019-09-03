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

matplotlib.rcParams['legend.numpoints'] = 1 
#plt.figure(figsize=(1,1))
#linestyles = ['_', '-', '--', ':']

from matplotlib import pyplot;
from pylab import genfromtxt;


# read simulation results
data1 = genfromtxt("h1_T.dat", names=True)
data2 = genfromtxt("h2_T.dat", names=True)

## same for both h1, h2
T = data1['T']


# lattice field exp values, in our case:
# h1 = \phi^\dagger\phi and h2 = 0.5 * \Sigma^a \Sigma^a
h1lat = data1['avg']
h2lat = data2['avg']

err1 = data1['error']
err2 = data2['error']

# remove lattice divergence: see hep-lat/9705003
# constants from lattice PT
Sigma = 3.1759115;
delta = 1.942130;
zeta = 0.08849;

# read required parameter from data file (these are the same in both h1, h2):
betaG = data1['betaG']; # SU2 beta
gsq = data1['gsq']; # 3d gauge coupling squared, units: GeV^1
RGscale = data1['RGscale'];

a = 4.0/(betaG * gsq); # lattice spacing

# vacuum counterterms (or their derivatives)
pi = math.pi
ct_phi = -2.0 * Sigma/(4.0*pi*a) - 3.0/(16.0*pi**2)*gsq*(np.log(6.0/(a*RGscale)) \
    + zeta + 0.25*Sigma**2 - delta)

ct_Sigma = -1.5 * Sigma/(4.0*pi*a) - 6.0/(16.0*pi**2)*gsq*(np.log(6.0/(a*RGscale)) \
    + zeta + 0.25*Sigma**2 - delta)


h1cont = h1lat / (1.0*a) + ct_phi;
h2cont = h2lat / (1.0*a) + ct_Sigma;
err1 = err1 / (1.0*a)
err2 = err2 / (1.0*a)

###############################################
######## Transform to 4d fields ###############
###############################################

## note: only works if there are equal number of rows in couplings4d as in input_couplings

# first read in the 4d couplings
#data4d = genfromtxt("couplings4d");
#ytsq = data4d[:,1];
#gsq = data4d[:,4];
#g1sq = data4d[:,5];
#Lb = data4d[:,14];
#Lf = data4d[:,15];

# field normalization relations from DR
#h1_4d = h1cont * T / (1 - (1/(16*math.pi**2)) * (0.75 * (3*gsq + g1sq) * Lb ));
#h2_4d = h2cont * T / (1 - (1/(16*math.pi**2)) * (0.75 * (3*gsq + g1sq) * Lb - 3*ytsq*Lf));

###############################################

# plot: error bars included. remember to change the errors to match what we plot
# also: multiply by 2 to compare with usual VEV
plt.errorbar(T, 2.0*h1cont/T, yerr=2.0*err1/T, fmt='o', label = r"2\langle\phi^\dagger\phi\rangle / T")
plt.errorbar(T, h2cont/T, yerr=err2/T, fmt='o', label = r"\langle\text{Tr}\Sigma^2\rangle / T")


# draw horizontal line at y = 0
plt.axhline(y=0, color='black', linestyle='-')

plt.title(r'$\beta_g = 16$',fontsize=12)
plt.xlabel(r'$T$ / GeV')
#h = plt.ylabel(r'$2\langle\phi^\dagger\phi\rangle/T ')
#h.set_rotation(0)
#plt.xlim([0,0.2])
#plt.ylim([-0.5,2.0])
#plt.xticks(np.arange(min(T), max(T)+1, 1.0)) # set denser ticks
plt.grid();
pyplot.legend();
pyplot.savefig('test_triplet.pdf')
plt.gcf().clear()

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

#plt.figure(figsize=(1,1))
#linestyles = ['_', '-', '--', ':']

from matplotlib import pyplot;
from pylab import genfromtxt;


# read simulation results
data1 = genfromtxt("h1_T.dat");
data2 = genfromtxt("h2_T.dat");

# first column is T
T = data1[:,0];

# lattice field squared exp values of lattice fields,
# defined as phiSq_lat = phiSq_cont / T + div
h1lat = data1[:,2];
h2lat = data2[:,2];

# remove lattice divergence: see hep-lat/9705003
# constants from lattice PT
Sigma = 3.1759115;
delta = 1.942130;
zeta = 0.08849;

# read from data file (these are the same in both h1, h2):
#betaG = data1[:,7]; # SU2 beta
#g3 = data1[:,8];  # units: GeV^1/2
#RGscale = data1[:,9];

#a = 4/(betaG * g3**2); # lattice spacing

# substract div: eq. (24) in the Rajantie&Laine paper
#h1cont = h1lat * T - Sigma/(2*math.pi*a) - 3*g3**2/(16*math.pi**2) * (np.log(6/(a*RGscale)) + zeta + 0.25 * Sigma**2 - delta );
#h2cont = h2lat * T - Sigma/(2*math.pi*a) - 3*g3**2/(16*math.pi**2) * (np.log(6/(a*RGscale)) + zeta + 0.25 * Sigma**2 - delta );

h1cont = h1lat*T;
h2cont = h2lat*T;

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

# FOR NOW ONLY PLOT THE CONTINUUM 3D FIELDS
# plot: error bars included. remember to change the errors to match what we plot
# also: multiply by 2 to compare with BG field thing
plt.errorbar(data1[:,0], 2*h1cont/T, yerr=2*data1[:,3], fmt='o', label = r"\phi_1^2")
plt.errorbar(data2[:,0], 2*h2cont/T, yerr=2*data2[:,3], fmt='o', label = r"\phi_2^2")


# draw horizontal line at y = 0
plt.axhline(y=0, color='black', linestyle='-')

plt.title(r'BM2: $\beta_g = 4$. Note: fields are the lattice fields',fontsize=12)
plt.xlabel(r'$T$')
h = plt.ylabel(r'$2\langle\phi^\dagger\phi\rangle/T ')
#h.set_rotation(0)
#plt.xlim([0,0.2])
#plt.ylim([-0.5,2.0])
#plt.xticks(np.arange(min(T), max(T)+1, 1.0)) # set denser ticks
pyplot.legend();
pyplot.savefig('test_triplet.pdf')
plt.gcf().clear()

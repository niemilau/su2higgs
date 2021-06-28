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
from pylab import genfromtxt;


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


kwargs = {'markeredgecolor':'black', 'markeredgewidth':0.5, 'markersize':5, 'elinewidth':1, 'capsize':3}

fig = plt.figure()
ax = fig.add_subplot(1,1,1)


# read simulation results
data_phisq = genfromtxt("h1_T.dat", names=True)
data_S = genfromtxt("h2_T.dat", names=True)
data_Ssq = genfromtxt("h3_T.dat", names=True)

## same for all observables
T = data_phisq['T']


# lattice field exp values
phisq_lat = data_phisq['avg']
S_lat = data_S['avg']
Ssq_lat = data_Ssq['avg']

phisq_err = data_phisq['err']
S_err = data_S['err']
Ssq_err = data_Ssq['err']

# remove lattice divergence: see hep-lat/9705003
# constants from lattice PT
Sigma = 3.1759115;
delta = 1.942130;
zeta = 0.08849;

# read required parameter from data file (same in all files):
betaG = data_phisq['betaG']; # SU2 beta
gsq = data_phisq['gsq']; # 3d gauge coupling squared, units: GeV^1
gpsq = data_phisq['gpsq']; # U(1) gauge coupling squared
RGscale = data_phisq['RGscale']; # MS-bar scale

a = 4.0/(betaG * gsq); # lattice spacing

logPlusZeta = np.log(6.0/(a*RGscale)) + zeta;

# derivatives of the vacuum counterterm. no b1 dependence and mSsq appears only at 1loop.
pi = math.pi

ct_phi = -2.0 * Sigma/(4.0*pi*a) \
    + 1.0/(16.0*pi**2)* (3*gsq + gpsq) * (-logPlusZeta + delta - Sigma**2 / 4.0)
    
# dV/d(1/2 m_S^2)
ct_SS = -Sigma/(4.0*pi*a)


phisq = phisq_lat / (1.0*a) + ct_phi;
phisq_err = phisq_err / (1.0*a)

Ssq = Ssq_lat / (1.0*a) + ct_SS;
Ssq_err = Ssq_err / (1.0*a)

S = S_lat / (1.0 * np.sqrt(a))
S_err = S_err / (1.0 * np.sqrt(a))


## plot: error bars included. multiply by 2 to compare with usual VEV
plt.errorbar(T, 2.0*phisq / T, yerr=2.0*phisq_err / T, fmt='o', label = r"$2\langle\phi^\dagger\phi\rangle / T$", color="green", **kwargs)
plt.errorbar(T, 2.0*Ssq / T, yerr=2.0*Ssq_err / T, fmt='o', label = r"$2\langle S^2 \rangle / T$", color="red", **kwargs)
plt.errorbar(T, S / np.sqrt(T), yerr= S_err / np.sqrt(T), fmt='o', label = r"$\langle S \rangle / \sqrt{T}$", color="blue", **kwargs)

# draw horizontal line at y = 0
plt.axhline(y=0, color='black', linestyle='-')

#plt.title(r'$\beta_G = 10, V = 24^3$. NB! starting from cold configuration',fontsize=12)
plt.xlabel(r'$T$ / GeV')
#h = plt.ylabel(r'$2\langle\phi^\dagger\phi\rangle/T ')
#h.set_rotation(0)
#plt.xlim([105,200])
#plt.ylim([-1.0,8.0])
#plt.xticks(np.arange(min(T), max(T)+1, 1.0)) # set denser ticks
plt.grid(linestyle=':', linewidth=0.75);
plt.legend();
plt.savefig('cond_singlet.pdf')
plt.gcf().clear()

## plot also susceptibility
h1susc = data1['susc']
h2susc = data2['susc']
suscErr1 = data1['suscErr']
suscErr2 = data2['suscErr']

plt.errorbar(T, h1susc, yerr=suscErr1, fmt='o', label = r"$\phi$")
plt.errorbar(T, h2susc, yerr=suscErr2, fmt='o', label = r"$\Sigma$")

plt.xlabel(r'$T$ / GeV')
plt.ylabel(r'$\chi$')
plt.grid();
pyplot.legend();
pyplot.savefig('susceptibility.pdf')
plt.gcf().clear()

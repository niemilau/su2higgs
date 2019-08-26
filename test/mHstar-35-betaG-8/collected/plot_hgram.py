import numpy as np
import os
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# from matplotlib.mlab import griddata
import scipy.interpolate
from scipy.interpolate import griddata
import scipy.optimize
import matplotlib
import sys
import string, math


lines = itertools.cycle(('--', ':','-.'))

# setup latex plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# make font size bigger
matplotlib.rcParams.update({'font.size': 26})
# but make legend smaller
matplotlib.rcParams.update({'legend.fontsize': 16})

#plt.figure(figsize=(1,1))
#linestyles = ['_', '-', '--', ':']

from matplotlib import pyplot;
from pylab import genfromtxt;

import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1


if not len(sys.argv) in [1,2]:
    sys.stderr.write('Usage: %s <datafile> \n' %
                     sys.argv[0])
    sys.stderr.write('If argument not specified, attempts to plot all .dat files in the folder.\n')

    sys.exit(1)

def plotHgram(datafile,beta, vol):
    data = genfromtxt(datafile);

    # first column is measured field squared value
    vev = data[:,0];
    # then its probability distribution + errorbars
    prob = data[:,1];
    # probErr = data[:,2];

    ## normalize the histogram by the area integral (total probability)
    dx = vev[1] - vev[0]
    norm = dx * sum(prob)
    prob = prob / norm

    # note that logs blow up if the hgram contains 0's, but python can handle these

    # some special cases for plot styling
    if (beta == 20):
        plt.plot(vev, np.log10(prob),dashes=[8,2,8,6,2,6],label=r'$\beta_G=$'+str(beta))
    elif (beta == 24):
        plt.plot(vev, np.log10(prob), linestyle='-',label=r'$\beta_G=$'+str(beta),linewidth=2.0, c='black')
    elif (beta == 32):
        plt.plot(vev, np.log10(prob), dashes=[12,6,3,6,3,6],label=r'$\beta_G=$'+str(beta),linewidth=1.0)
    else:
    #plt.plot(vev, np.log10(prob), linestyle=lines.next(),label=r'$\beta_G=$'+str(beta),dashes=[12,6,12,6,3,6])
        #plt.plot(vev, np.log10(prob), linestyle=lines.next(),label=r'$\beta_G=$'+str(beta)+", "+str(vol))
        plt.plot(vev, np.log10(prob), linestyle=lines.next(),label=str(vol))
    #plt.errorbar(vev, np.log10(prob), fmt='o', markersize=0, c="green")

if len(sys.argv) == 2:

    hgram_data = sys.argv[1]

    plotHgram(hgram_data, xx)



## else: plot all histograms in one picture
else:
    datafiles = []
    betas = []
    volumes = []
    for filename in os.listdir(os.path.dirname(os.path.realpath(__file__))):
        if filename.endswith(".dat"):
            # pick beta ad volume, assuming a certain pattern for file names
            smt = filename.split('beta');
            smt = smt[1].split('_')
            beta = smt[0]
            smt = filename.split('vol');
            smt = smt[1].split('.')
            vol = smt[0]

            datafiles.append(filename)
            betas.append(int(beta))
            volumes.append(vol)

    # sort lists
    files_sorted = []
    volumes_sorted = []
    betas.sort()
    for i in xrange(0,len(betas)):
        betastr = 'beta'+str(betas[i])
        for j in xrange(0,len(datafiles)):
            if betastr in datafiles[j] and datafiles[j] not in files_sorted:
                files_sorted.append(datafiles[j])
                volumes_sorted.append(volumes[j])


## filenames are now sorted by beta
for i in range(0,len(files_sorted)):
    plotHgram(files_sorted[i], betas[i], volumes_sorted[i])


plt.xlabel(r'$\langle \frac12 a\text{Tr} \Phi^\dagger\Phi \rangle$')
h = plt.ylabel(r'$log_{10}(P)$')
#h.set_rotation(0)
#plt.xlim([1.1,1.8])
#plt.ylim([-6,0])
plt.grid()
plt.tick_params(labelsize=20)
plt.xticks(np.arange(0, 5, 0.5)) # set denser ticks
pyplot.legend();
plt.tight_layout()
pyplot.savefig('hgrams.pdf')
plt.gcf().clear()

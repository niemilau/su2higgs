#!/usr/bin/env python

import sys
import numpy as np
from scipy.optimize import fsolve
import math
from pylab import genfromtxt


def nearest_indexOf(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index


def nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return array[index]


def lin_int(x, b, k):
	return b + k*x


## read in temperature T and beta_G
if not len(sys.argv) in [4,5]:
    sys.stderr.write('Usage: %s <input file> <T> <beta_G> <volume (optional)>\n' %
                     sys.argv[0])
    sys.exit(1)

datafile = sys.argv[1]
T = float(sys.argv[2])
beta = float(sys.argv[3])

if len(sys.argv) == 5:
    volume = int(sys.argv[4])

params = genfromtxt(datafile, names=True)

temps = params['T']


## find parameter labelled 'name' in the input file at given T
def get_param(name, T):
    if T not in temps:
        print("Error in get_param()!!\n")
        sys.exit(2)
    line = params[np.where(temps == T)]
    par = line[name][0]
    return par

## get all parameters
def getAllParams(T):

    gsq = get_param('gsq', T)
    gpsq = get_param('gpsq', T)
    mphisq = get_param('mphisq', T)
    lam = get_param('lam', T)
    mSsq = get_param('mSsq', T)
    b1 = get_param('b1', T)
    b3 = get_param('b3', T)
    b4 = get_param('b4', T)
    a1 = get_param('a1', T)
    a2 = get_param('a2', T)
    RGscale = get_param('RGscale', T)

    ## return as a list. Use this ordering from now on!!
    return [gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, RGscale]


## write params and T to a file
def write_params(fname, plist, T, islat):

	## if islat = 1, last parameter in plist is the scaling a, otherwise it is the MSbar scale
    if (islat):
        scale_id = 'spacing'
    else:
        scale_id = 'RGscale'

    [gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, scale] = plist

    f = open(fname, "w+")

    f.write('T %.16f\n' % T)
    f.write('%s %.16f\n' % (scale_id, scale))

    f.write('gsq %.16f\n' % gsq)
    f.write('gpsq %.16f\n' % gpsq)
    f.write('mphisq %.16f\n' % mphisq)
    f.write('lambda %.16f\n' % lam)
    f.write('mSsq %.16f\n' % mSsq)
    f.write('b1 %.16f\n' % b1)
    f.write('b3 %.16f\n' % b3)
    f.write('b4 %.16f\n' % b4)
    f.write('a1 %.16f\n' % a1)
    f.write('a2 %.16f\n' % a2)

    f.close()


#### convert to lattice parameters
def convert_lattice(p_cont):

    [gsq, gpsq, muphisq, lam, muSsq, b1, b3, b4, a1, a2, RGscale] = p_cont

    # generic constants
    Sigma = 3.17591153625
    zeta = 0.08849
    delta = 1.942130
    rho = -0.313964
    k1 = 0.958382
    k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
    k3 = 0.751498
    k4 = 1.204295

    a = 4.0/(gsq * beta)

    ## logarithm that appears frequently, plus zeta
    logPlusZeta = math.log(6.0/(a*RGscale)) + zeta

    ## TODO add U(1) contributions!!

    ## singlet tadpole term
    b1_ct1 = -Sigma/(4.0*a*math.pi) * (b3 + a1);
    b1_ct2 = (1.0)/(16.0*math.pi**2) * (3.0/8.0 * a1*gsq *(4.0*delta - Sigma**2) \
        + (2.0*b3*b4 + a1*a2 - 1.5*a1*gsq) * logPlusZeta);

    b1_lat = (b1 + b1_ct1 + b1_ct2) * a**(5.0/2.0)

    ## singlet mass counterterms
    mSsq_ct1 = -Sigma/(4.0*a*math.pi) * (2.0*a2 + 3.0*b4);
    mSsq_ct2 = (1.0)/(16.0*math.pi**2) * (0.75*a2*gsq * (4.0*delta - Sigma**2) \
        + (2.0*a2**2 + 6.0*b4**2 - 3.0*a2*gsq) * logPlusZeta);

    mSsq_lat = (muSsq + mSsq_ct1 + mSsq_ct2) * a**2.0

    ## Higgs mass counterterms. First contributions from SU(2) + Higgs only
    mphisq_ct1 = -Sigma/(4.0*a*math.pi) * (1.5*gsq + 6.0*lam);
    mphisq_ct2 = -1.0/(16*math.pi**2.0) * ( (51.0/16.0 * gsq**2 + 9.0*lam*gsq - \
        12.0*lam**2) * logPlusZeta + 9.0*lam*gsq * (0.25 *Sigma**2 - delta) \
        + 0.75*gsq**2 * (15.0/16.0 * Sigma**2 + math.pi/3.0 * Sigma + 1.25 \
        - 7.0/2.0 *delta - 4.0*rho + 4.0*k1 - k2 - k3 - 3.0*k4) );

    ## then contributions from the singlet
    mphisq_ct1 += -0.5*a2*Sigma/(4.0*a*math.pi);
    mphisq_ct2 += 1.0/(16*math.pi**2.0) * 0.5*a2**2 * logPlusZeta;


    mphisq_lat = (muphisq + mphisq_ct1 + mphisq_ct2) * a**2.0

    ## couplings: these are trivial
    gsq_lat = a * gsq
    gpsq_lat = a * gpsq
    lam_lat = a * lam
    b3_lat = b3 * a**(3.0/2.0)
    b4_lat = b4 * a
    a1_lat = a1 * a**(3.0/2.0)
    a2_lat = a2 * a

    return [gsq_lat, gpsq_lat, mphisq_lat, lam_lat, mSsq_lat, b1_lat, b3_lat, b4_lat, a1_lat, a2_lat, a]


####### end function definitions #########


##### Find input parameters at the given T, interpolate if needed
inputs = []
if T not in temps:
    print('\nNo direct match for the temperature found; interpolating...\n')
    index = nearest_indexOf(temps, T)
    T1 = temps[index-1]
    T2 = temps[index+1]
    p1 = getAllParams(T1)
    p2 = getAllParams(T2)
    for i in range(0, len(p1)):
        k = 1.0*(p2[i]-p1[i])/(T2 - T1)
        b = p1[i] - k * T1
        inputs.append( lin_int(T, b, k) ) ## uses same order as getAllParams()
else:
    inputs = getAllParams(T)


[gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, RGscale] = inputs

print
print('---- Read in continuum 3d parameters at T = '+str(T)+' ----')

print('gsq %g, gpsq %g, mphisq %g, lambda %g, mSsq %g, b1 %g, b3 %g, b4 %g, a1 %g, a2 %g, RGscale %g\n\n'
                     % (gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, RGscale))


## write these into a file
write_params('params_continuum.dat', inputs, T, 0)


## lattice parameters
paramsLat = convert_lattice(inputs)
[gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, a] = paramsLat

print('---- Lattice parameters for beta_G = '+str(beta)+' ----')

print('gsq %g, gpsq %g, mphisq %g, lambda %g, mSsq %g, b1 %g, b3 %g, b4 %g, a1 %g, a2 %g, a %g\n\n'
                     % (gsq, gpsq, mphisq, lam, mSsq, b1, b3, b4, a1, a2, a))


write_params('params_lattice.dat', paramsLat, T, 1)

### Calculate reweight string if the volume was given ###
if len(sys.argv) == 5:

    ## linearize in T
    T1 = nearest(temps, T-0.1)
    T2 = nearest(temps, T+0.1)

    if (T1 == T or T2 == T):
    	print('!!! Unable to calculate reweight string; include more values of T in the input file\n')
    	sys.exit(2)

    p1 = getAllParams(T1)
    p2 = getAllParams(T2)
    p1lat = convert_lattice(p1)
    p2lat = convert_lattice(p2)

    k = []
    for i in range(0, len(p1lat)):
    	k.append( 1.0*volume * (p2lat[i]-p1lat[i])/(T2-T1) )

    ## k now contains V * d/dT for each parameter in the action (V = number of sites):
    [gsq, gpsq, mphisq, lam, mSsq, b1, b2, b3, a1, a2, a] = k

    print('---- Reweight string ----\n')
    print('%.16f * #6 + %.16f * #7 + + %.16f * #8 + %.16f * #9 \
        + %.16f * #10 + %.16f * #11 + %.16f * #12 + %.16f * #13' %
        (mphisq, lam, b1, mSsq, b3, b4, a1, a2))

    #  #n is the nth column in measure file.

## end if len(sys.argv) == 5 ###

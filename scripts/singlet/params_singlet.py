#!/usr/bin/env python3

## This script reads in 3d continuum parameters from 'input file' at a given temperature,
## then converts and prints out the corresponding lattice parameters.
## Logic is such that spacing is fixed by specifying T and beta = 4/(a g_3^2(T)).
## If also the lattice volume is given, the script also calculates the 'reweight string',
## which is the temperature dependence of lattice action multiplied by volume (= number of lattice sites).
## Here also the T-dependence of beta is included, while spacing 'a' is kept fixed.

import sys
import numpy as np
import math
from pylab import genfromtxt


## Find index of nearest value in array
def nearest_indexOf(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index

## Find nearest value in array
def nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return array[index]


def lin_int(x, b, k):
	return b + k*x


## Find parameter labeled 'name' in the input array at given T. Will interpolate if needed
def get_param(paramList, name, T):
    temperatures = paramList['T']
    p = paramList[name]
    res = 0
    if T in temperatures:
        res = p[np.where(temperatures == T)][0]
    else: 
        ## linear interpolation around T
        i = nearest_indexOf(temperatures, T)
        T1 = temperatures[i-1]
        T2 = temperatures[i+1]
        p1 = p[np.where(temperatures == T1)][0]
        p2 = p[np.where(temperatures == T2)][0]
        k = (p2 - p1) / (T2 - T1)
        b = p1 - k * T1
        res = b + k * T
    return res


## get all parameters from the input file. Returns a dictionary
def getAllParams(paramList, T):
    
    res = {}
    res["T"] = T
    res["gsq"] = get_param(paramList, 'gsq', T)
    res["gpsq"] = get_param(paramList, 'gpsq', T)
    res["mphisq"] = get_param(paramList, 'mphisq', T)
    res["lambda"] = get_param(paramList, 'lam', T)
    res["mSsq"] = get_param(paramList, 'mSsq', T)
    res["b1"] = get_param(paramList, 'b1', T)
    res["b3"] = get_param(paramList, 'b3', T)
    res["b4"] = get_param(paramList, 'b4', T)
    res["a1"] = get_param(paramList, 'a1', T)
    res["a2"] = get_param(paramList, 'a2', T)
    res["RGscale"] = get_param(paramList, 'RGscale', T)

    return res

## Calculate dp/dT ~ (p2 - p1) / (T2 - T1) for each parameter p in the input dictionaries. 
## Assumes that both dictionaries contain element named 'T', the temperature
def CalcTDerivatives(paramDict1, paramDict2):
    dT = 0
    T1 = paramDict1["T"]
    T2 = paramDict2["T"]
    dT = T2 - T1
    if (dT == 0): 
        print("!!! Error in CalcTDerivatives(), invalid input given")
        exit(13)

    res = {}
    for key, value in paramDict1.items():
        p1 = value
        p2 = paramDict2[key] 
        res[key] = (p2 - p1) / dT

    return res



## Write a parameter dictionary to a file
def write_params(fname, pdict):

    with open(fname, 'w') as f: 
        for key, value in pdict.items(): 
            f.write('%s %s\n' % (key, value))


#### Convert to lattice parameters. p_cont and return value are dictionaries
def convert_lattice(p_cont, spacing, r_u1):

    gsq = p_cont["gsq"]
    gpsq = p_cont["gpsq"]
    muphisq = p_cont["mphisq"]
    lam = p_cont["lambda"]
    b1 = p_cont["b1"]
    muSsq = p_cont["mSsq"]
    b3 = p_cont["b3"]
    b4 = p_cont["b4"]
    a1 = p_cont["a1"]
    a2 = p_cont["a2"]
    RGscale = p_cont["RGscale"]

    # generic constants
    Sigma = 3.17591153625
    zeta = 0.08849
    delta = 1.942130
    rho = -0.313964
    k1 = 0.958382
    #k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
    #k3 = 0.751498
    k4 = 1.204295

    a = spacing

    ## logarithm that appears frequently, plus zeta
    logPlusZeta = math.log(6.0/(a*RGscale)) + zeta

    ### Counterterms. Here I assume the U(1) coupling gpsq to be in Y = 1/2 normalization


    ## Higgs mass counterterms. First contributions from SU(2) + Higgs only
    mphisq_ct1 = -Sigma/(8.0*math.pi*a) * (3*gsq + gpsq + 12*lam)

    mphisq_ct2 = 1.0/(16*math.pi**2.0) * ( \
        ( -51.0/16.0*gsq**2 + 9.0/8.0*gsq*gpsq \
        + 5.0/16.0*gpsq**2 - 3*lam * (3.0*gsq + gpsq) \
        + 12*lam**2) * logPlusZeta \
        + 3*lam * (3*gsq + gpsq) * (delta - 0.25*Sigma**2) \
        + gsq**2 * (-15.0/16.0 - 45.0/64.0 *Sigma**2 - math.pi/4.0*Sigma \
        + 33.0/8.0 *delta + 9.0/2.0*rho - 3*k1 + 3.0/2.0*k4) \
        + gpsq**2 * (1.0/16.0 - 1.0/64.0*Sigma**2 - math.pi*Sigma*r_u1**2 / 6.0 \
        + 1.0/8.0*delta + 0.5*rho) + gsq*gpsq * (3.0/8.0 - 3.0/32.0*Sigma**2 + 3.0/4.0*delta) \
        )

    ## then contributions from the singlet
    mphisq_ct1 += -a2*Sigma/(8.0*a*math.pi)
    mphisq_ct2 += 1.0/(16*math.pi**2.0) * 0.5*a2**2 * logPlusZeta;

    mphisq_lat = (muphisq + mphisq_ct1 + mphisq_ct2) * a**2


    ## singlet tadpole term
    b1_ct1 = -Sigma/(4.0*a*math.pi) * (b3 + a1)
    b1_ct2 = 1.0/(16.0*math.pi**2) * (  \
        (2.0*b3*b4 + a1*a2 - 0.5*a1* (3*gsq + gpsq)) * logPlusZeta \
        + 0.5*a1 * (3*gsq + gpsq)*(delta - Sigma**2 / 4.0) \
    )

    b1_lat = (b1 + b1_ct1 + b1_ct2) * a**(5.0/2.0)

    ## singlet mass counterterms
    mSsq_ct1 = -Sigma/(4.0*a*math.pi) * (2*a2 + 3*b4)
    mSsq_ct2 = 1.0/(16.0*math.pi**2) * ( \
        (2*a2**2 + 6*b4**2 - a2* (3*gsq + gpsq)) * logPlusZeta \
        + a2 * (3*gsq + gpsq)*(delta - Sigma**2 / 4.0) \
    )

    mSsq_lat = (muSsq + mSsq_ct1 + mSsq_ct2) * a**2.0


    ## couplings: these are trivial
    gsq_lat = a * gsq
    beta = 4.0 / gsq_lat
    gpsq_lat = a * gpsq
    betau1 = 4.0 / (gpsq_lat * r_u1**2) 
    ## factor 4 appears here because I want to convert continuum Higgs hypercharge Y = 1/2 to Y = 1 on lattice

    lam_lat = a * lam
    b3_lat = b3 * a**(3.0/2.0)
    b4_lat = b4 * a
    a1_lat = a1 * a**(3.0/2.0)
    a2_lat = a2 * a

    res = {
        "spacing" : a,
        "beta" : beta,
        "betaU1" : betau1,
        "r_U1" : r_u1,
        "mphisq" : mphisq_lat,
        "lambda" : lam_lat,
        "mSsq" : mSsq_lat,
        "b1" : b1_lat,
        "b3" : b3_lat,
        "b4" : b4_lat,
        "a1" : a1_lat,
        "a2" : a2_lat,
        "T" : p_cont["T"]
    }

    return res


####### Begin main #########

def main(): 
    ## read in temperature T and beta_G
    if not len(sys.argv) in [5,6]:
        sys.stderr.write('Usage: %s <input file> <T> <beta_G> <r_u1> <volume (optional)>\n' %
                        sys.argv[0])
        sys.exit(1)

    datafile = sys.argv[1]
    T = float(sys.argv[2])
    betaIn = float(sys.argv[3])
    r_u1 = float(sys.argv[4])

    if len(sys.argv) == 6:
        volume = int(sys.argv[5])

    if not (r_u1.is_integer()):
        print("!!! r_u1 = %.16f is not an integer, does not define U(1) irrep\n" % r_u1)

    paramList = genfromtxt(datafile, names=True)

    ## Will interpolate if the requested T was not found in file
    temperatures = paramList['T']
    bInterpolate = (T not in temperatures)

    if (bInterpolate):
        print('\n   ! No direct match for the temperature found; interpolating...\n')

    index = nearest_indexOf(temperatures, T)
    if (index-1 < 0 or index+1 >= len(temperatures)):
        print("!!! Found too few T-values around T=%g, cannot calculate T-derivatives. Use a file with broader T-range." % T)
        exit(11)

    ## Read the continuum parameters
    params_cont = getAllParams(paramList, T)

    print('\n---- Read in continuum 3d parameters at T = '+str(T)+' ----')
    print(params_cont)
    ## write these into a file
    write_params('params_continuum.dat', params_cont)


    ## Calculate lattice spacing here, keep fixed when computing T-derivatives
    gsq_cont = params_cont["gsq"]


    spacing = 4.0 / (betaIn * gsq_cont)

    
    ## lattice parameters
    params_lat = convert_lattice(params_cont, spacing, r_u1)


    print('\n---- Lattice parameters for beta_G = '+str(betaIn)+' ----')
    print(params_lat)
    write_params('params_lattice.dat', params_lat)


    print('\nScales (GeV): g^2 T %g, gT %g, a2 T %g, lattice cutoff %g' % (gsq_cont, math.sqrt(gsq_cont) * math.sqrt(T), params_cont["a2"], math.pi/spacing))

    ### Calculate reweight string if the volume was given ###
    if len(sys.argv) == 6:
        sidelen = volume**(1.0/3.0)
        print('Length scales (1/GeV): 1/(g^2 T) %g, lattice side %g' % (1.0/gsq_cont, spacing*sidelen))
        print('\n')

        ## First get continuum parameters at two nearby temperatures
        i = nearest_indexOf(temperatures, T)
        T1 = temperatures[i-1]
        T2 = temperatures[i+1]

        p1_cont = getAllParams(paramList, T1)
        p2_cont = getAllParams(paramList, T2)

        ## And convert them to lattice params
        p1_lat = convert_lattice(p1_cont, spacing, r_u1)
        p2_lat = convert_lattice(p2_cont, spacing, r_u1)

        dpdT = CalcTDerivatives(p1_lat, p2_lat)


        ## reweight string is dS/dT, but we measure volume averages instead of sums of type sum_x O(x) => multiple by volume here
        print('---- Reweight string ----\n')
        print( """%g * ( %.16g*#4 + %.16g*#5 + %.16g*#7 + %.16g*#8 +%.16g*#9 + %.16g*#10 + %.16g*#11 + %.16g*#12 + %.16g*#13 + %.16g*#14)""" % 
            (volume, dpdT["beta"], dpdT["betaU1"], dpdT["mphisq"], dpdT["lambda"], dpdT["b1"], dpdT["mSsq"], 
                dpdT["b3"], dpdT["b4"], dpdT["a1"], dpdT["a2"]))

        #  #n is the nth column in measure file.

    ## end if len(sys.argv) == 6 ###



main()
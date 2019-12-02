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


## eventually we will linearize the lattice parameters in T; used for latent heat and reweighting. this is done below but check here already that we have enough values

T1 = nearest(temps, T-0.5)
T2 = nearest(temps, T+0.5)

if (T1 == T or T2 == T):
	print('Linearization went out of bounds! Include more values of T in the input file.')
	print('Terminating...\n')
	sys.exit(2)



## pick the parameters corresponding to the given temperature

def getParams(T):
	if T in temps:
		inputs = params[np.where(temps == T)]
		gsq = inputs['gsq'][0]
		gpsq = inputs['gpsq'][0]
		muphisq = inputs['muphisq'][0]
		muSigmasq = inputs['muSigmasq'][0]
		lam = inputs['lam'][0]
		a2 = inputs['a2'][0]
		b4 = inputs['b4'][0]
		RGscale = inputs['RGscale3d'][0]
		return [T, gsq, gpsq, muphisq, muSigmasq, lam, a2, b4, RGscale]

	else:
		print('\nNo direct match for the temperature found; interpolating...\n')
		index = nearest_indexOf(temps, T)
		p1 = params[index]
		p2 = params[index+1]

		inputs = []
		for i in range(0, len(p1)):
			k = (p2[i]-p1[i])/(temps[index+1]-temps[index])
			b = p1[i] - k * temps[index]
			inputs.append( lin_int(T, b, k) )

		return inputs




# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
k3 = 0.751498
k4 = 1.204295



[T0, gsq, gpsq, muphisq, muSigmasq, lam, a2, b4, RGscale] = getParams(T)

print
print('---- Read in continuum 3d parameters at T = '+str(T)+' ----')
#print('Units of GeV are assumed (GeV^2 for mass squared).\n')

print('gsq %g, gpsq %g, muphisq %g, muSigmasq %g, lambda %g, b4 %g, a2 %g, RGscale %g\n\n'
                     % (gsq, gpsq, muphisq, muSigmasq, lam, b4, a2, RGscale))


## write these into a file
f = open('params_MSbar.dat',"w+")

f.write('T gsq gpsq muphisq muSigmasq lambda b4 a2 RGscale\n')

paramsMS = [T, gsq, gpsq, muphisq, muSigmasq, lam, b4, a2, RGscale]

for x in paramsMS:
	f.write(repr(x)+' ')

f.close()


## limitations on lattice size and spacing? correlation lengths \xi must satisfy
#   a << xi_min < xi_max << N*a, N = lattice side length
# assume here that Higgs v = T in the broken phase and calculate masses (tree level)
# note that this rough estimate only makes sense near the critical temperature

v = math.sqrt(T) # 3d units
mW = 0.5 * math.sqrt(gsq) * v
mh = math.sqrt(abs(2.0 * muphisq))
mSigma = math.sqrt(abs(muSigmasq + 0.5 * a2 * v**2))
## some symmetric phase masses
mSigmaSymm = math.sqrt(abs(muSigmasq))
masses = [mW, mh, mSigma, mSigmaSymm]
xi_min = 1.0/max(masses)
xi_max = 1.0/min(masses)

print('---- For phase transition with v ~ T, need:  ----')
print(' a << '+str(xi_min)+'; beta_G >> '+ str(4.0/(xi_min * gsq)) + ' and ' + 'L >> ' + str(xi_max / xi_min) +'\n')



#### convert to lattice parameters in David's convention

def convert_lattice(pars):
	[T0, gsq, gpsq, muphisq, muSigmasq, lam, a2, b4, RGscale] = pars

	a = 4.0/(gsq * beta)

	### lattice mass counterterms. NB! U(1) parts are missing from these

	deltamuphi1 = -(1.5*gsq + 6*lam) *float(Sigma)/(4.0*a*math.pi) \
        - 1.5*a2*float(Sigma)/(4.0*a*math.pi)


	deltamuphi2 = - 1.0/(16*math.pi**2.0) * ((51.0/16.0 * gsq**2 + 9.0*lam*gsq - \
        12.0*lam**2)*(math.log(6.0/(a*RGscale)) + zeta) + 9.0*lam*gsq * (0.25 *Sigma**2 - \
        delta) + 0.75*gsq**2 * (15.0/16.0 * Sigma**2 + math.pi/3.0 * Sigma + 1.25 - \
        7.0/2.0 *delta - 4.0*rho + 4.0*k1 - k2 - k3 - 3.0*k4) ) - 1.0/(16.0*math.pi**2) * \
        ( (-0.75*gsq**2 + 6.0*gsq*a2 - 1.5*a2**2) * (math.log(6.0/(a*RGscale)) + zeta) \
        + 6.0*gsq*a2 * (0.25*Sigma**2 - delta) - 3.0*gsq**2 * rho )


	mphisqLat = a**2 * (muphisq + deltamuphi1 + deltamuphi2)


	deltamuSigma1 = -(4.0*gsq + 5.0*b4 + 2.0*a2)*float(Sigma)/(4.0*math.pi*a)

	deltamuSigma2 = -1.0/(16.0*math.pi**2) * ( (20.0*b4*gsq - 10.0*b4**2) * \
        (math.log(6.0/(a*RGscale)) + zeta) + 20.0*b4*gsq * (0.25*Sigma**2 - delta) \
        + 2.0*gsq**2 *(1.25*Sigma**2 + math.pi/3.0 * Sigma - 6.0*delta - 6.0*rho + 4.0*k1 \
        - k2 - k3 - 3.0*k4) ) - 1.0/(16.0*math.pi**2) * ( (-gsq**2 + 3.0*a2*gsq - 2.0*a2**2) \
         * (math.log(6.0/(a*RGscale)) + zeta) + 3.0*a2*gsq * (0.25*Sigma**2 - delta) \
         - 4.0*gsq**2*rho )


	mSigmasqLat = a**2 * (muSigmasq + deltamuSigma1 + deltamuSigma2)


	## couplings: these are trivial

	gLat = math.sqrt(a * gsq)
	gpLat = math.sqrt(a * gpsq)
	lamLat = a * lam
	a2Lat = a * a2
	b4Lat = a * b4

	return [gLat, gpLat, mphisqLat, mSigmasqLat, lamLat, a2Lat, b4Lat, a]



[gLat, gpLat, mphisqLat, mSigmasqLat, lamLat, a2Lat, b4Lat, a] \
    = convert_lattice([T0, gsq, gpsq, muphisq, muSigmasq, lam, a2, b4, RGscale])

print('---- Lattice parameters for beta_G = '+str(beta)+' ----')

print('g %g, gp %g, muphisq %g, muSigmasq %g, lambda %g, b4 %g, a2 %g, a %g\n\n'
                     % (gLat, gpLat, mphisqLat, mSigmasqLat, lamLat, b4Lat, a2Lat, a))


#print('gsq: '+repr(gsqLat)+' gpsq: '+repr(gpsqLat)+' muphisq: '+repr(mphisqLat)+' muSigmasq: '+repr(mSigmasqLat)+' lambda: '+repr(lamLat)+' a2: '+repr(a2Lat)+' b4: '+repr(b4Lat)+' a: '+repr(a)\n)


## write these into a file
f = open('params_lattice.dat',"w+")

f.write('T beta g gp muphisq muSigmasq lambda b4 a2 a\n')

paramsLat = [T, beta, gLat, gpLat, mphisqLat, mSigmasqLat, lamLat, b4Lat, a2Lat, a]

for x in paramsLat:
	f.write(repr(x)+' ')

f.close()


#print('Written the lattice parameters at T = ' +str(T) + ' to params_lattice.dat.\n')


### Perform linearization of the lattice parameters.
# Note that the spacing a does NOT stay constant with the temperature if beta_G is fixed
p1 = getParams(T1)
p2 = getParams(T2)

p1lat = convert_lattice(p1)
p2lat = convert_lattice(p2)

k = []
for i in range(0, len(p1lat)):
	k.append( (p2lat[i]-p1lat[i])/(T2-T1) )


print('---- Temperature dependence of lattice parameters near T = '+str(T)+' ----\n')
print('g %g, gp %g, muphisq %g, muSigmasq %g, lambda %g, a2 %g, b4 %g, a %g\n\n'
                     % (k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]))


if len(sys.argv) == 5:
    ## if volume was given, print reweight string for Kari's FSH program.
    ## Note that we do multiply by volume here.
    ## Wilson action for SU(2) is not included, because beta_G does not depend on T by construction.
    ## also, more precision here just in case.
    print('---- Reweight string ----\n')
    print('%.16f * #5 + %.16f * #6 + %.16f * #8 + %.16f * #9 + %.16f * #10 \n' %
        (volume*k[2], volume*k[4], volume*k[3], volume*k[6], volume*k[5]))

    print('---- Same without spaces, for aa: ----\n')
    print('%.16f*#5+%.16f*#6+%.16f*#8+%.16f*#9+%.16f*#10\n' %
        (volume*k[2], volume*k[4], volume*k[3], volume*k[6], volume*k[5]))

# explanation: #n is the n.th column in measure file.
# So Kari's aa knows how to interpret this string

import sys
import numpy as np
from scipy.optimize import fsolve
import math
from pylab import genfromtxt


## read in temperature T and beta_G
if not len(sys.argv) in [6]:
    sys.stderr.write('Usage: %s <beta_G> <lambda> <msq> <beta_U1> <gamma_U1>\n' %
                     sys.argv[0])
    sys.exit(1)

beta = float(sys.argv[1])
lam = float(sys.argv[2])
msq = float(sys.argv[3])
betau1 = float(sys.argv[4])
gamma = float(sys.argv[5])

# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
k3 = 0.751498
k4 = 1.204295

print('---- Read in: beta_G = %g, lambda = %g, msq = %g, beta_U1 = %g, gamma_U1 = %g (gamma is in my convention!) ----\n'
    % (beta, lam, msq, betau1, gamma))

#### convert to x. Note that this is same in continuum and on (unimproved) lattice

gsq = 4.0/beta
x = lam / gsq
lam = x * gsq

if (betau1 != 0):
    gpsq = 1.0 / betau1
else:
    gpsq = 0

z = gpsq / gsq

RGscale = gsq

# The counterterms are written in my convention for gamma

# units of a^2 for the mass
mass_ct1 = -Sigma/(8.0*math.pi) * (3*gsq + gpsq + 12*lam);

mass_ct2 = 1.0/(16*math.pi**2.0) * ( \
    ( -51.0/16.0*gsq**2 + 9.0/8.0*gsq*gpsq \
    + 5.0/16.0*gpsq**2 - 3*lam * (3.0*gsq + gpsq) \
    + 12*lam**2)*(math.log(6.0/(RGscale)) + zeta) \
    + 3*lam * (3*gsq + gpsq) * (delta - 0.25*Sigma**2) \
    + gsq**2 * (-15.0/16.0 - 45.0/64.0 *Sigma**2 - math.pi/4.0*Sigma \
    + 33.0/8.0 *delta + 9.0/2.0*rho - 3*k1 + 3.0/2.0*k4) \
    + gpsq**2 * (1.0/16.0 - 1.0/64.0*Sigma**2 - 2*math.pi/(3.0*gamma**2)*Sigma \
    + 1.0/8.0*delta + 0.5*rho) + gsq*gpsq * (3.0/8.0 - 3.0/32.0*Sigma**2 + 3.0/4.0*delta) \
    )

## y in continuum:
y = 1.0*(msq - mass_ct1 - mass_ct2) / gsq**2

print('---- Continuum parameters: x = %g, y = %g, z = %g ----\n'
    % (x, y, z))

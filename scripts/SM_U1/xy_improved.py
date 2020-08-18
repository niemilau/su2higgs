import sys
import numpy as np
from scipy.optimize import fsolve
import math
from pylab import genfromtxt


## read in continuum x, y and beta_G (no U(1) implemented!!)
# these are also the "improved" lattice x, y, beta
if not len(sys.argv) in [4]:
    sys.stderr.write('Usage: %s <x> <y> <beta_G>\n' %
                     sys.argv[0])
    sys.exit(1)

x_imp = float(sys.argv[1])
y_imp = float(sys.argv[2])
beta_imp = float(sys.argv[3])
z = 0

print('---- Read in: x = %g, y = %g, z = %g, beta_G = %g ----\n'
    % (x_imp, y_imp, z, beta_imp))

## Calculate "naive" parameters in the lattice action
# see hep-lat/9709053, esp appendix B.1

ZgInv = 1.0 + 0.667394974 / beta_imp
Zm = 1 + (-0.23584293 + 0.29193980*x_imp) / beta_imp

beta = beta_imp * ZgInv
x = ZgInv * (x_imp + (0.01824624 - 0.47168586*x_imp + 0.58387959*x_imp**2) / beta_imp )
y = ZgInv**2 * Zm * y_imp

### now use naive lattice-continuum relations to calculate
## the lattice mass which should model the continuum theory with these naive values

# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
k3 = 0.751498
k4 = 1.204295

#### convert to lattice parameters, in units of a.
# i.e. the standard relations but treat a=1.

gsq = 4.0/beta
gpsq = gsq * z
lam = x * gsq

if (z != 0):
    betau1 = beta * (1.0 / z)
else:
    betau1 = 0

RGscale = gsq
# units of a^2 for the mass
mass_ct1 = -(1.5*gsq + 0.5*gpsq + 6.0*lam) *float(Sigma)/(4.0*math.pi)

mass_ct2 = - 1.0/(16*math.pi**2.0) * ((51.0/16.0 * gsq**2 - 9.0/8.0 * gsq*gpsq \
    - 5.0/16.0 * gpsq**2 + 3.0*lam * (3.0*gsq + gpsq) \
    - 12.0*lam**2)*(math.log(6.0/(RGscale)) + zeta) \
    + 3.0*lam * (3.0*gsq + gpsq) * (0.25 *Sigma**2 - delta) \
    + 0.75*gsq**2 * (15.0/16.0 * Sigma**2 + math.pi/3.0 * Sigma + 1.25 \
    - 7.0/2.0 *delta - 4.0*rho + 4.0*k1 - k2 - k3 - 3.0*k4) \
    - 0.9 * gsq * gpsq  + (0.01 + 1.7) * gpsq**2 + 1.7 * lam * gpsq  )

msq = gsq**2 * y + mass_ct1 + mass_ct2

print('---- Lattice parameters: msq = %g, lam = %g, betau1 = %g ----\n'
    % (msq, lam, betau1))

print('---- Naive lattice x, y, beta: x = %g, y = %g, beta = %g, Zg^(-2)*Zm = %g ----\n'
    % (x, y, beta, ZgInv**2 * Zm))

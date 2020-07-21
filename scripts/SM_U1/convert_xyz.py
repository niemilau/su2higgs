import sys
import numpy as np
from scipy.optimize import fsolve
import math
from pylab import genfromtxt


## read in temperature T and beta_G
if not len(sys.argv) in [5]:
    sys.stderr.write('Usage: %s <x> <y> <z> <beta_G>\n' %
                     sys.argv[0])
    sys.exit(1)

x = float(sys.argv[1])
y = float(sys.argv[2])
z = float(sys.argv[3])
beta = float(sys.argv[4])


# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
k3 = 0.751498
k4 = 1.204295

print('---- Read in: x = %g, y = %g, z = %g, beta_G = %g ----\n'
    % (x, y, z, beta))

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

print('---- Lattice parameters: msq = %.12lf, lam = %.12lf, betau1 = %.12lf ----\n'
    % (msq, lam, betau1))

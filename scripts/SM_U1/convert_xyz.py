import sys
import numpy as np
from scipy.optimize import fsolve
import math
from pylab import genfromtxt


## read in temperature T and beta_G
if not len(sys.argv) in [6]:
    sys.stderr.write('Usage: %s <x> <y> <z> <r_u1> <beta_G>\n' %
                     sys.argv[0])
    sys.exit(1)

x = float(sys.argv[1])
y = float(sys.argv[2])
z = float(sys.argv[3])
r_u1 = float(sys.argv[4])
beta = float(sys.argv[5])


# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
#k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
#k3 = 0.751498
k4 = 1.204295

print('---- Read in: x = %g, y = %g, z = %g, r_u1 = %g, beta_G = %g ----\n'
    % (x, y, z, r_u1, beta))

if not (r_u1.is_integer()):
    print("!!! r_u1 = %.16f is not an integer, does not define U(1) irrep\n" % r_u1)

#### convert to lattice parameters, in units of a.
# i.e. the standard relations but treat a=1.

## z = gpsq / gsq, where the U(1) coupling gpsq is in Y = 1/2 normalization (Higgs hypercharge)

gsq = 4.0/beta
gpsq = gsq * z
lam = x * gsq

if (z != 0):
    betau1 = 4.0 / (gpsq * r_u1**2)
else:
    betau1 = 0

RGscale = gsq

# units of a^2 for the mass
mass_ct1 = -Sigma/(8.0*math.pi) * (3*gsq + gpsq + 12*lam);

mass_ct2 = 1.0/(16*math.pi**2.0) * ( \
    ( -51.0/16.0*gsq**2 + 9.0/8.0*gsq*gpsq \
    + 5.0/16.0*gpsq**2 - 3*lam * (3.0*gsq + gpsq) \
    + 12*lam**2)*(math.log(6.0/(RGscale)) + zeta) \
    + 3*lam * (3*gsq + gpsq) * (delta - 0.25*Sigma**2) \
    + gsq**2 * (-15.0/16.0 - 45.0/64.0 *Sigma**2 - math.pi/4.0*Sigma \
    + 33.0/8.0 *delta + 9.0/2.0*rho - 3*k1 + 3.0/2.0*k4) \
    + gpsq**2 * (1.0/16.0 - 1.0/64.0*Sigma**2 - math.pi*Sigma*r_u1**2 / 6.0 \
    + 1.0/8.0*delta + 0.5*rho) + gsq*gpsq * (3.0/8.0 - 3.0/32.0*Sigma**2 + 3.0/4.0*delta) \
    )

msq = gsq**2 * y + mass_ct1 + mass_ct2

print('---- Lattice parameters: msq = %.12lf, lam = %.12lf, betau1 = %.12lf, r_u1 = %.12lf ----\n'
    % (msq, lam, betau1, r_u1))

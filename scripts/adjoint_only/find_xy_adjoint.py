#!/usr/bin/env python

import sys
import numpy as np
from scipy.optimize import fsolve
import math


## read in lattice parameters msq, b4, beta_G (all dimensionless)
if not len(sys.argv) in [4]:
    sys.stderr.write('Usage: %s <msq> <b4> <beta_G> (all dimensionless lattice parameters)\n' %
                     sys.argv[0])
    sys.exit(1)

msq = float(sys.argv[1]) # msq = a^2 (g^4 y + ct)
b4 = float(sys.argv[2]) # b4 = a g^2 x 
beta = float(sys.argv[3]) # beta = 4/(a g^2)

print
print('---- Read in parameters msq=%g, b4=%g, beta_G=%g ----\n\n' % (msq, b4, beta))


# generic constants
Sigma = 3.17591153625
zeta = 0.08849
delta = 1.942130
rho = -0.313964
k1 = 0.958382
k2 = 0.25*Sigma**2 - 0.5 * delta - 0.25
k3 = 0.751498
k4 = 1.204295

## lattice g_lat^2 = a g^2
gsq = 4.0/beta

RGscale = gsq ## MSbar scale = continuum g^2 (logarithm below has a*g^2)

## Counterterms for the adjoint Higgs
# these are now in units of a^2 (so just change continuum coupling -> lattice coupling)

msq_ct1 = -(4.0*gsq + 5.0*b4)*Sigma/(4.0*math.pi)

## 2-loop contribution from SU(2) + adjoint Higgs
msq_ct2 = -1.0/(16.0*math.pi**2) * ( (20.0*b4*gsq - 10.0*b4**2) * \
    (math.log(6.0/(RGscale)) + zeta) + 20.0*b4*gsq * (0.25*Sigma**2 - delta) \
    + 2.0*gsq**2 *(1.25*Sigma**2 + math.pi/3.0 * Sigma - 6.0*delta - 6.0*rho + 4.0*k1 \
    - k2 - k3 - 3.0*k4) )


## 
y = (msq - msq_ct1 - msq_ct2) / (gsq**2)

x = b4 / gsq
 

print('---- Continuum parameters for beta_G=%g ----' % (beta))
print('x %g   y %g' % (x, y) )



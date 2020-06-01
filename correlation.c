/** @file correlation.c
*
* Routines for calculating correlation lengths.
* This is done along a fixed direction ("z direction"),
* so the flag MEASURE_Z should be specified in Makefile.
*
* //TODO blocking
*/

#ifdef MEASURE_Z

#include "su2.h"

/* Calculate H(l) = (1/V) \sum_z h(z) h(z+l),
* where h(z) = \sum_{x,y} Tr\Phi(z)\Phi(z)
*/
double higgs_correlator(params* p, fields* f, long z) {

}

/* Notes for blocking (note: need side lengths to be even):
* 1. can create a new field struct f_block and copy the fields there.
*    this would be the k=0 blocking.
* 2. for k=1 blocking, use eqs 6.4 - 6.6 in hep-lat/9612006.
*    can probably use routines already defined in su2u1.c to calculate these.
     The blocked fields are effectively defined on a lattice where spacing
     in xy plane is doubled.
    probably need to create new lattice tables too...
    move these all to a new struct "lattice"?
*/


#endif

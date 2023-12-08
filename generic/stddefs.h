#ifndef STDDEFS_H
#define STDDEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#include "mersenne.h"

// Globals
int myRank; // MPI rank. Also 'rank' in lattice struct 
int MPISize; // also 'l.size'

/* Check compatibility of makefile flags */

#if (NHIGGS > 2)
	#warning !!! Higgs potential not implemented for N>2 doublets !!!
#elif (NHIGGS > 1) && defined (TRIPLET)
	#warning !!! TRIPLET with N>1 Higgs doublets not implemented !!!
#elif (NHIGGS > 1) && defined (SINGLET)
	#warning !!! SINGLET with N>1 Higgs doublets not implemented !!!
#elif ( (NHIGGS > 0) || defined (SINGLET)) && defined (GRADFLOW)
	#warning !!! Gradient flow not implemented for all fields !!!
#elif ( (NHIGGS > 0) || defined (SINGLET)) && defined (BLOCKING)
	#warning !!! Blocking not implemented for all fields !!!
#endif

/* If using blocking, need larger halos because of link smearing in su2u1.c */
#ifdef BLOCKING
	#define HALOWIDTH 2

#else
 #define HALOWIDTH 1 // no blocking, halo extends just one site in each direction

#endif

// degrees of freedom per site for different fields
#define SU2DB 4
#define SU2LINK 4
#define SU2TRIP 3

// update algorithms
#define METROPOLIS 1 // Metropolis
#define HEATBATH 2 // Heatbath
#define OVERRELAX 3 // Overrelaxation

// parity identifiers
#define EVEN 0
#define ODD 1
#define EVENODD 2

// multicanonical order parameters
#define PHISQ 1
#define SIGMASQ 2
#define PHI2MINUSSIGMA2 3
#define PHI2SQ 4

// different modes for updating the multicanonical weight function
#define READONLY 0
#define FAST 1
#define SLOW 2 // these are protected by the SLOW_MUCA flag

// Some convenient global variables:
// Keeping track of evaluation time
double waittime;
double Global_comms_time, Global_total_time;
double Global_current_action; // for debugging gradient flows etc

// complex numbers
typedef struct {
	double re, im;
}	complex;


/* inlines */

// Flip parity
inline char otherparity(char parity) {
	if (parity == EVEN)
		return ODD;
	else if (parity == ODD) {
		return EVEN;
	} else {
		printf("!!! Error in function otherparity !!!\n");
	}
}

// multiply two complex numbers
inline complex cmult(complex z1, complex z2) {
	complex res;
	res.re = z1.re*z2.re - z1.im*z2.im;
	res.im = z2.re*z1.im + z1.re*z2.im;
	return res;
}

#endif // ifndef STDDEFS_H

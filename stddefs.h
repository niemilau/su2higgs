#ifndef STDDEFS_H
#define STDDEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#ifdef MPI
	#include <mpi.h>
#else
	// No MPI, define dummy communicator (not actually used in serial)
	typedef struct {} MPI_Comm;
#endif


/* Check compatibility of makefile flags */

#if (NHIGGS > 2)
	#warning !!! Higgs potential not implemented for N>2 doublets !!!
#elif (NHIGGS > 1) && defined (TRIPLET)
	#warning !!! TRIPLET with multiple Higgs doublets not implemented !!!
#elif (NHIGGS > 0) && defined (GRADFLOW)
	#warning !!! Gradient flow with Higgs not implemented !!!
#elif (NHIGGS > 0) && defined (BLOCKING)
	#warning !!! Blocking not implemented for the Higgs !!!
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
#define PHI2SIGMA2 4

// nasty globals for keeping track of evaluation time
double waittime;
double Global_comms_time, Global_total_time;
double Global_current_action; // for debugging gradient flows etc

// complex numbers
typedef struct {
	double re, im;
}	complex;


/* inlines */
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

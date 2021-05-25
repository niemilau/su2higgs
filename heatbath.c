/** @file heatbath.c
*
* Routines for implementing Kennedy-Pendleton heatbath algorithm
*	for gauge links (Phys.Lett. 156B (1985) 393-399).
*
*/

#include "su2.h"

/*
* Update a single SU(2) link using KP heatbath.
*/
int heatbath_su2link(lattice const* l, fields* f, params const* p, long i, int dir) {

	// only used if accept/reject step is needed
	double oldlink[4];
	double oldact = 0.0;

	#ifdef TRIPLET
		oldlink[0] = f->su2link[i][dir][0];
		oldlink[1] = f->su2link[i][dir][1];
		oldlink[2] = f->su2link[i][dir][2];
		oldlink[3] = f->su2link[i][dir][3];
		oldact += hopping_triplet_forward(l, f, p, i, dir);
	#endif


	// calculate the staple and its determinant
	double V[4];
	su2link_staple(l, f, p, i, dir, V);

	double a = sqrt( su2sqr(V) );

	// normalize V to produce an SU(2) matrix
	for (int k=0; k<4; k++) {
		V[k] *= -1.0/a; // minus sign here because the weight is exp(-Tr U.S), S=full staple
	}

	// normalization factor alpha for K-P algorithm.
	// factor of 2 here since Tr U = 2 a_0
	a *= 2.0;

	double a0 = 0.0;

	double r1,r2,r3,d;

	int loop = 1;
	int maxloops = 200;
	do {
		r1 = -1.0*log(1.0 - dran())/a; //dran() is between [0.0, 1.0)
	  r2 = -1.0*log(1.0 - dran())/a;
	  r3 = cos(2*M_PI*dran());
	  r3 *= r3;
	  d = r1*r3 + r2; // this is the delta in K-P paper
		r3 = 1.0 - dran();

		loop++;
	} while (r3*r3 > 1 - 0.5*d && loop <= maxloops);

	if (loop > maxloops) {
		fprintf(stderr,
	  "Exceeded loop limit in SU(2) heatbath update! Was %d\n", maxloops);
		a0 = 1e-9;
	} else {
		// now we have a[0]:
		a0 = 1.0 - d;
	}

	// generate uniform a[1],a[2],a[3] on S_2 sphere with radius sqrt(1 - a0^2)
	double rad = 1.0 - a0*a0;
	if (rad < 0.0) {
		fprintf(stderr,
	  "Negative radius in SU(2) heatbath update, value %g\n", rad);
		return 0;
	}
	rad = sqrt(rad);

	// careful with the a[i] here, they can be zero
	d = 0.0;
	while (!(d > 0.0)) {
		r1 = 1.0 - 2.0*dran();
		r2 = 1.0 - 2.0*dran();
		r3 = 1.0 - 2.0*dran();
		d = sqrt(r1*r1 + r2*r2 + r3*r3);
	}
	d = 1.0/d;


	f->su2link[i][dir][0] = a0;
	f->su2link[i][dir][1] = r1*rad*d;
	f->su2link[i][dir][2] = r2*rad*d;
	f->su2link[i][dir][3] = r3*rad*d;

	// new link value is obtained by rotating this from the right
	// with hermitian conjugate of the normalized staple.
	V[1] = -V[1];
	V[2] = -V[2];
	V[3] = -V[3];
	su2rot(f->su2link[i][dir], V);

	#ifdef TRIPLET
		double newact = 0.0;
		newact += hopping_triplet_forward(l, f, p, i, dir);
		double diff = oldact - newact;
		if (dran() < exp(diff)) {
			return 1;
		} else {
			f->su2link[i][dir][0] = oldlink[0];
			f->su2link[i][dir][1] = oldlink[1];
			f->su2link[i][dir][2] = oldlink[2];
			f->su2link[i][dir][3] = oldlink[3];
			return 0;
		}
	#endif

	return 1;

}

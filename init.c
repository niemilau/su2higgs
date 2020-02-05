/** @file init.c
*
* Routines for initializing a simulation.


//TODO
*
*/

#include "su2.h"


/*
****************** Field initializations ******************
*/

// set SU(2) links to unity
void setsu2(fields f, params p) {

	for (long i=0; i<p.sites_total; i++) {
		for (int dir=0; dir<p.dim; dir++) {
			f.su2link[i][dir][0] = 1.0;
			f.su2link[i][dir][1] = 0.0;
			f.su2link[i][dir][2] = 0.0;
			f.su2link[i][dir][3] = 0.0;
		}
	}
}

// set U(1) links to unity
void setu1(fields f, params p) {

	for (long i=0; i<p.sites_total; i++) {
		for (int dir=0; dir<p.dim; dir++) {
			f.u1link[i][dir] = 0.0;
		}
	}
}

/* generate a random SU(2) matrix and store it in the argument
* Original implementation by David Weir
*/
void random_su2link(double *su2) {

	double u[4];
	// this number can be adjusted for e.g. metro update
	u[0] = 10.0 * (drand48() - 0.5);
	u[1] = 1.0 * (drand48() - 0.5);
	u[2] = 1.0 * (drand48() - 0.5);
	u[3] = 1.0 * (drand48() - 0.5);

	double norm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);

	su2[0] = u[0]/norm;
	su2[1] = u[1]/norm;
	su2[2] = u[2]/norm;
	su2[3] = u[3]/norm;
}

// Initialize SU(2) doublets to unity * 1/sqrt(2).
void setdoublets(fields f, params p) {
	for (long i=0; i<p.sites_total; i++) {
		f.su2doublet[i][0] = p.phi0;
		f.su2doublet[i][1] = 0.0;
		f.su2doublet[i][2] = 0.0;
		f.su2doublet[i][3] = 0.0;
	}
}


// Initialize SU(2) adjoint scalars
void settriplets(fields f, params p) {
	for (long i=0; i<p.sites_total; i++) {
		f.su2triplet[i][0] = p.sigma0;
		f.su2triplet[i][1] = 0.0;
		f.su2triplet[i][2] = 0.0;
	}
}


// Initialize all fields used in the simulation
void setfields(fields f, params p) {

	setsu2(f, p);
	#ifdef U1
		setu1(f, p);
	#endif

	#ifdef HIGGS
		setdoublets(f, p);
	#endif
	#ifdef TRIPLET
		settriplets(f, p);
	#endif
}


// Initialize accept/reject counters and time
void init_counters(counters* c) {

	c->comms_time = 0.0;
	c->total_time = 0.0;
	c->iter = 1;

	c->higgs_sweeps = 0;
	c->triplet_sweeps = 0;
	c->accepted_su2link = 0;
	c->accepted_u1link = 0;
	c->accepted_doublet = 0;
	c->accepted_triplet = 0;
	c->acc_overrelax_doublet = 0;
	c->acc_overrelax_triplet = 0;
	c->total_su2link = 0;
	c->total_u1link = 0;
	c->total_doublet = 0;
	c->total_triplet = 0;
	c->total_overrelax_doublet = 0;
	c->total_overrelax_triplet = 0;
	c->accepted_muca = 0;
	c->total_muca = 0;
}

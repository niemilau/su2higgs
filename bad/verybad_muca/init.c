/** @file init.c
*
* Routines for initializing a simulation.


//TODO
*
*/

#include "su2.h"



/*
* Function that loops over the lattice index and stores the neighboring site indices
* of each site into p.next and p.prev. Also calculate parity for the sites.
*/
/* OUTDATED
void calculate_neighbors(params* p) {

	uint* x = malloc(p->dim * sizeof(x));

	if (!p->rank)
		printf("Constructing lookup table for lattice sites...\n");

	for (long i=0; i<p->vol; i++) {

		// get the physical coordinates of site i and store in x
		indexToCoords(*p, i, x);
		// calculate and store the parity of site i
		uint coord = 0;
		for (int j=0; j<p->dim; j++) {
			coord += x[j];
		}
		p->parity[i] = coord % 2;


		for (int dir=0; dir < p->dim; dir++) {
			x[dir]++;
			p->next[i][dir] = coordsToIndex(*p,x);
			x[dir] -= 2;
			p->prev[i][dir] = coordsToIndex(*p,x);
			// return to the original value
			x[dir]++;
		}

		#ifdef DEBUG
		if (i != coordsToIndex(*p, x)) {
			printf("Sanity check failed for site i = %lu!! Coordinates do not reproduce the correct site index!\n", i);
		}
		#endif
	}

	free(x);
	if (!p->rank)
		printf("Lookup table constructed succesfully.\n");
}

*/



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
	if (!p.rank)
		printf("SU(2) links set to unity.\n");
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
	
	c->accepted_su2link = 0;
	c->accepted_doublet = 0;
	c->accepted_triplet = 0;
	c->acc_overrelax_doublet = 0;
	c->acc_overrelax_triplet = 0;
	c->total_su2link = 0;
	c->total_doublet = 0;
	c->total_triplet = 0;
	c->total_overrelax_doublet = 0;
	c->total_overrelax_triplet = 0;
	c->accepted_muca = 0;
	c->total_muca = 0;
}

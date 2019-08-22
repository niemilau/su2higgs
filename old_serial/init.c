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
void calculate_neighbors(params* p) {

	uint* x = malloc(p->dim * sizeof(x));

	printf("Constructing lookup table for lattice sites...\n");
	
	for (ulong i=0; i<p->vol; i++) {

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

	printf("Lookup table constructed succesfully.\n");
}


/*
* 	Routines for fetching neighbours on a periodic lattice
*/

// Calculate product L[0]*L[1]*...*L[max-1] of lattice lengths
inline ulong Lprod(params p, int max) {
	ulong L = 1;

	if (!(p.dim < max)) {
		for (int i=0; i<max; i++) { L *= p.L[i]; }
	}

	return L;
}

// From the site index, get the physical coordinates (x, y, z,...) and store in x
// see Mathematica notebook coordinates.nb for analytical relations
void indexToCoords(params p, ulong i, uint* x) {

	x[0] = i % p.L[0];

	// we want the floor() of integer division here. This is automatic in C, but I'm being explicit here.
	// possible problem here: floor(0.9999) returns 0, when the analytical value should be 1?
	x[p.dim-1] =  (int)floor(i / (Lprod(p, p.dim-1)) );


	for (int dir=p.dim-2; dir>0; dir--) {
		int a=0;
			for (int k=dir+1; k<p.dim; k++) {
				a += x[k] * Lprod(p, k);
  		}
		x[dir] = (uint)floor((i - a) / Lprod(p, dir));
	}


}


// from coordinates (x, y, z...) get the site index
ulong coordsToIndex(params p, uint* x) {

		ulong res = 0;
		for (int dir=0; dir<p.dim; dir++) {
			res += ( (x[dir] + p.L[dir]) % p.L[dir] ) * Lprod(p, dir);
		}
		return res;
}


/*
****************** Field initializations ******************
*/

// set SU(2) links to unity
void setsu2(fields f, params p) {

	for (ulong i=0; i<p.vol; i++) {
		for (int dir=0; dir<p.dim; dir++) {
			f.su2link[i][dir][0] = 1.0;
			f.su2link[i][dir][1] = 0.0;
			f.su2link[i][dir][2] = 0.0;
			f.su2link[i][dir][3] = 0.0;
		}
	}
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
	for (ulong i=0; i<p.vol; i++) {
		f.su2doublet[i][0] = p.phi0;
		f.su2doublet[i][1] = 0.0;
		f.su2doublet[i][2] = 0.0;
		f.su2doublet[i][3] = 0.0;
	}
}


// Initialize SU(2) adjoint scalars
void settriplets(fields f, params p) {
	for (ulong i=0; i<p.vol; i++) {
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

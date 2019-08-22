/** @file metropolis.c
*
* Routines for implementing the Metropolis algorithm for fields.
*	Actual updating of the lattice is performed in update.c
*
* TODO
*
*/

#include "su2.h"

/*
* Update a single SU(2) link using Metropolis.
* Returns 1 if update was accepted and 0 if rejected.
*/
int metro_su2link(fields f, params p, long i, int dir) {

	double oldlink[4];
	memcpy(oldlink, f.su2link[i][dir], SU2LINK*sizeof(double));

	double linkact_old = localact_su2link(f, p, i, dir);

	// generate new link randomly
	double newlink[4];
	random_su2link(newlink);

	su2rot(f.su2link[i][dir], newlink);

	double linkact_new = localact_su2link(f, p, i, dir);

	double diff = linkact_new - linkact_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > drand48() )) {
		return 1;
	}
	else {
		memcpy(f.su2link[i][dir], oldlink, SU2LINK*sizeof(double));
		return 0;
	}
}


/*
* Update a single SU(2) scalar doublet using Metropolis.
* If transverse == 1, forces a transverse update, i.e.
* keeps Tr \Phi^\dagger \Phi constant. This is useful
* for improving ergodicity especially for multicanonical
* runs that use phisq as the order parameter.
* Returns 1 if update was accepted and 0 if rejected.
*/
int metro_doublet(fields f, params p, long i) {

	double oldfield[4];
	oldfield[0] = f.su2doublet[i][0];
	oldfield[1] = f.su2doublet[i][1];
	oldfield[2] = f.su2doublet[i][2];
	oldfield[3] = f.su2doublet[i][3];


	double act_old = localact_doublet(f, p, i);

	// modify the old field by random values
	f.su2doublet[i][0] += 1.0*(drand48() - 0.5);
	f.su2doublet[i][1] += 1.0*(drand48() - 0.5);
	f.su2doublet[i][2] += 1.0*(drand48() - 0.5);
	f.su2doublet[i][3] += 1.0*(drand48() - 0.5);


	double act_new = localact_doublet(f, p, i);

	int accept;
	double diff = act_new - act_old;
	if (diff < 0) {
		accept = 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > drand48() )) {
		accept = 1;
	}
	else {
		f.su2doublet[i][0] = oldfield[0];
		f.su2doublet[i][1] = oldfield[1];
		f.su2doublet[i][2] = oldfield[2];
		f.su2doublet[i][3] = oldfield[3];
		accept = 0;
	}

	return accept;
}


/*
* Update a single SU(2) scalar triplet using Metropolis.
* If transverse == 1, forces a transverse update, i.e.
* keeps Tr A^2 constant.
* Returns 1 if update was accepted and 0 if rejected.
*/
int metro_triplet(fields f, params p, long i) {

	double oldfield[3];
	oldfield[0] = f.su2triplet[i][0];
	oldfield[1] = f.su2triplet[i][1];
	oldfield[2] = f.su2triplet[i][2];

	double act_old = localact_triplet(f, p, i);

	// modify the old field by random values
	f.su2triplet[i][0] += 1.0*(drand48() - 0.5);
	f.su2triplet[i][1] += 1.0*(drand48() - 0.5);
	f.su2triplet[i][2] += 1.0*(drand48() - 0.5);


	double act_new = localact_triplet(f, p, i);

	double diff = act_new - act_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > drand48() )) {
		return 1;
	}
	else {
		f.su2triplet[i][0] = oldfield[0];
		f.su2triplet[i][1] = oldfield[1];
		f.su2triplet[i][2] = oldfield[2];
		return 0;
	}
}

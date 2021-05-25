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
* Returns 1 if update was accepted and 0 if rejected. */
int metro_su2link(lattice const* l, fields* f, params const* p, long i, int dir) {

	double oldlink[4];
	memcpy(oldlink, f->su2link[i][dir], SU2LINK*sizeof(double));

	double linkact_old = localact_su2link(l, f, p, i, dir);

	// generate new link randomly
	double newlink[4];
	random_su2link(newlink);

	su2rot(f->su2link[i][dir], newlink);

	double linkact_new = localact_su2link(l, f, p, i, dir);

	double diff = linkact_new - linkact_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > dran() )) {
		return 1;
	}
	else {
		memcpy(f->su2link[i][dir], oldlink, SU2LINK*sizeof(double));
		return 0;
	}
}

/*
* Update a single U(1) link using Metropolis.
* Remember that our links are U_i(x) = exp(i a_i(x))
* and a_i(x) is in f.u1link
* Returns 1 if update was accepted and 0 if rejected. */
int metro_u1link(lattice const* l, fields* f, params const* p, long i, int dir) {

	double oldlink = f->u1link[i][dir];

	double linkact_old = localact_u1link(l, f, p, i, dir);

	// multiply link by a random phase (can adjust the overall number here)
	f->u1link[i][dir] += 1.0*(dran() - 0.5);

	double linkact_new = localact_u1link(l, f, p, i, dir);

	double diff = linkact_new - linkact_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > dran() )) {
		return 1;
	}
	else {
		f->u1link[i][dir] = oldlink;
		return 0;
	}
}


#if (NHIGGS > 0)
/* Update an SU(2) scalar doublet using Metropolis.
* Returns 1 if update was accepted and 0 if rejected */
int metro_doublet(lattice const* l, fields* f, params const* p, long i, int higgs_id) {

	double **phi = f->su2doublet[higgs_id];
	double oldfield[4];
	memcpy(oldfield, phi[i], SU2DB * sizeof(phi[i][0]));

	double act_old = localact_doublet(l, f, p, i, higgs_id);

	// modify the old field by random values
	for (int k=0; k<SU2DB; k++) {
		phi[i][k] += 1.0*(dran() - 0.5);
	}

	double act_new = localact_doublet(l, f, p, i, higgs_id);

	int accept;
	double diff = act_new - act_old;
	if (diff < 0) {
		accept = 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > dran() )) {
		accept = 1;
	}
	else {
		memcpy(phi[i], oldfield, SU2DB * sizeof(phi[i][0]));
		accept = 0;
	}

	return accept;
}

#endif // if (NHIGGS > 0)


/* Update a single SU(2) scalar triplet using Metropolis.
* Returns 1 if update was accepted and 0 if rejected. */
int metro_triplet(lattice const* l, fields* f, params const* p, long i) {

	double oldfield[3];
	oldfield[0] = f->su2triplet[i][0];
	oldfield[1] = f->su2triplet[i][1];
	oldfield[2] = f->su2triplet[i][2];

	double act_old = localact_triplet(l, f, p, i);

	// modify the old field by random values
	f->su2triplet[i][0] += 1.0*(dran() - 0.5);
	f->su2triplet[i][1] += 1.0*(dran() - 0.5);
	f->su2triplet[i][2] += 1.0*(dran() - 0.5);


	double act_new = localact_triplet(l, f, p, i);

	double diff = act_new - act_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > dran() )) {
		return 1;
	}
	else {
		f->su2triplet[i][0] = oldfield[0];
		f->su2triplet[i][1] = oldfield[1];
		f->su2triplet[i][2] = oldfield[2];
		return 0;
	}
}

#ifdef SINGLET

/* Metropolis update for a singlet field at site i */
int metro_singlet(lattice const* l, fields* f, params const* p, long i) {
	double oldfield = f->singlet[i][0];
	double act_old = localact_singlet(l, f, p, i);

	f->singlet[i][0] += 1.0*(dran() - 0.5);
	double act_new = localact_singlet(l, f, p, i);

	double diff = act_new - act_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > dran() )) {
		return 1;
	}
	else {
		f->singlet[i][0] = oldfield;
		return 0;
	}

}

#endif

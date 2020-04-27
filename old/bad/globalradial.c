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

/*  NB! This should NOT be used; global updates scale very badly with volume. 
* Kari ditched this update ages ago...
*
* Global radial update for the Higgs as described in hep-lat/9510020.
* The doublet Phi can be written as Phi = R U, where U is an SU(2) matrix.
* For flat potential the "slowness" of the Higgs is related to the radial modes,
* so we want to boost updating in those. Our Phi is written as 
* 	Phi = 1/sqrt(2) (phi_0 I + i phi_i sigma_i), where phi_a = f.su2doublet[i][a], 
* so clearly our R = 1/sqrt(2) * sqrt(phi_a phi_a). Thus we just update 
*		phi_a -> exp(xi) phi_a
*
* Accept/reject is first checked in root node, that then tells
* others to do the update if accepted. If rejected, no comms needed.
// TODO optimize
*/
int globalradial_doublet(params p, fields f) {
	// write Higgs action as \sum_x A * R^2 + B R^4 + hopping,
	// where R^2 = 0.5 Tr\Phi^\dagger\Phi

	double phisq = 0.0, phi4 = 0.0, hop = 0.0;
	for (long i=0; i<p.sites; i++) {
		double mod = doubletsq(f.su2doublet[i]);
		phisq += mod;
		phi4 += mod * mod;
		for (int dir=0; dir<p.dim; dir++) {
			hop += hopping_doublet_forward(f, p, i, dir);
		}
	}
	phisq = reduce_sum(phisq);
	phi4 = reduce_sum(phi4);
	hop = reduce_sum(hop);

	double a = (p.msq_phi + 2.0 * p.dim) * phisq + hop;
	double b = p.lambda_phi * phi4;

	#ifdef TRIPLET
		double phi2Sigma2 = 0.0;
		for (long i=0; i<p.sites; i++) {
			phi2Sigma2 += doubletsq(f.su2doublet[i]) * tripletsq(f.su2triplet[i]);
		}
		phi2Sigma2 = reduce_sum(phi2Sigma2);
		a += p.a2 * phi2Sigma2;
	#endif


	int accept = 0;
	double xi;
	if (p.rank == 0) {
	 	xi = drand48()*2.0 - 1.0; // [-1.0, 1.0)
		// aim for acceptance rate 60% - 70%, so the range should DECREASE with volume
		double eps = 0.01;
		xi *= eps; // [-eps, eps)


		// note that the Haarr measure changes [d\phi] -> exp(4*xi*p.vol) [d\phi],
		// because phi has 4 components and there are p.vol integral measures
		double diff = SU2DB*p.vol*xi - (a*(exp(2.0*xi)-1.0) + b*(exp(4.0*xi)-1.0));
		//printf0(p, "xi = %lf, diff = %lf \n", xi, diff);
		if (diff >= 0) {
			accept = 1;
		}
		else if ( exp(diff) > drand48() ) {
			accept = 1;
		}
		else {
			accept = 0;
		}
	}

	bcast_int(&accept);
	if (accept) {
		// OK, so go ahead and update doublet everywhere (also halos!)
		bcast_double(&xi);
		for (long i=0; i<p.sites_total; i++) {
			for (int k=0; k<SU2DB; k++) {
				f.su2doublet[i][k] *= exp(xi);
			}
		}
		return 1;
	} else {
		return 0;
	}
}

/* Local transverse update on the doublet; leaves phisq invariant.
* DOES NOT WORK?!?!
*/
int metro_doublet_transverse(fields f, params p, long i) {

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

	double norm_old = doubletsq(oldfield);
	double norm_new = doubletsq(f.su2doublet[i]);

	// keep phisq unchanged
	for (int k=0; k<SU2DB; k++) {
		f.su2doublet[i][k] *= norm_old / norm_new;
	}

	double act_new = localact_doublet(f, p, i);

	double diff = act_new - act_old;
	if (diff < 0) {
		return 1;
	}
	else if (diff > 0 && ( exp(-(diff)) > drand48() )) {
		return 1;
	}
	else {
		f.su2doublet[i][0] = oldfield[0];
		f.su2doublet[i][1] = oldfield[1];
		f.su2doublet[i][2] = oldfield[2];
		f.su2doublet[i][3] = oldfield[3];
		return 0;
	}

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

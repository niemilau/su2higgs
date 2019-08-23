/** @file update.c
*
* Routines for updating and sweeping over the lattice
*
* TODO
*
*/

#include "su2.h"


/* Sweep over the lattice in a checkerboard layout and update half of the links.
* Gauge links are updated only in direction specified by dir.
* Parity of a site is assumed to be stored in params.
*
* Future optimization: loop only over sites with the correct parity?
* I have tried storing odd / even sites in predetermined arrays and using those
* for the loop, but it did not result in consistent improvement in computation time in short test runs...
*/
void checkerboard_sweep_su2link(fields f, params p, counters* c, char parity, int dir) {
	// EVEN sites come before ODD 
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p.evensites;
	} else {
		offset = p.evensites; max = p.sites;
	}

	for (long i=offset; i<max; i++) {
		if (p.algorithm_su2link == HEATBATH) {
			c->accepted_su2link += heatbath_su2link(f, p, i, dir);
			c->total_su2link++;
		} else if (p.algorithm_su2link == METROPOLIS) {
			c->accepted_su2link += metro_su2link(f, p, i, dir);
			c->total_su2link++;
		}
	}
	
}


/* Sweep over the lattice in a checkerboard layout and update half of the doublets.
* Parity of a site is assumed to be stored in params.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
*/
void checkerboard_sweep_su2doublet(fields f, params p, weight* w, counters* c, char parity, char metro) {

	double** oldfield;
	double muca_param_old, muca_param_new;

	if (p.multicanonical) {
		// multicanonical order parameter BEFORE updates
		muca_param_old = calc_orderparam(p, f);

		// Store the old field values, in case the muca step is rejected
		long paritysites;
		if (parity == EVEN) {
			paritysites = p.evensites;
		} else {
			paritysites = p.oddsites;
		}

		oldfield = make_field(paritysites, SU2DB);
		long k = 0;
		for (long i=0; i<p.sites; i++) {
			if (p.parity[i] == parity) {
				for (int d=0; d<SU2DB; d++) {
					oldfield[k][d] = f.su2doublet[i][d];
				}
				k++;
			}
		}
	}
	
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p.evensites;
	} else {
		offset = p.evensites; max = p.sites;
	}
	
	for (long i=offset; i<max; i++) {
		if (p.algorithm_su2doublet == OVERRELAX && (metro != 1)) {
			c->acc_overrelax_doublet += overrelax_doublet(f, p, i);
			c->total_overrelax_doublet++;

		} else if (p.algorithm_su2doublet == METROPOLIS || (metro == 1)) {
			c->accepted_doublet += metro_doublet(f, p, i);
			c->total_doublet++;
		}
	}

	// sweep done, now recalculate the muca order parameter and acc/rej
	if (p.multicanonical) {

		muca_param_new = calc_orderparam(p, f);
		int accept = multicanonical_acceptance(p, w, muca_param_old, muca_param_new);
		// multicanonical_acceptance() also does weight update if accepted

		if (!accept) {
			// multicanonical step was rejected, undo field changes
			long k = 0;
			for (long i=0; i<p.sites; i++) {
				if (p.parity[i] == parity) {
					for (int d=0; d<SU2DB; d++) {
						f.su2doublet[i][d] = oldfield[k][d];
					}
					k++;
				}
			}

		}

		c->accepted_muca += accept;
		c->total_muca++;
		free_field(oldfield);
	}

}


/* Sweep over the lattice in a checkerboard layout and update half of the triplet.
* Parity of a site is assumed to be stored in params.
*/
void checkerboard_sweep_su2triplet(fields f, params p, counters* c, char parity, char metro) {
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p.evensites;
	} else {
		offset = p.evensites; max = p.sites;
	}
	
	for (long i=offset; i<max; i++) {
		if (p.algorithm_su2triplet == OVERRELAX && (metro != 1)) {
			c->acc_overrelax_triplet += overrelax_triplet(f, p, i);
			c->total_overrelax_triplet++;
		} else if (p.algorithm_su2triplet == METROPOLIS || (metro == 1)) {
			c->accepted_triplet += metro_triplet(f, p, i);
			c->total_triplet++;
		}
	
	}
}


/* Full update on all sites + halo communication.
* Following, hep-lat/9804019, we first update the gauge links and then scalars.
* Each link direction is updated separately, I.E. first dir1 with even and odd,
* then dir2 etc. This is necessary because the links in "positive" directions are not independent
* because of the Wilson staple, which for U_1(x) depends on U_2(x-i+j), for example.
*/
void update_lattice(fields f, params p, comlist_struct* comlist, weight* w, counters* c, char metro) {

	// update each link direction separately.
	// I.E. first dir1 with even and odd, then dir2 etc
	// This is necessary because the links in "positive" directions are not independent
	// because of the Wilson staple, which for U_1(x) depends on U_2(x-i+j), for example
	for (int dir=0; dir<p.dim; dir++) {
		checkerboard_sweep_su2link(f, p, c, EVEN, dir);
		c->comms_time += update_gaugehalo(comlist, EVEN, f.su2link, SU2LINK, dir);
		checkerboard_sweep_su2link(f, p, c, ODD, dir);
		c->comms_time += update_gaugehalo(comlist, ODD, f.su2link, SU2LINK, dir);
	}

	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
		checkerboard_sweep_su2doublet(f, p, w, c, EVEN, metro);
		c->comms_time += update_halo(comlist, EVEN, f.su2doublet, SU2DB);

		checkerboard_sweep_su2doublet(f, p, w, c, ODD, metro);
		c->comms_time += update_halo(comlist, ODD, f.su2doublet, SU2DB);
	}

	#endif
	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
		checkerboard_sweep_su2triplet(f, p, c, EVEN, metro);
		c->comms_time += update_halo(comlist, EVEN, f.su2triplet, SU2TRIP);

		checkerboard_sweep_su2triplet(f, p, c, ODD, metro);
		c->comms_time += update_halo(comlist, ODD, f.su2triplet, SU2TRIP);
	}
	#endif


}

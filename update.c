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
void checkerboard_sweep_su2doublet(fields f, params p, counters* c, char parity, char metro) {

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
void update_lattice(fields* f, params p, comlist_struct* comlist, counters* c, char metro) {

	// update each link direction separately.
	// I.E. first dir1 with even and odd, then dir2 etc
	// This is necessary because the links in "positive" directions are not independent
	// because of the Wilson staple, which for U_1(x) depends on U_2(x-i+j), for example
	for (int dir=0; dir<p.dim; dir++) {
		checkerboard_sweep_su2link(*f, p, c, EVEN, dir);
		c->comms_time += update_gaugehalo(comlist, EVEN, f->su2link, SU2LINK, dir);
		checkerboard_sweep_su2link(*f, p, c, ODD, dir);
		c->comms_time += update_gaugehalo(comlist, ODD, f->su2link, SU2LINK, dir);
	}

	// EVEN sweeps for scalars
	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
		checkerboard_sweep_su2doublet(*f, p, c, EVEN, metro);
		c->comms_time += update_halo(comlist, EVEN, f->su2doublet, SU2DB);
	}
	#endif

	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
		checkerboard_sweep_su2triplet(*f, p, c, EVEN, metro);
		c->comms_time += update_halo(comlist, EVEN, f->su2triplet, SU2TRIP);
	}
	#endif

	// ODD sweeps for scalars
	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
		checkerboard_sweep_su2doublet(*f, p, c, ODD, metro);
		c->comms_time += update_halo(comlist, ODD, f->su2doublet, SU2DB);
	}
	#endif
	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
		checkerboard_sweep_su2triplet(*f, p, c, ODD, metro);
		c->comms_time += update_halo(comlist, ODD, f->su2triplet, SU2TRIP);
	}
	#endif

}

/* Same as update_lattice(), but takes multicanonical weight
* into account. Note that multicanonical_acceptance() also does
* weight update if accepted. Here the timing of the call to
* multicanonical_acceptance() is crucial: we need to call it
* only after ALL fields contributing to the order parameter
* have been updated. Note that some halo updates may need to be undone
* as well; this is implemented in store_muca_fields() and reset_muca_fields(),
* and assume that the Higgs is updated LAST.
*/
void update_lattice_muca(fields* f, params p, comlist_struct* comlist, weight* w, counters* c, char metro) {

	int accept;
	double muca_param_old = w->param_value[EVEN] + w->param_value[ODD];
	store_muca_fields(p, f, w);

	for (int dir=0; dir<p.dim; dir++) {
		checkerboard_sweep_su2link(*f, p, c, EVEN, dir);
		c->comms_time += update_gaugehalo(comlist, EVEN, f->su2link, SU2LINK, dir);
		checkerboard_sweep_su2link(*f, p, c, ODD, dir);
		c->comms_time += update_gaugehalo(comlist, ODD, f->su2link, SU2LINK, dir);
	}

	// parity loop for scalars. EVEN = 0, ODD = 1; defined in su2.h
	for (char par=0; par<=1; par++) {

		#ifdef TRIPLET
		for (int k=0; k<p.update_su2triplet; k++) {
			checkerboard_sweep_su2triplet(*f, p, c, par, metro);

			// multicanonical step if the order parameter ONLY depends on the adjoint
			if (w->orderparam == SIGMASQ) {
				double muca_param_new = calc_orderparam(p, *f, w, par); // this also updates w->param_value[par]
				accept = multicanonical_acceptance(p, w, muca_param_old, muca_param_new);

				if (!accept) {
					// rejected, undo field changes and w->param_value
					reset_muca_fields(p, f, w, par);
					w->param_value[par] = muca_param_old - w->param_value[otherparity(par)];
				}
			}

			c->comms_time += update_halo(comlist, par, f->su2triplet, SU2TRIP);
		}
		#endif

		#ifdef HIGGS
		for (int k=0; k<p.update_su2doublet; k++) {
			checkerboard_sweep_su2doublet(*f, p, c, par, metro);

			// multicanonical step if the order parameter depends on the Higgs
			if (w->orderparam == PHISQ || w->orderparam == PHI2SIGMA2) {
				double muca_param_new = calc_orderparam(p, *f, w, par); // this also updates w->param_value[par]
				accept = multicanonical_acceptance(p, w, muca_param_old, muca_param_new);

				if (!accept) {
					// rejected, undo field changes and w->param_value
					reset_muca_fields(p, f, w, par);
					w->param_value[par] = muca_param_old - w->param_value[otherparity(par)];
				}
			}

			c->comms_time += update_halo(comlist, par, f->su2doublet, SU2DB);
		}
		#endif
		
		// update muca_param_old for the next parity
		muca_param_old = w->param_value[EVEN] + w->param_value[ODD];
		c->accepted_muca += accept;
		c->total_muca++;
	} // end parity loop

}

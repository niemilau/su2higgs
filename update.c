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
*/
void checkerboard_sweep_su2link(fields* f, params const* p, counters* c, char parity, int dir) {
	// EVEN sites come before ODD
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	// check algorithm before looping: this should provide some optimization.
	if (p->algorithm_su2link == HEATBATH) {

		for (long i=offset; i<max; i++) {
			c->accepted_su2link += heatbath_su2link(f, p, i, dir);
			c->total_su2link++;
		}
	} else if (p->algorithm_su2link == METROPOLIS) {

		for (long i=offset; i<max; i++) {
			c->accepted_su2link += metro_su2link(f, p, i, dir);
			c->total_su2link++;
		}
	}

}

/* Same as checkerboard_sweep_su2link(), but for U(1) links instead.
*/
void checkerboard_sweep_u1link(fields* f, params const* p, counters* c, char parity, int dir) {
	// EVEN sites come before ODD
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	for (long i=offset; i<max; i++) {
		if (p->algorithm_u1link == HEATBATH) {
			//c->accepted_u1link += heatbath_su2link(f, p, i, dir);
		} else if (p->algorithm_u1link == METROPOLIS) {
			c->accepted_u1link += metro_u1link(f, p, i, dir);
		}
		c->total_u1link++;
	}

}


/* Sweep over the lattice in a checkerboard layout and update half of the doublets.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
* Return value is 1 if the entire sweep was accepted by multicanonical, 0 otherwise.
*/
int checkerboard_sweep_su2doublet(fields* f, params const* p, counters* c, weight* w, char parity, char metro) {

	int accept = 1;
	double muca_param_old = 0;
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	// multicanonical preparations if the order parameter depends on the doublet
	int do_muca = 0;
	if (p->multicanonical) {
		if (w->orderparam == PHISQ || w->orderparam == PHI2MINUSSIGMA2 || w->orderparam == PHI2SIGMA2) {
			muca_param_old = w->param_value[EVEN] + w->param_value[ODD];
			store_muca_fields(p, f, w);
			do_muca = 1;
		}
	}

	// then the update sweep
	if (p->algorithm_su2doublet == OVERRELAX && (metro == 0)) {
		for (long i=offset; i<max; i++) {

			c->acc_overrelax_doublet += overrelax_doublet(f, p, i);
			c->total_overrelax_doublet++;
		}

	} else if (p->algorithm_su2doublet == METROPOLIS || (metro != 0)) {
		for (long i=offset; i<max; i++) {
			c->accepted_doublet += metro_doublet(f, p, i);
			c->total_doublet++;
		}
	}

	// now the global multicanonical step
	if (do_muca) {
		double muca_param_new = calc_orderparam(p, f, w, parity); // this also updates w->param_value[par]
		accept = multicanonical_acceptance(p, w, muca_param_old, muca_param_new);

		if (!accept) {
			// rejected, undo field changes and w->param_value
			reset_muca_fields(p, f, w, parity);
			w->param_value[parity] = muca_param_old - w->param_value[otherparity(parity)];
		}

		c->accepted_muca += accept;
		c->total_muca++;
	}

	return accept;

}


/* Sweep over the lattice in a checkerboard layout and update half of the triplets.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
* Return value is 1 if the entire sweep was accepted by multicanonical, 0 otherwise.
*/
int checkerboard_sweep_su2triplet(fields* f, params const* p, counters* c, weight* w, char parity, char metro) {

	int accept = 1;
	double muca_param_old = 0;
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	// multicanonical preparations if the order parameter depends on the triplet
	int do_muca = 0;
	if (p->multicanonical) {
		if (w->orderparam == SIGMASQ || w->orderparam == PHI2MINUSSIGMA2 || w->orderparam == PHI2SIGMA2) {
			muca_param_old = w->param_value[EVEN] + w->param_value[ODD];
			store_muca_fields(p, f, w);
			do_muca = 1;
		}
	}

	// then the update sweep
	if (p->algorithm_su2triplet == OVERRELAX && (metro == 0)) {
		for (long i=offset; i<max; i++) {
			c->acc_overrelax_triplet += overrelax_triplet(f, p, i);
			c->total_overrelax_triplet++;
		}
	} else if (p->algorithm_su2triplet == METROPOLIS || (metro != 0)) {
		for (long i=offset; i<max; i++) {
			c->accepted_triplet += metro_triplet(f, p, i);
			c->total_triplet++;
		}
	}

	// now the global multicanonical step
	if (do_muca) {
		double muca_param_new = calc_orderparam(p, f, w, parity); // this also updates w->param_value[par]
		accept = multicanonical_acceptance(p, w, muca_param_old, muca_param_new);

		if (!accept) {
			// rejected, undo field changes and w->param_value
			reset_muca_fields(p, f, w, parity);
			w->param_value[parity] = muca_param_old - w->param_value[otherparity(parity)];
		}

		c->accepted_muca += accept;
		c->total_muca++;
	}

	return accept;

}


/* Full update on all sites + halo communication.
* Following, hep-lat/9804019, we first update the gauge links and then scalars.
* Each link direction is updated separately, I.E. first dir1 with even and odd,
* then dir2 etc. This is necessary because the links in "positive" directions are not independent
* because of the Wilson staple, which for U_1(x) depends on U_2(x-i+j), for example.
*/
void update_lattice(fields* f, params const* p, comlist_struct* comlist, counters* c, weight* w) {

	int accept;

	// EVEN and ODD sweeps for gauge links
	for (int k=0; k<p->update_links; k++) {
		for (int dir=0; dir<p->dim; dir++) {
			checkerboard_sweep_su2link(f, p, c, EVEN, dir);
			c->comms_time += update_gaugehalo(comlist, EVEN, f->su2link, SU2LINK, dir);
			checkerboard_sweep_su2link(f, p, c, ODD, dir);
			c->comms_time += update_gaugehalo(comlist, ODD, f->su2link, SU2LINK, dir);
		}
	}
	#ifdef U1
	// here I use update_halo() instead of update_gaugehalo(), so halo is
	// actually updated for all directions after updating just one direction.
	// this is of course OK, but could be optimized.
	for (int k=0; k<p->update_links; k++) {
		for (int dir=0; dir<p->dim; dir++) {
			checkerboard_sweep_u1link(f, p, c, EVEN, dir);
			c->comms_time += update_halo(comlist, EVEN, f->u1link, p->dim);
			checkerboard_sweep_u1link(f, p, c, ODD, dir);
			c->comms_time += update_halo(comlist, ODD, f->u1link, p->dim);
		}
	}
	#endif


	/* scalar updates. We update all scalar fields p->scalar_sweeps times
	* per iteration. Ordering is such that each field gets a full sweep before
	* moving on to the next field. Additionally, each individual field
	* can be sweeped over k times per "scalar sweep", where k is specified in
	* config, e.g. p->update_su2doublet = 2 will update Higgs 2 times for each
	* scalar sweep. The point here is that p->scalar_sweeps controls how the scalar
	* evolve on an equal footing, while the extra sweeps for different fields are
	* reserved for special use (such as keeping the Higgs constant).
	* parity loop uses EVEN = 0, ODD = 1; defined in su2.h */

	// how many sweeps before forcing metropolis, for ergodicity
	int metro_interval = 5;
	int metro; // no effect if p.algorithm is set to METROPOLIS


	for (int s=0; s<p->scalar_sweeps; s++) {

		#ifdef HIGGS
		for (int k=0; k<p->update_su2doublet; k++) {

			if (c->higgs_sweeps >= metro_interval-1) {
				metro = 1;
				c->higgs_sweeps = 0;
			} else {
				metro = 0;
				c->higgs_sweeps++;
			}

			for (char par=0; par<=1; par++) {
				accept = checkerboard_sweep_su2doublet(f, p, c, w, par, metro);

				// if the sweep was rejected, no need to sync halos
				if (accept) {
					c->comms_time += update_halo(comlist, par, f->su2doublet, SU2DB);
				}

			}
		}
		#endif

		#ifdef TRIPLET
		for (int k=0; k<p->update_su2triplet; k++) {

			if (c->triplet_sweeps >= metro_interval-1) {
				metro = 1;
				c->triplet_sweeps = 0;
			} else {
				metro = 0;
				c->triplet_sweeps++;
			}

			for (char par=0; par<=1; par++) {
				accept = checkerboard_sweep_su2triplet(f, p, c, w, par, metro);
				// if the sweep was rejected, no need to sync halos
				if (accept) {
					c->comms_time += update_halo(comlist, par, f->su2triplet, SU2TRIP);
				}

			}
		}
		#endif
	}
}


/* Routine sync_halos()
* This is to be called before starting the simulation, to ensure that
* all halo fields are in sync.
*/
void sync_halos(fields* f, params const* p, comlist_struct* comlist) {

	for (int parity=0; parity<=1; parity++) {

		for (int dir=0; dir<p->dim; dir++) {
			update_gaugehalo(comlist, parity, f->su2link, SU2LINK, dir);
			#ifdef U1
				update_halo(comlist, parity, f->u1link, p->dim);
			#endif
		}

		#ifdef HIGGS
			update_halo(comlist, parity, f->su2doublet, SU2DB);
		#endif

		#ifdef TRIPLET
			update_halo(comlist, parity, f->su2triplet, SU2TRIP);
		#endif

	}
}

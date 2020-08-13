/** @file update.c
*
* Routines for updating and sweeping over the lattice
*
* TODO
*
*/

#include "su2.h"


/* Sweep over the lattice in a checkerboard layout and update half of the links.
* Gauge links are updated only in direction specified by dir. */
void checkerboard_sweep_su2link(lattice const* l, fields* f, params const* p, counters* c, int parity, int dir) {
	// EVEN sites come before ODD
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = l->evensites;
	} else {
		offset = l->evensites; max = l->sites;
	}

	for (long i=offset; i<max; i++) {
		if (p->algorithm_su2link == HEATBATH) {
			c->accepted_su2link += heatbath_su2link(l, f, p, i, dir);
		} else if (p->algorithm_su2link == METROPOLIS) {
			c->accepted_su2link += metro_su2link(l, f, p, i, dir);
		}
		c->total_su2link++;
	}

}

/* Same as checkerboard_sweep_su2link(), but for U(1) links instead. */
void checkerboard_sweep_u1link(lattice const* l, fields* f, params const* p, counters* c, int parity, int dir) {
	// EVEN sites come before ODD
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = l->evensites;
	} else {
		offset = l->evensites; max = l->sites;
	}

	for (long i=offset; i<max; i++) {
		if (p->algorithm_u1link == HEATBATH) {
			//c->accepted_u1link += heatbath_su2link(f, p, i, dir);
		} else if (p->algorithm_u1link == METROPOLIS) {
			c->accepted_u1link += metro_u1link(l, f, p, i, dir);
		}
		c->total_u1link++;
	}

}


#if (NHIGGS > 0)
/* Sweep over the lattice in a checkerboard layout and update half of the doublets.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
* Return value is 0 if ALL local updates were rejected by multicanonical, nonzero otherwise */
int checkerboard_sweep_su2doublet(lattice const* l, fields* f, params const* p, counters* c,
			weight* w, int parity, int metro, int higgs_id) {

	int accept = 1;
	long muca_interval = l->sites_total; // initialize to some large value to avoid bugs
	long muca_count = 0;
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = l->evensites;
	} else {
		offset = l->evensites; max = l->sites;
	}

	// check if multicanonical is to be used
	int do_muca = 0;
	if (w->do_acceptance) {
		if (w->orderparam == PHISQ || w->orderparam == PHI2MINUSSIGMA2 || w->orderparam == PHI2SIGMA2) {
			store_muca_fields(l, f, w);
			do_muca = 1;
			muca_interval = (max - offset) / w->checks_per_sweep; // takes floor if not integer
			if (muca_interval <= 0) muca_interval = 1;
			accept = 0; // the sweep may be rejected by multicanonical
		}
	}

	// then the update sweep, doing a global muca acc/rej every muca_interval sites
	for (long i=offset; i<max; i++) {

		if (p->algorithm_su2doublet == OVERRELAX && (metro == 0)) {

			#if (NHIGGS == 2)
				c->acc_overrelax_doublet[higgs_id] += overrelax_higgs2(l, f, p, i, higgs_id);
			#else
				c->acc_overrelax_doublet[higgs_id] += overrelax_doublet(l, f, p, i);
			#endif
			c->total_overrelax_doublet[higgs_id]++;

		} else if (p->algorithm_su2doublet == METROPOLIS || (metro != 0)) {
			c->accepted_doublet[higgs_id] += metro_doublet(l, f, p, i, higgs_id);
			c->total_doublet[higgs_id]++;
		}

		if (do_muca) {
			muca_count++;
			if (muca_count % muca_interval == 0) {
				// do the global muca acc/rej step, and take new backups unless the sweep is finished
				int make_backups = (i < max-1);
				accept += muca_check(l, f, p, c, w, parity, make_backups);
				/* if rejected, the field changes have now been undone. */
			} // end muca check
		} // end do muca

	} // end site loop

	// if the muca interval did not add up, do a final check here without taking new backups
	if (do_muca && (max - offset) % w->checks_per_sweep != 0) {
		accept += muca_check(l, f, p, c, w, parity, 0);
	}

	return accept; // return is nonzero if at least one muca check was accepted
}

#endif // if (NHIGGS > 0)

/* */
int muca_check(lattice const* l, fields* f, params const* p, counters* c, weight* w, int parity, int make_backups) {

	// local field updates do not recalculate the muca order parameter, so the old value is still in w->param_value:
	double orderparam_old = w->param_value[EVEN] + w->param_value[ODD];
	// recalculate with the new fields; this also updates w->param_value[par]
	double orderparam_new = calc_orderparam(l, f, p, w, parity);

	int accept = multicanonical_acceptance(l, w, orderparam_old, orderparam_new);
	if (!accept) {
		// rejected, undo field changes and w->param_value
		reset_muca_fields(l, f, w, parity);
		w->param_value[parity] = orderparam_old - w->param_value[otherparity(parity)];
	} else if (make_backups) {
		// accepted, make new field backups
		store_muca_fields(l, f, w);
	}
	c->accepted_muca += accept;
	c->total_muca++;

	return accept;
}


/* Sweep over the lattice in a checkerboard layout and update half of the triplets.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
* Return value is 0 if ALL local updates were rejected by multicanonical, nonzero otherwise */
int checkerboard_sweep_su2triplet(lattice const* l, fields* f, params const* p, counters* c, weight* w, int parity, int metro) {

	int accept = 1;
	double muca_param_old = 0;
	long offset, max;
	if (parity == EVEN) {
		offset = 0; max = l->evensites;
	} else {
		offset = l->evensites; max = l->sites;
	}

	// multicanonical preparations if the order parameter depends on the triplet
	int do_muca = 0;
	if (w->do_acceptance) {
		if (w->orderparam == SIGMASQ || w->orderparam == PHI2MINUSSIGMA2 || w->orderparam == PHI2SIGMA2) {
			muca_param_old = w->param_value[EVEN] + w->param_value[ODD];
			store_muca_fields(l, f, w);
			do_muca = 1;
		}
	}

	// then the update sweep
	for (long i=offset; i<max; i++) {
		if (p->algorithm_su2triplet == OVERRELAX && (metro == 0)) {
			c->acc_overrelax_triplet += overrelax_triplet(l, f, p, i);
			c->total_overrelax_triplet++;
		} else if (p->algorithm_su2triplet == METROPOLIS || (metro != 0)) {
			c->accepted_triplet += metro_triplet(l, f, p, i);
			c->total_triplet++;
		}
	}

	// now the global multicanonical step
	if (do_muca) {
		double muca_param_new = calc_orderparam(l, f, p, w, parity); // this also updates w->param_value[par]
		accept = multicanonical_acceptance(l, w, muca_param_old, muca_param_new);

		if (!accept) {
			// rejected, undo field changes and w->param_value
			reset_muca_fields(l, f, w, parity);
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
* However for realtime heatbath, we also allow for a random ordering of sweeps.
*/
void update_lattice(lattice* l, fields* f, params const* p, counters* c, weight* w) {

	int accept;

	/* Sweeps for gauge links. Arrays par_a (parity) and dir_a (directions)
	* specify the order of updates. These have length 2*dim. Defaults are (with dim=3)
	* par_a = [0,1,0,1,0,1], dirs_a = [0,0,1,1,2,2], i.e.
	* first EVEN and ODD links in direction 0, then EVEN/ODD in dir 1 etc. */
	for (int k=0; k<p->update_links; k++) { // repeat the whole process p->update_links times

		const int NV = 2*l->dim;
		int par_a[NV], dir_a[NV]; // parities, directions
		int dir = 0;
		// default ordering, non-random
		for (int j=0; j<NV; j+=2) {
			par_a[j] = EVEN;
			par_a[j+1] = ODD;
			dir_a[j] = dir;
			dir_a[j+1] = dir;
			dir++;
		}

		if (p->random_sweeps) {
			// randomize order in root node and broadcast to others
			if (!l->rank) {
				shuffle(par_a, NV);
				shuffle(dir_a, NV);
			}
			bcast_int_array(par_a, NV, l->comm);
			bcast_int_array(dir_a, NV, l->comm);
		}

		// now update in the specified order
		for (int j=0; j<NV; j++) {
			int dir = dir_a[j];
			int par = par_a[j];
			checkerboard_sweep_su2link(l, f, p, c, par, dir);
			update_gaugehalo(l, par, f->su2link, SU2LINK, dir);

			#ifdef U1
				checkerboard_sweep_u1link(l, f, p, c, par, dir);
				// here I use update_halo() instead of update_gaugehalo(), so halo is
				// actually updated for all directions after updating just one direction.
				// Could be optimized. (is it OK to update U1 together with SU2 like here?)
				update_halo(l, par, f->u1link, l->dim);
			#endif
		}
	} // gauge links done


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

		// Parity ordering, default par_a[EVEN, ODD].
		int par_a[2];
		if (p->random_sweeps) {
			par_a[0] = (drand48() < 0.5) ? EVEN : ODD;
			par_a[1] = otherparity(par_a[0]);
			bcast_int_array(par_a, 2, l->comm);
		} else {
			par_a[0] = EVEN; par_a[1] = ODD;
		}

		#if (NHIGGS > 0)
		for (int k=0; k<p->update_su2doublet; k++) {

			if (c->higgs_sweeps >= metro_interval-1) {
				metro = 1;
				c->higgs_sweeps = 0;
			} else {
				metro = 0;
				c->higgs_sweeps++;
			}

			for (int j=0; j<=1; j++) {
				int par = par_a[j];

				for (int db=0; db<NHIGGS; db++) {
					accept = checkerboard_sweep_su2doublet(l, f, p, c, w, par, metro, db);
					if (accept) update_halo(l, par, f->su2doublet[db], SU2DB);
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

			for (int j=0; j<=1; j++) {
				int par = par_a[j];
				accept = checkerboard_sweep_su2triplet(l, f, p, c, w, par, metro);
				// if the sweep was rejected, no need to sync halos
				if (accept) update_halo(l, par, f->su2triplet, SU2TRIP);

			}
		}
		#endif
	}
}


/* Routine sync_halos(): just syncs all halo fields. */
void sync_halos(lattice* l, fields* f) {

	for (int parity=0; parity<=1; parity++) {

		for (int dir=0; dir<l->dim; dir++) {
			update_gaugehalo(l, parity, f->su2link, SU2LINK, dir);
			#ifdef U1
				update_halo(l, parity, f->u1link, l->dim);
			#endif
		}

		#if (NHIGGS > 0)
			for (int db=0; db < NHIGGS; db++) update_halo(l, parity, f->su2doublet[db], SU2DB);
		#endif

		#ifdef TRIPLET
			update_halo(l, parity, f->su2triplet, SU2TRIP);
		#endif

	}
}

// Shuffle an array
void shuffle(int *arr, int len) {
    int i, temp;
    for(i = len-1; i > 0; i--) {
        int j = rand() % (i+1);
				temp = arr[i];
				arr[i] = arr[j];
				arr[j] = temp;
    }
}

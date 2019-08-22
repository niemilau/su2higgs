/** @file update.c
*
* Routines for updating and sweeping over the lattice
*
* TODO
*
*/

#include "su2.h"


/* Sweep over the lattice in a checkerboard layout and update half of the sites.
* Gauge links are updated only in direction specified by dir.
* Parity of a site is assumed to be stored in params.
*/
// TODO figure out if this is even allowed
// todo figure out possible optimizations
// NOT USED ATM
void checkerboard_sweep(fields f, params p, counters* c, int dir, char parity) {

	// figure out a faster loop here? only loop over sites with the correct parity
	for (ulong i=0; i<p.vol; i++) {
		if (p.parity[i] == parity) {

			// SU(2) links
			if (p.algorithm_su2link == HEATBATH) {
				heatbath_su2link(f, p, i, dir);
			} else if (p.algorithm_su2link == METROPOLIS) {
				c->accepted_su2link += metro_su2link(f, p, i, dir);
				c->total_su2link++;
			}
			#ifdef HIGGS
			// SU(2) doublets
			for (int k=0; k<p.update_su2doublet; k++) {
				if (p.algorithm_su2doublet == OVERRELAX) {
					c->acc_overrelax_doublet += overrelax_doublet(f, p, i);
				} else if (p.algorithm_su2doublet == METROPOLIS) {
					c->accepted_doublet += metro_doublet(f, p, i);
				}
				c->total_doublet++;
			}
			#endif
		}
	}

}


/* Sweep over the lattice in a checkerboard layout and update half of the links.
* Gauge links are updated only in direction specified by dir.
* Parity of a site is assumed to be stored in params.
*/
void checkerboard_sweep_su2link(fields f, params p, counters* c, char parity, int dir) {

	// figure out a faster loop here? only loop over sites with the correct parity
	for (ulong i=0; i<p.vol; i++) {
		if (p.parity[i] == parity) {

			if (p.algorithm_su2link == HEATBATH) {
				heatbath_su2link(f, p, i, dir);
			} else if (p.algorithm_su2link == METROPOLIS) {
				c->accepted_su2link += metro_su2link(f, p, i, dir);
				c->total_su2link++;
			}
		}
	}
}

/* Sweep over the lattice in a checkerboard layout and update half of the doublets.
* Parity of a site is assumed to be stored in params.
* Last argument metro is 1 if we force a metropolis update and 0 otherwise.
*/
void checkerboard_sweep_su2doublet(fields f, params p, counters* c, char parity, char metro) {

	for (ulong i=0; i<p.vol; i++) {
		if (p.parity[i] == parity) {

			// SU(2) links
			if (p.algorithm_su2doublet == OVERRELAX && (metro != 1)) {
					c->acc_overrelax_doublet += overrelax_doublet(f, p, i);
					c->total_overrelax_doublet++;
				} else if (p.algorithm_su2doublet == METROPOLIS || (metro == 1)) {
					c->accepted_doublet += metro_doublet(f, p, i);
					c->total_doublet++;
				}
		}
	}
}

/* Sweep over the lattice in a checkerboard layout and update half of the triplet.
* Parity of a site is assumed to be stored in params.
*/
void checkerboard_sweep_su2triplet(fields f, params p, counters* c, char parity, char metro) {
	for (ulong i=0; i<p.vol; i++) {
		if (p.parity[i] == parity) {
			if (p.algorithm_su2triplet == OVERRELAX && (metro != 1)) {
				c->acc_overrelax_triplet += overrelax_triplet(f, p, i);
				c->total_overrelax_triplet++;
			} else if (p.algorithm_su2triplet == METROPOLIS || (metro == 1)) {
				c->accepted_triplet += metro_triplet(f, p, i);
				c->total_triplet++;
			}
		}
	}
}


/* Full update on all sites
*/
void update_lattice(fields f, params p, counters* c, char metro) {
	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
		checkerboard_sweep_su2doublet(f, p, c, 0, metro);
		checkerboard_sweep_su2doublet(f, p, c, 1, metro);
	}
	#endif
	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
		checkerboard_sweep_su2triplet(f, p, c, 0, metro);
		checkerboard_sweep_su2triplet(f, p, c, 1, metro);
	}
	#endif
	// update each link direction separately.
	// I.E. first dir1 with even and odd, then dir2 etc
	// This is needed because the links in "positive" directions are not independent
	for (int dir=0; dir<p.dim; dir++) {
		checkerboard_sweep_su2link(f, p, c, 0, dir);
		checkerboard_sweep_su2link(f, p, c, 1, dir);
	}

}

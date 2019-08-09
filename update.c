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

	for (long i=0; i<p.sites; i++) {
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
void checkerboard_sweep_su2doublet(fields f, params p, counters* c, char parity, char metro, char transverse) {

	for (long i=0; i<p.sites; i++) {
		if (p.parity[i] == parity) {
			if (p.algorithm_su2doublet == OVERRELAX && (metro != 1)) {
				c->acc_overrelax_doublet += overrelax_doublet(f, p, i);
				c->total_overrelax_doublet++;
			} else if (p.algorithm_su2doublet == METROPOLIS || (metro == 1)) {
				c->accepted_doublet += metro_doublet(f, p, i, transverse);
				c->total_doublet++;
			}
		}
	}
}

/* Sweep over the lattice in a checkerboard layout and update half of the triplet.
* Parity of a site is assumed to be stored in params.
*/
void checkerboard_sweep_su2triplet(fields f, params p, counters* c, char parity, char metro, char transverse) {

	for (long i=0; i<p.sites; i++) {
		if (p.parity[i] == parity) {
			if (p.algorithm_su2triplet == OVERRELAX && (metro != 1)) {
				c->acc_overrelax_triplet += overrelax_triplet(f, p, i);
				c->total_overrelax_triplet++;
			} else if (p.algorithm_su2triplet == METROPOLIS || (metro == 1)) {
				c->accepted_triplet += metro_triplet(f, p, i, transverse);
				c->total_triplet++;
			}
		}
	}
}


/* Full update on all sites + halo communication
*/
void update_lattice(fields f, params p, comlist_struct comlist, counters* c, char metro) {

	char trans = 0;
	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
		checkerboard_sweep_su2doublet(f, p, c, EVEN, metro, trans);
		c->comms_time += update_halo(&comlist, EVEN, f.su2doublet, SU2DB);
		checkerboard_sweep_su2doublet(f, p, c, ODD, metro, trans);
		c->comms_time += update_halo(&comlist, ODD, f.su2doublet, SU2DB);
	}
	#endif
	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
		checkerboard_sweep_su2triplet(f, p, c, EVEN, metro, trans);
		c->comms_time += update_halo(&comlist, EVEN, f.su2triplet, SU2TRIP);
		checkerboard_sweep_su2triplet(f, p, c, ODD, metro, trans);
		c->comms_time += update_halo(&comlist, ODD, f.su2triplet, SU2TRIP);
	}
	#endif
	// update each link direction separately.
	// I.E. first dir1 with even and odd, then dir2 etc
	// This is necessary because the links in "positive" directions are not independent
	// because of the Wilson staple, which for U_1(x) depends on U_2(x-i+j), for example
	for (int dir=0; dir<p.dim; dir++) {
		checkerboard_sweep_su2link(f, p, c, EVEN, dir);
		c->comms_time += update_gaugehalo(&comlist, EVEN, f.su2link, SU2LINK, dir);
		checkerboard_sweep_su2link(f, p, c, ODD, dir);
		c->comms_time += update_gaugehalo(&comlist, ODD, f.su2link, SU2LINK, dir);
	}

	// finish with a transverse (metro) update on the scalars, TEST
	trans = 1;
	metro = 1;
	#ifdef HIGGS
	for (int k=0; k<p.update_su2doublet; k++) {
	 checkerboard_sweep_su2doublet(f, p, c, EVEN, metro, trans);
	 c->comms_time += update_halo(&comlist, EVEN, f.su2doublet, SU2DB);
	 checkerboard_sweep_su2doublet(f, p, c, ODD, metro, trans);
	 c->comms_time += update_halo(&comlist, ODD, f.su2doublet, SU2DB);
	}
	#endif
	#ifdef TRIPLET
	for (int k=0; k<p.update_su2triplet; k++) {
	 checkerboard_sweep_su2triplet(f, p, c, EVEN, metro, trans);
	 c->comms_time += update_halo(&comlist, EVEN, f.su2triplet, SU2TRIP);
	 checkerboard_sweep_su2triplet(f, p, c, ODD, metro, trans);
	 c->comms_time += update_halo(&comlist, ODD, f.su2triplet, SU2TRIP);
	}
	#endif

}

/** @file checkpoint.c
*
* Routines for storing lattice configuration in a file,
* and printing useful information about the simulation.
*
* TODO
*
*/

#include "su2.h"

/* Prints acceptance rates of relevant algorithms
* in a compact form. Called at each checkpoint.
*/
void print_acceptance(params p, counters c) {
	int print_gauge = 0;
	#ifdef TRIPLET
		print_gauge = 1;
	#endif
	printf("\n-- Acceptance rates --\n	");
	if (p.algorithm_su2link == METROPOLIS || print_gauge) {
		printf("SU(2) link %.2lf%%, ",
			100.0*c.accepted_su2link/c.total_su2link);
	}

	#ifdef HIGGS
		if (p.algorithm_su2doublet == METROPOLIS) {
			printf("Higgs %.2lf%%, ",
				100.0*c.accepted_doublet/c.total_doublet);

		} else if (p.algorithm_su2doublet == OVERRELAX) {
			printf("Higgs overrelax %.2lf%%, Higgs Metropolis %.2lf%%, ",
				100.0*c.acc_overrelax_doublet/c.total_overrelax_doublet,
				100.0*c.accepted_doublet/c.total_doublet);
		}
	#endif

	#ifdef TRIPLET
		if (p.algorithm_su2triplet == METROPOLIS) {
			printf("Triplet %.2lf%%, ",
				100.0*c.accepted_triplet/c.total_triplet);
		} else if (p.algorithm_su2triplet == OVERRELAX) {
			printf("Triplet overrelax %.2lf%%, Triplet Metropolis %.2lf%%, ",
				100.0*c.acc_overrelax_triplet/c.total_overrelax_triplet,
				100.0*c.accepted_triplet/c.total_triplet);
		}
	#endif
	if (p.multicanonical) {
		printf("multicanonical %.2lf%%, ",
				100.0*c.accepted_muca/c.total_muca);
	}
	printf("\b \b\n");

}


#include "su2.h"

int main(int argc, char *argv[]) {


	// structures for holding data
	params p;
	fields f;
	counters c;
	// initialize counters. move to a separate file?
	c.accepted_su2link = 0;
	c.accepted_doublet = 0;
	c.accepted_triplet = 0;
	c.acc_overrelax_doublet = 0;
	c.acc_overrelax_triplet = 0;
	c.total_su2link = 0;
	c.total_doublet = 0;
	c.total_triplet = 0;
	c.total_overrelax_doublet = 0;
	c.total_overrelax_triplet = 0;

	// read in the config file.
	// This needs to be done before allocating anything since we don't know the dimensions otherwise
	get_parameters(argv[1], &p);
	print_parameters(p);

	// initialize the site neighbor lookup pointers
	alloc_neighbors(&p);
	calculate_neighbors(&p);

	// initialize all fields
	alloc_fields(p, &f);
	setfields(f, p);

	// main iteration loop
	ulong iter = 1;
	ulong start = 1;
	// metro = 1 if we force metropolis sweep, 0 otherwise
	char metro = 0;
	// how often to force metropolis
	int metro_interval = 5;

	// labels for results file (TODO only if new simulation)
	print_labels();

	// possible optimization: keep track of changes to action and other observables
	// when doing updates on the lattice and use those in measurements, instead of
	// calculating everything from the scratch each time.
	// Drawback: the variables used to keep track accumulate errors from limited float precision.

	printf("Starting simulation!\n");
	while (iter <= p.iterations) {

		if ((iter % p.interval == 0) && (iter != start)) {
			measure(f, p);
		}
		if ((iter % p.checkpoint == 0) && (iter != start)) {
			printf("Checkpointing at iteration %lu. \n", iter);

			printf("---------------------- Metropolis acceptance rates ----------------------\n");
			if (p.algorithm_su2link == METROPOLIS) {
				printf("SU(2) link %lu/%lu (%.2lf%%), ", c.accepted_su2link, c.total_su2link,
					100.0*c.accepted_su2link/c.total_su2link);
			}
			#ifdef HIGGS
				if (p.algorithm_su2doublet == METROPOLIS) {
					printf("Higgs %lu/%lu (%.2lf%%), ", c.accepted_doublet, c.total_doublet,
						100.0*c.accepted_doublet/c.total_doublet);

				} else if (p.algorithm_su2doublet == OVERRELAX) {
					printf("Higgs overrelax %lu/%lu (%.2lf%%), Higgs Metropolis %lu/%lu (%.2lf%%), ", c.acc_overrelax_doublet, c.total_overrelax_doublet,
						100.0*c.acc_overrelax_doublet/c.total_overrelax_doublet, c.accepted_doublet, c.total_doublet,
						100.0*c.accepted_doublet/c.total_doublet);
				}
			#endif
			#ifdef TRIPLET
				if (p.algorithm_su2triplet == METROPOLIS) {
					printf("Triplet %lu/%lu (%.2lf%%), ", c.accepted_triplet, c.total_triplet,
						100.0*c.accepted_triplet/c.total_triplet);
				} else if (p.algorithm_su2doublet == OVERRELAX) {
					printf("Triplet overrelax %lu/%lu (%.2lf%%), Triplet Metropolis %lu/%lu (%.2lf%%), ", c.acc_overrelax_triplet, c.total_overrelax_triplet,
						100.0*c.acc_overrelax_triplet/c.total_overrelax_triplet, c.accepted_triplet, c.total_triplet,
						100.0*c.accepted_triplet/c.total_triplet);
				}
			#endif
			printf("\b\b \b.\n");
		}
		if (iter % metro_interval == 0) {
			metro = 1;
		} else {
			metro = 0;
		}

		update_lattice(f, p, &c, metro);
		iter++;
	}

	fclose(p.resultsfile);

	// free memory and finish
	free_fields(p, &f);
	free_neighbors(&p);

	free(p.L);
	free(p.parity);

  return 0;
}

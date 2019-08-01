
#include "su2.h"
#include "comms.h"
	

int main(int argc, char *argv[]) {


	// temp
	waittime = 0.0;

	// standard data structures
	params p;
	fields f;
	counters c;
	comlist_struct comlist;
	

	double time = 0;
	clock_t start_time, end_time;


	#ifdef MPI

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p.size); // how many MPI threads

  if(!p.rank)
    printf("MPI support enabled\n");

	printf0(p, "Starting %d MPI processes\n", p.size);
  MPI_Barrier(MPI_COMM_WORLD);

	#else // no MPI
	p.rank = 0;
	p.size = 1;
	#endif

	// read in the config file.
	// This needs to be done before allocating anything since we don't know the dimensions otherwise
	get_parameters(argv[1], &p);
	// print parameters from master node only (master = rank 0)
	if (!p.rank) {
		print_parameters(p);
	}

	// initialize parallel layout and lookup tables
	start_time = clock();

  layout(&p, &comlist);

	end_time = clock();
	time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	
	// initialize accept/reject/etc counters
	init_counters(&c);
	
	printf0(p, "Initialization done! Took %lf seconds.\n", time);


	// initialize all fields
	alloc_fields(p, &f);
	setfields(f, p);

	// main iteration loop
	long iter = 1;
	long start = 1;
	// metro = 1 if we force metropolis sweep, 0 otherwise
	char metro = 0;
	// how often to force metropolis
	int metro_interval = 5;

	// labels for results file (TODO only if new simulation)
	if (!p.rank) {
		print_labels();
	}
	// possible optimization: keep track of changes to action and other observables
	// when doing updates on the lattice and use those in measurements, instead of
	// calculating everything from the scratch each time.
	// Drawback: the variables used to keep track accumulate errors from limited float precision.
	if (!p.rank) {
		printf("Starting simulation!\n");
	}
	

	start_time = clock();
	
	while (iter <= p.iterations) {

		if (iter % p.interval == 0) {
			measure(f, p, &c);
		}
		if ((iter % p.checkpoint == 0) && (iter != start)) {
			// Checkpoint time
			end_time = clock();
			time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
			if (!p.rank) {
				// MPI barrier before checkpointing??
				printf0(p, "\nCheckpointing at iteration %lu. Total time so far: %.1lfs, of which %.2lf%% comms.\n", iter, time, 100.0*c.comms_time/time);

				printf("\n--Metropolis acceptance rates--\n	");
				if (p.algorithm_su2link == METROPOLIS) {
					printf0(p, "SU(2) link %lu/%lu (%.2lf%%), ", c.accepted_su2link, c.total_su2link,
						100.0*c.accepted_su2link/c.total_su2link);
				}
				#ifdef HIGGS
					if (p.algorithm_su2doublet == METROPOLIS) {
						printf0(p, "Higgs %lu/%lu (%.2lf%%), ", c.accepted_doublet, c.total_doublet,
							100.0*c.accepted_doublet/c.total_doublet);

					} else if (p.algorithm_su2doublet == OVERRELAX) {
						printf0(p, "Higgs overrelax %lu/%lu (%.2lf%%), Higgs Metropolis %lu/%lu (%.2lf%%), ", c.acc_overrelax_doublet, c.total_overrelax_doublet,
							100.0*c.acc_overrelax_doublet/c.total_overrelax_doublet, c.accepted_doublet, c.total_doublet,
							100.0*c.accepted_doublet/c.total_doublet);
					}
				#endif
				#ifdef TRIPLET
					if (p.algorithm_su2triplet == METROPOLIS) {
						printf0(p, "Triplet %lu/%lu (%.2lf%%), ", c.accepted_triplet, c.total_triplet,
							100.0*c.accepted_triplet/c.total_triplet);
					} else if (p.algorithm_su2doublet == OVERRELAX) {
						printf0(p, "Triplet overrelax %lu/%lu (%.2lf%%), Triplet Metropolis %lu/%lu (%.2lf%%), ", c.acc_overrelax_triplet, c.total_overrelax_triplet,
							100.0*c.acc_overrelax_triplet/c.total_overrelax_triplet, c.accepted_triplet, c.total_triplet,
							100.0*c.accepted_triplet/c.total_triplet);
					}
				#endif
				printf0(p, "\b\b \b\n");
			}	// end printing from master node
		}
		if (iter % metro_interval == 0) {
			metro = 1;
		} else {
			metro = 0;
		}

		update_lattice(f, p, comlist, &c, metro);
		iter++;
		
	} // end main loop

	if (p.resultsfile != NULL)	{
		fclose(p.resultsfile);
	}

	// free memory and finish
	free_fields(p, &f);
	free_comlist(&comlist);
	free_lattice_arrays(&p);

	end_time = clock();
	time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
	
  printf0(p, "Reached end! Total time taken: %.1lfs, of which %.2lf%% comms.\n", time, 100.0*c.comms_time/time);
	printf("Node %d: Time spent waiting: %.1lfs \n", p.rank, waittime);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

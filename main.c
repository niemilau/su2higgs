
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
	weight w;

	clock_t start_time, end_time;
	double time = 0.0;


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

	// read stuff for multicanonical. if non-multicanonical run, just sets dummy weight
	get_weight_parameters(argv[1], &p, &w);

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
	int metro_interval = 6;

	// labels for results file (TODO only if new simulation)
	if (!p.rank) {
		print_labels();
	}

	if (p.multicanonical) {
		// initialize multicanonical
		load_weight(p, &w);
		if (!w.readonly) {
			printf0(p, "readonly not 0, so multicanonical weight WILL be modified!\n");
		} else {
			printf0(p, "Read-only run, will not modify weight. \n");
		}
	}

	// possible optimization: keep track of changes to action and other observables
	// when doing updates on the lattice and use those in measurements, instead of
	// calculating everything from the scratch each time.
	// Drawback: the variables used to keep track accumulate errors from limited float precision.


	if (!p.rank) {
		printf("\nStarting simulation!\n");
	}

	#ifdef MPI
	// make sure weight is not written before all nodes get the initial weight.
	// should not happen, but just to be sure
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	start_time = clock();

	while (iter <= p.iterations) {

		if (iter % p.interval == 0) {
			measure(f, p, &c, &w);
		}
		if ((iter % p.checkpoint == 0) && (iter != start)) {
			// Checkpoint time
			end_time = clock();
			time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
			if (!p.rank) {
				printf("\nCheckpointing at iteration %lu. Total time: %.1lfs, %.2lf%% comms.\n", iter, time, 100.0*c.comms_time/time);
				print_acceptance(p, c);
			}

		// update max iterations if the config file has been changed by the user
		get_parameters(argv[1], &p);

		} // end checkpoint
		if (iter % metro_interval == 0) {
			metro = 1;
		} else {
			metro = 0;
		}

		update_lattice(f, p, &comlist, &w, &c, metro);
		iter++;

	} // end main loop

	if (p.resultsfile != NULL)	{
		fclose(p.resultsfile);
	}

	// free memory and finish
	free_fields(p, &f);
	free_comlist(&comlist);
	free_lattice_arrays(&p);
	if (p.multicanonical) {
		free_muca_arrays(&w);
	}

	end_time = clock();
	time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	printf("Node %d ready, time spent waiting: %.1lfs \n", p.rank, waittime);
	#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	#endif
	printf0(p, "Reached end! Total time taken: %.1lfs, of which %.2lf%% comms. Time per iteration: %.3lfs \n", time, 100.0*c.comms_time/time, time/p.iterations);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

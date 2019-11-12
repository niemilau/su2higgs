
#include "su2.h"
#include "comms.h"
#ifdef WALL
	#include "wallprofile.h"
#endif

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
	double timing = 0.0;


	#ifdef MPI

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p.size); // how many MPI threads

	printf0(p, "Starting %d MPI processes\n", p.size);
  MPI_Barrier(MPI_COMM_WORLD);

	#else // no MPI
	p.rank = 0;
	p.size = 1;
	#endif

	// set RNG seed: different for each node!
	long seed = (long) (time(NULL) * (p.rank + 1.0));
	srand48(seed);

	if (!p.rank) {
		printf("Seed in root node: %ld\n", seed);
	}

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
	timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	// initialize accept/reject/etc counters
	init_counters(&c);

	printf0(p, "Initialization done! Took %lf seconds.\n", timing);

	// initialize all fields
	alloc_fields(p, &f);

	// check if p.latticefile exists and load it; if not, call setfields()
	if (access(p.latticefile,R_OK) == 0) {
		// ok
		printf0(p, "\nLoading latticefile: %s\n", p.latticefile);
		load_lattice(p, f, &c);
	} else {
		printf0(p, "No latticefile found; starting with cold configuration.\n");
		setfields(f, p);
		p.reset = 1;
	}

	#ifdef WALL
		// setup wall. NB! this overrides any other field initializations
		prepare_wall(&f, p);
		measure_wall(&f, p);
	#endif

	// labels for results file
	if (!p.rank) {
		print_labels();
	}

	if (p.multicanonical) {
		// initialize multicanonical. Needs to come after field initializations
		load_weight(p, &w);
		alloc_backup_arrays(p, &f, w);
		calc_orderparam(p, f, &w, EVEN);
		calc_orderparam(p, f, &w, ODD);
		if (!w.readonly) {
			printf0(p, "readonly not 0, so multicanonical weight WILL be modified!\n");
		} else {
			printf0(p, "Read-only run, will not modify weight. \n");
		}
	}

	// main iteration loop
	long iter;
	if (p.reset) {
		iter = 1;
		printf0(p, "\nStarting new simulation!\n");
	} else {
		iter = c.iter + 1;
		printf0(p, "\nContinuing from iteration %ld!\n", iter-1);
	}
	// metro = 1 if we force metropolis sweep, 0 otherwise
	char metro = 0;
	// how often to force metropolis
	int metro_interval = 5;

	#ifdef MPI
	// make sure weight is not written before all nodes get the initial weight.
	// should not happen, but just to be sure
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	start_time = clock();

	while (iter <= p.iterations) {

		// measure & update fields first, then checkpoint if needed
		if (iter % p.interval == 0) {
			measure(f, p, &c, &w);
			#ifdef WALL
				measure_wall(&f, p);
			#endif
		}

		if (iter % metro_interval == 0) {
			metro = 1;
		} else {
			metro = 0;
		}

		// keep halos in sync...
		barrier();
		if (!p.multicanonical) {
			update_lattice(&f, p, &comlist, &c, metro);
		} else {
			update_lattice_muca(&f, p, &comlist, &w, &c, metro);
		}


		if ((iter % p.checkpoint == 0)) {
			// Checkpoint time; print acceptance and save fields to latticefile
			end_time = clock();
			timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

			start_time = clock(); // restart timer

			c.total_time += timing; c.iter = iter; // store for I/O
			if (!p.rank) {
				printf("\nCheckpointing at iteration %lu. Total time: %.1lfs, %.2lf%% comms.\n", iter, c.total_time, 100.0*c.comms_time/c.total_time);
				print_acceptance(p, c);
			}

			save_lattice(p, f, c);
			// update max iterations if the config file has been changed by the user
			read_updated_parameters(argv[1], &p);
		} // end checkpoint

		iter++;

	} // end main loop

	if (p.resultsfile != NULL)	{
		fclose(p.resultsfile);
	}

	end_time = clock();
	timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	// save final configuration
	c.total_time += timing; c.iter = iter;
	save_lattice(p, f, c);

	// free memory and finish
	free_fields(p, &f);
	free_comlist(&comlist);
	free_lattice_arrays(&p);
	if (p.multicanonical) {
		free_muca_arrays(&f, &w);
	}

	#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	#endif
	printf("Node %d ready, time spent waiting: %.1lfs \n", p.rank, waittime);
	printf0(p, "Reached end! Total time taken: %.1lfs, of which %.2lf%% comms. time per iteration: %.3lfs \n", c.total_time, 100.0*c.comms_time/c.total_time, c.total_time/p.iterations);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

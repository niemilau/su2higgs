
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

	// print usage if the arguments are invalid
	if (argc != 2) {
		printf0(p, "Usage: ./<program name> <config file>\n");
		die(0);
	}

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

	sync_halos(&f, &p, &comlist);

	#ifdef WALL
		if (p.reset) {
			// setup wall. This overrides any other field initializations!!
			prepare_wall(&f, &p, &comlist); // halos re-synced here
		}
		measure_wall(&f, &p);
	#endif

	// labels for results file
	if (!p.rank) {
		print_labels();
	}

	if (p.multicanonical) {
		// initialize multicanonical. Needs to come after field initializations
		load_weight(&p, &w);
		alloc_backup_arrays(&p, &f, &w);
		calc_orderparam(&p, &f, &w, EVEN);
		calc_orderparam(&p, &f, &w, ODD);
		if (!w.readonly) {
			printf0(p, "readonly not 0, so multicanonical weight WILL be modified!\n");
		} else {
			printf0(p, "Read-only run, will not modify weight. \n");
		}
	}

	// metro = 1 if we force metropolis sweep, 0 otherwise
	char metro = 0;
	// how often to force metropolis
	int metro_interval = 5;

	/* if no lattice file was given or if reset=1 in config,
	* start by thermalizing without multicanonical,
	* except if the WALL flag is set, in which case do
	* multicanonical here too to prevent initial wall from collapsing. */
	long iter = 1;
	int is_muca = 0;
	if (p.reset) {


		if (p.multicanonical) {
			is_muca = 1;
			#ifndef WALL
				p.multicanonical = 0;
			#endif
		}


		printf0(p, "Thermalizing %ld iterations\n", p.n_thermalize);
		start_time = clock();
		while (iter <= p.n_thermalize) {

			if (iter % metro_interval == 0) {
				metro = 1;
			} else {
				metro = 0;
			}

			barrier();
			update_lattice(&f, &p, &comlist, &c, &w, metro);
			iter++;
		}

		end_time = clock();
		timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
		printf0(p, "Thermalization done, took %lf seconds.\n", timing);

		// now reset iteration counter and turn muca back on, if necessary
		iter = 1;
		if (is_muca) {
			p.multicanonical = 1;
		}

		printf0(p, "\nStarting new simulation!\n");
	} else {
		iter = c.iter + 1;
		printf0(p, "\nContinuing from iteration %ld!\n", iter-1);
	}

	// make sure weight is not written before all nodes get the initial weight.
	// should not happen, but just to be sure
	barrier();

	start_time = clock();

	// main iteration loop
	while (iter <= p.iterations) {

		// measure & update fields first, then checkpoint if needed
		if (iter % p.interval == 0) {
			measure(&f, &p, &c, &w);
		}

		if (iter % metro_interval == 0) {
			metro = 1;
		} else {
			metro = 0;
		}

		// keep nodes in sync. without this things can go wrong if some node is much
		// faster than others, so that it sends new fields to others before they managed
		// to receive the earlier ones. Todo: optimize this
		barrier();
		// update all fields once. multicanonical checks are contained in sweep routines
		update_lattice(&f, &p, &comlist, &c, &w, metro);

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

			#ifdef WALL
				measure_wall(&f, &p);
			#endif

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

	c.total_time += timing; c.iter = iter;
	// save final configuration
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
	//printf("Node %d ready, time spent waiting: %.1lfs \n", p.rank, waittime);
	printf0(p, "\nReached end! Total time taken: %.1lfs, of which %.2lf%% comms. time per iteration: %.6lfs \n", c.total_time, 100.0*c.comms_time/c.total_time, c.total_time/p.iterations);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

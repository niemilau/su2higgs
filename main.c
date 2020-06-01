
#include "su2.h"
#include "comms.h"
#include "hb_trajectory.h"

int main(int argc, char *argv[]) {

	// temp
	waittime = 0.0;

	// initialize global time variables
	Global_comms_time = 0.0;
	Global_total_time = 0.0;

	// standard data structures
	params p;
	fields f;
	counters c;
	comlist_struct comlist;
	weight w;
	#ifdef HB_TRAJECTORY
		trajectory traj;
	#endif

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
	srand(seed+1); // seed also rand(), used when shuffling arrays


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
	alloc_fields(&p, &f);
	if (!p.rank)
		printf("Allocated memory for fields.\n");

	// check if p.latticefile exists and load it; if not, call setfields()
	if (access(p.latticefile,R_OK) == 0) {
		// ok
		printf0(p, "\nLoading latticefile: %s\n", p.latticefile);
		load_lattice(&p, &f, &c, &comlist); // also calls sync_halos()
		printf0(p, "Fields loaded succesfully.\n");
	} else {
		printf0(p, "No latticefile found; starting with cold configuration.\n");
		setfields(f, p);
		sync_halos(&f, &p, &comlist);
		p.reset = 1;
	}

	// by default, update ordering is not randomized
	p.random_sweeps = 0;

	#ifdef WALL
		if (p.reset) {
			// setup wall. This overrides any other field initializations!!
			prepare_wall(&f, &p, &comlist); // halos re-synced here
			measure_along_z(&f, &p, 0); // initial measurements, identifier=0
		}
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
			printf0(p, "Increment reduction factor: %lf\n", w.reduction_factor);
		} else {
			printf0(p, "Read-only run, will not modify weight. \n");
		}
	}

	#ifdef HB_TRAJECTORY
		if(!p.multicanonical) {
			 printf0(p, "\nTurn multicanonical on for realtime trajectories! Exiting...\n");
			 die(-44);
		} else {
			printf0(p, "\nReal time simulation: reading file \"realtime_config\"\n");
			read_realtime_config("realtime_config", &traj);
			printf0(p, "Heatbath trajectory mode every %ld iterations, %ld trajectories each\n", traj.mode_interval, traj.n_traj);
			printf0(p, "--- min %lf, %max %lf, measure interval %ld ---\n", traj.min, traj.max, traj.interval);

			if (p.algorithm_su2link != HEATBATH) {
				printf0(p, "\n--- WARNING: SU(2) update not using heatbath!\n");
			}
			#ifdef U1
			if (p.algorithm_u1link != HEATBATH) {
				printf0(p, "\n--- WARNING: U(1) update not using heatbath!\n");
			}
			#endif
		}
	#endif

	/* if no lattice file was given or if reset=1 in config,
	* start by thermalizing without updating multicanonical weight
	*/
	long iter = 1;
	int modify_weight = 0;
	if (p.reset) {

		if (p.multicanonical) {
			if (!w.readonly) {
				modify_weight = 1;
				w.readonly = 1;
			}
		}


		printf0(p, "\nThermalizing %ld iterations\n", p.n_thermalize);
		start_time = clock();
		while (iter <= p.n_thermalize) {

			barrier();
			update_lattice(&f, &p, &comlist, &c, &w);
			iter++;
		}

		end_time = clock();
		timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
		printf0(p, "Thermalization done, took %lf seconds.\n", timing);
		Global_total_time += timing;

		// now reset iteration and time counters and turn weight updating back on, if necessary
		iter = 1;
		init_counters(&c);
		if (modify_weight) {
			w.readonly = 0;
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

	int flow_id = 1;
	// main iteration loop
	while (iter <= p.iterations) {

		// measure & update fields first, then checkpoint if needed
		if (iter % p.interval == 0) {
			measure(p.resultsfile, &f, &p, &w);

			#ifdef GRADFLOW
			int do_flow = 1;
			if (do_flow) {
				double t_max = 1.0;
				double dt = 0.2;
				grad_flow(&p, &f, &comlist, &w, t_max, dt, flow_id);
				flow_id++;
			}
			#endif

		}
		#ifdef MEASURE_Z
			if (iter % p.meas_interval_z == 0) {
				measure_along_z(&f, &p, iter / p.meas_interval_z);
			}
		#endif

		// update all fields. multicanonical checks are contained in sweep routines
		update_lattice(&f, &p, &comlist, &c, &w);

		if (iter % p.checkpoint == 0) {
			// Checkpoint time; print acceptance and save fields to latticefile
			end_time = clock();
			timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

			start_time = clock(); // restart timer

			Global_total_time += timing;
			c.iter = iter; // store for I/O
			if (!p.rank) {
				printf("\nCheckpointing at iteration %lu. Total time: %.1lfs, %.2lf%% comms.\n",
							iter, Global_total_time, 100.0*Global_comms_time/Global_total_time);
				print_acceptance(p, c);
			}

			save_lattice(p, f, c);
			// update max iterations if the config file has been changed by the user
			read_updated_parameters(argv[1], &p);
		} // end checkpoint

		#ifdef HB_TRAJECTORY
			if (iter % traj.mode_interval == 0) {
				make_realtime_trajectories(&p, &f, &comlist, &c, &w, &traj);
			}
		#endif

		iter++;

	} // end main loop

	if (p.resultsfile != NULL)	{
		fclose(p.resultsfile);
	}

	end_time = clock();
	timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	Global_total_time += timing;
	c.iter = iter;
	// save final configuration
	save_lattice(p, f, c);

	// free memory and finish
	free_fields(&p, &f);
	if (!p.rank)
		printf("Freed memory allocated for fields.\n");

	free_comlist(&comlist);
	free_lattice_arrays(&p);
	if (p.multicanonical) {
		free_muca_arrays(&f, &w);
	}
	// miscellaneous frees for stuff not allocated in alloc.c
	#ifdef MEASURE_Z
		free_latticetable(p.site_at_z);
	#endif

	barrier();

	//printf("Node %d ready, time spent waiting: %.1lfs \n", p.rank, waittime);
	printf0(p, "\nReached end! Total time taken: %.1lfs, of which %.2lf%% comms. time per iteration: %.6lfs \n",
				Global_total_time, 100.0*Global_comms_time/Global_total_time, Global_total_time/p.iterations);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

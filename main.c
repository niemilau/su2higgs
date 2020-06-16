
#include "su2.h"
#include "comms.h"
#include "hb_trajectory.h"

int main(int argc, char *argv[]) {

	// temp
	waittime = 0.0;

	// initialize global variables
	Global_comms_time = 0.0;
	Global_total_time = 0.0;

	// standard data structures
	lattice l;
	params p;
	fields f;
	counters c;
	weight w;
	#ifdef HB_TRAJECTORY
		trajectory traj;
	#endif

	clock_t start_time, end_time;
	double timing = 0.0;


	#ifdef MPI

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &l.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &l.size); // how many MPI threads

	printf0(l, "Starting %d MPI processes\n", l.size);
  MPI_Barrier(MPI_COMM_WORLD);

	#else // no MPI
		l.rank = 0;
		l.size = 1;
	#endif

	// print usage if the arguments are invalid
	if (argc != 2) {
		printf0(l, "Usage: ./<program name> <config file>\n");
		die(0);
	}

	// set RNG seed: different for each node!
	long seed = (long) (time(NULL) * (l.rank + 1.0));
	srand48(seed);
	srand(seed+1); // seed also rand(), used when shuffling arrays


	if (!l.rank) {
		printf("Seed in root node: %ld\n", seed);
	}

	// read in the config file.
	// This needs to be done before allocating anything since we don't know the dimensions otherwise
	get_parameters(argv[1], &l, &p);
	// print parameters from master node only (master = rank 0)
	if (!l.rank) {
		print_parameters(l, p);
	}

	// read stuff for multicanonical. if non-multicanonical run, just sets dummy weight
	get_weight_parameters(argv[1], &l, &p, &w);

	// initialize parallel layout and lookup tables
	start_time = clock();

	int do_prints = 1;
  layout(&l, do_prints, p.run_checks);

	#ifdef MEASURE_Z
		print_z_labels(&l, &p);
		printf0(l, "Measuring profiles along direction %d every %d iterations\n", l.z_dir+1, p.meas_interval_z);
	#endif

	end_time = clock();
	timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	// initialize accept/reject/etc counters
	init_counters(&c);

	printf0(l, "Initialization done! Took %lf seconds.\n", timing);

	// initialize all fields
	alloc_fields(&l, &f);
	if (!l.rank)
		printf("Allocated memory for fields.\n");

	// check if p.latticefile exists and load it; if not, call setfields()
	if (access(p.latticefile,R_OK) == 0) {
		// ok
		printf0(l, "\nLoading latticefile: %s\n", p.latticefile);
		load_lattice(&l, &f, &c, p.latticefile); // also calls sync_halos()
		printf0(l, "Fields loaded succesfully.\n");
	} else {
		printf0(l, "No latticefile found; starting with cold configuration.\n");
		setfields(f, l, p);
		sync_halos(&l, &f);
		p.reset = 1;
	}

	// by default, update ordering is not randomized
	p.random_sweeps = 0;

	#ifdef MEASURE_Z
		if (p.reset) {
			// setup wall. This overrides any other field initializations!!
			prepare_wall(&l, &f, &p); // halos re-synced here
			measure_along_z(&l, &f, &p, 0); // initial measurements, identifier=0
		}
	#endif

	// labels for results file
	if (!l.rank) {
		print_labels();
		if (p.do_local_meas) {

			print_labels_local(&l, "labels_local");

		}
	}

	if (p.multicanonical) {
		// initialize multicanonical. Needs to come after field initializations
		load_weight(&l, &w);
		alloc_backup_arrays(&l, &f, &w);
		calc_orderparam(&l, &f, &p, &w, EVEN);
		calc_orderparam(&l, &f, &p, &w, ODD);
		if (!w.readonly) {
			printf0(l, "readonly not 0, so multicanonical weight WILL be modified!\n");
			printf0(l, "Increment reduction factor: %lf\n", w.reduction_factor);
		} else {
			printf0(l, "Read-only run, will not modify weight. \n");
		}
	}

	#ifdef HB_TRAJECTORY
		if(!p.multicanonical) {
			 printf0(l, "\nTurn multicanonical on for realtime trajectories! Exiting...\n");
			 die(-44);
		} else {
			printf0(l, "\nReal time simulation: reading file \"realtime_config\"\n");
			read_realtime_config("realtime_config", &traj);
			printf0(l, "Heatbath trajectory mode every %ld iterations, %ld trajectories each\n", traj.mode_interval, traj.n_traj);
			printf0(l, "--- min %lf, %max %lf, measure interval %ld ---\n", traj.min, traj.max, traj.interval);

			if (p.algorithm_su2link != HEATBATH) {
				printf0(l, "\n--- WARNING: SU(2) update not using heatbath!\n");
			}
			#ifdef U1
			if (p.algorithm_u1link != HEATBATH) {
				printf0(l, "\n--- WARNING: U(1) update not using heatbath!\n");
			}
			#endif
		}
	#endif

	#ifdef GRADFLOW
		if (p.do_flow) {
			printf0(l, "\n----- Gradient flow every %d iterations -----\n", p.flow_interval);
			printf0(l, "dt %lf	t_max %lf		meas_interval %d \n", p.flow_dt, p.flow_t_max, p.flow_meas_interval);
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


		printf0(l, "\nThermalizing %ld iterations\n", p.n_thermalize);
		start_time = clock();
		while (iter <= p.n_thermalize) {

			barrier();
			update_lattice(&l, &f, &p, &c, &w);
			iter++;
		}

		end_time = clock();
		timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
		printf0(l, "Thermalization done, took %lf seconds.\n", timing);
		Global_total_time += timing;

		// now reset iteration and time counters and turn weight updating back on, if necessary
		iter = 1;
		init_counters(&c);
		if (modify_weight) {
			w.readonly = 0;
		}

		printf0(l, "\nStarting new simulation!\n");
	} else {
		iter = c.iter + 1;
		printf0(l, "\nContinuing from iteration %ld!\n", iter-1);
	}

	// make sure weight is not written before all nodes get the initial weight.
	// should not happen, but just to be sure
	barrier();

	start_time = clock();

	int flow_id = 1; // only used for gradient flows
	// main iteration loop
	while (iter <= p.iterations) {

		// measure & update fields first, then checkpoint if needed
		if (iter % p.interval == 0) {
			measure(p.resultsfile, &l, &f, &p, &w);
		}
		#ifdef MEASURE_Z
			if (iter % p.meas_interval_z == 0) {
				measure_along_z(&l, &f, &p, iter / p.meas_interval_z);
			}
		#endif

		#ifdef GRADFLOW
			if (p.do_flow && iter % p.flow_interval == 0) {
				grad_flow(&l, &f, &p, &w, p.flow_t_max, p.flow_dt, flow_id);
				flow_id++;
			}
		#endif

		// update all fields. multicanonical checks are contained in sweep routines
		update_lattice(&l, &f, &p, &c, &w);

		if (iter % p.checkpoint == 0) {
			// Checkpoint time; print acceptance and save fields to latticefile
			end_time = clock();
			timing = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

			start_time = clock(); // restart timer

			Global_total_time += timing;
			c.iter = iter; // store for I/O
			if (!l.rank) {
				printf("\nCheckpointing at iteration %lu. Total time: %.1lfs, %.2lf%% comms.\n",
							iter, Global_total_time, 100.0*Global_comms_time/Global_total_time);
				print_acceptance(p, c);
			}

			save_lattice(&l, f, c, p.latticefile);
			// update max iterations if the config file has been changed by the user
			read_updated_parameters(argv[1], &l, &p);
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
	save_lattice(&l, f, c, p.latticefile);

	// free memory and finish
	free_fields(&l, &f);
	if (!l.rank)
		printf("Freed memory allocated for fields.\n");

	free_lattice(&l);
	if (p.multicanonical) {
		free_muca_arrays(&f, &w);
	}
	// miscellaneous frees for stuff not allocated in alloc.c
	#ifdef MEASURE_Z
		free_latticetable(l.site_at_z);
	#endif

	barrier();

	//printf("Node %d ready, time spent waiting: %.1lfs \n", p.rank, waittime);
	printf0(l, "\nReached end! Total time taken: %.1lfs, of which %.2lf%% comms. time per iteration: %.6lfs \n",
				Global_total_time, 100.0*Global_comms_time/Global_total_time, Global_total_time/p.iterations);
	#ifdef MPI
  MPI_Finalize();
	#endif

  return 0;
}

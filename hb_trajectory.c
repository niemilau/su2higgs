/** @file hb_trajectory

*/

#ifdef HB_TRAJECTORY

#include "su2.h"

/* Do a set of realtime trajectories. 'weight' struct is used only for calculating the order parameter;
* otherwise multicanonical is turned off. 'id' is the identifier for current set of trajectories */
void make_realtime_trajectories(lattice* l, fields const* f, params* p, counters* c,
      weight* w, trajectory* traj, int id) {

  /* enter "heatbath trajectory mode" by turning off multicanonical weighting
  * and randomizing order of gauge link updates */
  int do_acc = w->do_acceptance;
  int rand_sweeps = p->random_sweeps;
  w->do_acceptance = 0;
  p->random_sweeps = 1;

  int iter = 1;
  int current_traj = 1;
  double muca_param = w->param_value[EVEN] + w->param_value[ODD];

  // copy field configuration to a new struct
  fields f_traj;
  alloc_fields(l, &f_traj);
  copy_fields(l, f, &f_traj);
  sync_halos(l, &f_traj);

  // open trajectory file in root node and write header
  if (!l->rank) {
    traj->trajectoryfile = fopen("trajectory", "a");
    fprintf(traj->trajectoryfile, "\n----- Begin set %d, start value %lf ----- \n", id, muca_param);
    fprintf(traj->trajectoryfile, "\n--- Trajectory %d --- \n", current_traj);
  }

  while (current_traj <= traj->n_traj) {

    // perform measurements every traj.interval iterations.
    if (iter % traj->interval == 0 || iter == 1) {

      /* recalculate order parameter for measure() and to check if the trajectory completed
      * Need to do this here because all multicanonical calls are skipped in update_lattice() */
      calc_orderparam(l, &f_traj, p, w, EVEN); // updates EVEN contribution only
      muca_param = calc_orderparam(l, &f_traj, p, w, ODD); // updates ODD and returns the full value

      // measure(traj->trajectoryfile, l, &f_traj, p, w);

      /* 'iter' counts the time: one full update on the gauge fields per update_lattice() call */
      if (!l->rank) {
        fprintf(traj->trajectoryfile, "%d %.8lf\n", iter, muca_param);
      }

      if (muca_param > traj->max || muca_param < traj->min) {
        current_traj++;
        // trajectory done, revert back to the initial configuration and repeat

        if (current_traj <= traj->n_traj) {
          copy_fields(l, f, &f_traj);
          sync_halos(l, &f_traj);

          // write header for the next trajectory
          if (!l->rank) {
            fflush(traj->trajectoryfile);
            fprintf(traj->trajectoryfile, "\n--- Trajectory %d --- \n", current_traj);
          }
          iter = 1;
          continue;
        } // else: will automatically exit the while loop

      }
    } // end measure if

    /* update fields as usual, but now there is no multicanonical weighting */
    update_lattice(l, &f_traj, p, c, w);

    // additional fflush every so often, so we see if the trajectory gets stuck
    if (iter % 1000 == 0 && !l->rank) {
      fflush(traj->trajectoryfile);
    }
    iter++;
  }

  // trajectories done, exit "heatbath trajectory mode".
  // lattice loaded already when exiting trajectory loop
  w->do_acceptance = do_acc;
  p->random_sweeps = rand_sweeps;

  if (!l->rank)	{
    fclose(traj->trajectoryfile);
  }
  free_fields(l, &f_traj);

}

/* Read input for the realtime simulation.
* This is much like read_config() in parameters.c */
void read_realtime_config(char *filename, trajectory* traj) {

	FILE *config;
	char key[100];
  char value[100];
	char total[200];
  int ret;

  int set_n_traj = 0;
  int set_min = 0, set_max = 0;
  int set_interval = 0, set_mode_interval = 0;

	if(access(filename,R_OK) == 0) {
		config = fopen(filename, "r");

		while(!feof(config)) {

			if(fgets(total,198,config) == NULL) {
				// end of file?
				break;
			}
			ret = sscanf(total,"%99s%99s",key,value);
			if(ret == EOF) {
				continue;
			}
			if(key[0] == '#') {
				continue;
			}

			if(!strcasecmp(key,"n_traj")) {
				traj->n_traj = strtol(value,NULL,10);
				set_n_traj = 1;
			}
      else if(!strcasecmp(key,"interval")) {
				traj->interval = strtol(value,NULL,10);
				set_interval = 1;
			}
      else if(!strcasecmp(key,"mode_interval")) {
				traj->mode_interval = strtol(value,NULL,10);
				set_mode_interval = 1;
			}
			else if(!strcasecmp(key,"min")) {
				traj->min = strtod(value,NULL);
				set_min = 1;
			}
      else if(!strcasecmp(key,"max")) {
				traj->max = strtod(value,NULL);
				set_max = 1;
			}

		}
		fclose(config);
	} else {
    fprintf(stderr, "Could not open realtime file \"%s\"! Exiting...\n", filename);
    die(-2);
  }

  check_set(set_interval, "interval");
  check_set(set_mode_interval, "mode_interval");
  check_set(set_min, "min");
  check_set(set_max, "max");
  check_set(set_n_traj, "n_traj");
}


void print_realtime_params(params p, trajectory traj) {

  printf("Heatbath trajectory mode every %d iterations, %d trajectories each\n", traj.mode_interval, traj.n_traj);
  printf("--- min %lf, max %lf, measure interval %d ---\n", traj.min, traj.max, traj.interval);

  if (p.algorithm_su2link != HEATBATH) {
    printf("\n--- WARNING: SU(2) update not using heatbath!\n");
  }
  #ifdef U1
    if (p.algorithm_u1link != HEATBATH) {
      printf("\n--- WARNING: U(1) update not using heatbath!\n");
    }
  #endif
  fflush(stdout);
}

#endif

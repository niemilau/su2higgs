/** @file hb_trajectory

*/

#ifdef HB_TRAJECTORY

#include "su2.h"
#include "hb_trajectory.h"

void write_trajectory_header(params* p, trajectory* traj, int current_traj, double init_val) {
  // file should be open in root node (only in root!)
  if (p->rank == 0 && traj->trajectoryfile == NULL) {
    printf("Error in hb_trajectory.c; trajectory file not open!\n");
    die(-321);
  }

  if(!p->rank) {
    fprintf(traj->trajectoryfile, "\n--- Trajectory %d starting from %lf --- \n", current_traj, init_val);
  }
}


void make_realtime_trajectories(params* p, fields* f, comlist_struct* comlist, counters* c, weight* w, trajectory* traj) {

  // enter "heatbath trajectory mode" by turning off multicanonical weighting
  // and randomizing order of gauge link updates
  w->do_acceptance = 0;
  p->random_sweeps = 1;
  // backup initial configuration
  printf0(*p, "");
  save_lattice(*p, *f, *c);

  // open trajectory file in root node
  if (!p->rank) {
    traj->trajectoryfile = fopen("trajectory", "a");
  }

  int iter_traj = 1;
  int current_traj = 1;
  double muca_param;
  // recalculate muca parameter
  calc_orderparam(p, f, w, EVEN); // updates EVEN contribution only
  double init_value = calc_orderparam(p, f, w, ODD); // updates ODD and returns the full value

  printf0(*p, "Generating %d heatbath trajectories starting at order parameter value %lf\n",
      traj->n_traj, init_value);

  // header for the first trajectory
  write_trajectory_header(p, traj, current_traj, init_value);

  while (current_traj <= traj->n_traj) {

    // perform measurements every traj.interval iterations.
    if (iter_traj % traj->interval == 0 || iter_traj == 1) {
      measure(traj->trajectoryfile, f, p, c, w);
      // after measurement, check if the trajectory completed
      muca_param = w->param_value[EVEN] + w->param_value[ODD];
      if (muca_param > traj->max || muca_param < traj->min) {
        // trajectory done, revert back to the initial configuration and repeat
        load_lattice(p, f, c, comlist);
        current_traj++;
        if (current_traj <= traj->n_traj) {
          // write header for the next trajectory
          write_trajectory_header(p, traj, current_traj, init_value);
          iter_traj = 1;
          continue;
        } // else: will automatically exit the while loop

      }
    }

    // update fields as usual, but now there is no multicanonical weighting
    barrier();
    update_lattice(f, p, comlist, c, w);

    iter_traj++;
  }

  // trajectories done, exit "heatbath trajectory mode".
  // lattice loaded already when exiting trajectory loop
  w->do_acceptance = 1;
  p->random_sweeps = 0;
  // close trajectory file
  if (!p->rank)	{
    fclose(traj->trajectoryfile);
    printf("Trajectories ready\n");
  }

}

/* Read input for the realtime simulation.
* This is much like read_config() in parameters.c
*/
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

#endif

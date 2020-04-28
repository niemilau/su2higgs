/** @file hb_trajectory.h
*
* TODO
*
*/

#include "su2.h"

#ifndef HB_TRAJECTORY_H
#define HB_TRAJECTORY_H

#ifdef HB_TRAJECTORY

typedef struct {
  FILE *trajectoryfile;
  int n_traj; // how many trajectories to generate from each initial configuration
  int interval; // how often to measure in trajectory mode
  int mode_interval; // how often to switch to trajectory mode
  double min, max; // min and max values of the order parameter before trajectory is considered complete
} trajectory;

// hb_trajectory.c
void write_trajectory_header(params* p, trajectory* traj, int current_traj, double init_val);
void make_realtime_trajectories(params* p, fields* f, comlist_struct* comlist, counters* c, weight* w, trajectory* traj);
void read_realtime_config(char *filename, trajectory* traj);

#endif

#endif

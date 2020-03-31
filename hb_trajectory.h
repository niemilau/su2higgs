/** @file hb_trajectory.h
*
* TODO
*
*/

#include "su2.h"

#ifndef HB_TRAJECTORY_H
#define HB_TRAJECTORY_H

typedef struct {

  FILE *trajectoryfile;
  int n_traj; // how many trajectories to generate from each initial configuration

} trajectory;

void write_trajectory_header(trajectory* traj, double init_val);

void trajectory_mode_on(params* p) {

  // store lattice config and disable multicanonical.
  // add argument in measure() for output file; in heatbath mode set this to trajectoryfile instead of p->measurefile
  //do these in main instead??

  p->multicanonical = 0;

}

void trajectory_mode_off(params* p, int is_muca) {
  // turn muca back on if necessary
  if (is_muca) {
    p->multicanonical = 1;
  }
  // undo changes made by trajectory_mode_on() and reload lattice config

}

#endif

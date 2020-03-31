/** @file hb_trajectory

*/

#include "su2.h"
#include "hb_trajectory.h"

void write_trajectory_header(params* p, trajectory* traj, double init_val) {

  // file should be open in root node (only in root!)
  if (p->rank == 0 && traj->trajectoryfile == NULL) {
    printf("Error in hb_trajectory.c; trajectory file not open!\n");
    die(-321);
  }






}

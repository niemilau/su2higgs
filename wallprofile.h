/** @file wallprofile.h
*
* Includes stuff needed for wall profile routines.
* This file is intended to be left out from makefile if the routines are not needed.
*
* TODO
*
*/

#include "su2.h"
#include "comms.h"

#ifndef WALL_H
#define WALL_H

// define some nasty globals. these are initialized in layout()

/* wallcoord[z][x] gives the site index i of site with z-coord z,
* and x runs over all sites with that z coordinate. */
long** wallcoord;
long offset_z; // z coordinate of the MPI node times p.sliceL[p.dim-1]
long sites_per_z; // how many lattice sites for a given z coord
long wall_count; // how many times have we measured the wall


// wallprofile.c
void prepare_wall(fields* f, params p, comlist_struct* comlist);
void measure_wall(fields* f, params p);

#endif

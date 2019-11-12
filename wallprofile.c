/** @file wallprofile.c
*
* Routines for studying behavior of phase boundaries.
* Here we assume that for each z coordinate (last direction), there is always the
* same number of lattice sites.
*
* TODO polish this:
* 1) get rid of global variables
* 2) make it work for general lattice slicings
* 3) overall easier control over the z-coordinate
*
*/

#include "su2.h"
#include "wallprofile.h"


/* Create an initial configuration where the phase wall
*  is clearly visible. This needs to be adjusted depending on the theory.
* NB! need a long lattice in one direction, this routine assumes the longest
* direction is the "last" direction (z-coordinate).
*/
void prepare_wall(fields* f, params p) {

  wall_count = 0;
  long nz = p.sliceL[p.dim-1];
  for (long z=0; z<nz; z++) {

    for (long x=0; x<sites_per_z; x++) {

      long i = wallcoord[z][x];
      if (z + offset_z < 0.5 * p.L[p.dim-1]) {
        // for small z:
        f->su2doublet[i][0] = 0.2;
        f->su2doublet[i][1] = 0.0;
        f->su2doublet[i][2] = 0.0;
        f->su2doublet[i][3] = 0.0;
        #ifdef TRIPLET
        f->su2triplet[i][0] = 1.5;
        f->su2triplet[i][1] = 0.0;
        f->su2triplet[i][2] = 0.0;
        #endif
      } else {
        // for large z:
        f->su2doublet[i][0] = 1.5;
        f->su2doublet[i][1] = 0.0;
        f->su2doublet[i][2] = 0.0;
        f->su2doublet[i][3] = 0.0;
        #ifdef TRIPLET
        f->su2triplet[i][0] = 0.2;
        f->su2triplet[i][1] = 0.0;
        f->su2triplet[i][2] = 0.0;
        #endif
      }
    }
  }
  printf0(p, "Wall profile initialized.\n");

}

/* Measure some quantity along the z-axis and print to a file.
* This routine assumes that ONLY the z-direction is split into MPI nodes!!
*/
void measure_wall(fields* f, params p) {

  long nz = p.sliceL[p.dim-1];
  FILE* wallfile;

  barrier();

  char fname[200], temp[200];
  sprintf(temp, "%ld", wall_count);
  strcpy(fname, "wallprofile_");
  strcat(fname, temp);

  // temporary arrays for storing the observables as functions of z
  double* f1 = make_singletfield(nz);
  double* f2 = make_singletfield(nz);

  if (!p.rank)
    wallfile = fopen(fname, "a");

  for (long z=0; z<nz; z++) {
    double phi2 = 0.0;
    double Sigma2 = 0.0;

    for (long x=0; x<sites_per_z; x++) {
      long i = wallcoord[z][x];
      phi2 += doubletsq(f->su2doublet[i]);
      #ifdef TRIPLET
      Sigma2 += tripletsq(f->su2triplet[i]);
      #endif
    } // end x loop

    f1[z] = phi2 / ((double) sites_per_z);
    #ifdef TRIPLET
    f2[z] = Sigma2 / ((double) sites_per_z);
    #endif
  }

  #ifdef MPI
    // if not root, send
    if (p.rank != 0) {
      MPI_Send(f1, nz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(f2, nz, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
  #endif

  if (!p.rank) {
    // print to file from root node
    for (int rank=0; rank<p.size; rank++) {
      long offset;
      #ifdef MPI
      if (rank > 0) {
        MPI_Recv(f1, nz, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(f2, nz, MPI_DOUBLE, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        offset = rank * nz;
      } else {
        offset = 0;
      }
      #endif

      for (long z=0; z<nz; z++) {
        fprintf(wallfile, "%ld %g", z+offset, f1[z]);
        #ifdef TRIPLET
        fprintf(wallfile, " %g", f2[z]);
        #endif
        fprintf(wallfile, "\n");
      }
    }
  }


  free_singletfield(f1);
  free_singletfield(f2);

  if (!p.rank) {
    fclose(wallfile);
    wall_count++;
  }
  barrier();
}
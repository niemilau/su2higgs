/** @file z_coord.c
*
* Routines for measuring stuff along a given direction, denoted "z direction".
* Here I take z = longest direction on the lattice, but can be changed in init_z_coord()
*/

#ifdef MEASURE_Z // makefile flag, do nothing if not defined

#include "su2.h"

/* Initialize z-coord variables in params struct.
* Call in layout() after sitemap() and all other site remappings */
void init_z_coord(params* p) {
  // which dir?
  int longest=p->dim-1;
  for (int dir=p->dim-1; dir>=0; dir--) {
    if (p->L[dir] > p->L[longest]) longest = dir;
  }
  p->z_dir = longest;

  /* offset: each node normally loops over 0<= z < p.sliceL[p.z_dir].
  * Then p.offset_z + z should give the physical z coordinate */
  long xnode[p->dim]; // coordinates of the MPI node
  indexToCoords(p->dim, p->nslices, p->rank, xnode);

  p->offset_z = xnode[p->z_dir] * p->sliceL[p->z_dir];

  // site_at_z
  for (long nz=0; nz<p->sliceL[p->dim-1]; nz++) {
    long tot = 0;
    for (long i=0; i<p->sites; i++) {
      if (nz + offset_z == p->coords[i][p->dim-1]) {
        wallcoord[nz][tot] = i;
        tot++;
      }
    }
    if (tot != sites_per_z) {
      printf("Error in wall profile routines in layout.c!\n");
      die(420);
    }
  }


}

/* Create an initial configuration where the phase wall
*  is clearly visible. This needs to be adjusted depending on the theory.
* NB! need a long lattice in one direction, this routine assumes the longest
* direction is the "last" direction (z-coordinate).
*/
void prepare_wall(fields* f, params const* p, comlist_struct* comlist) {

  wall_count = 0;
  long nz = p->sliceL[p->dim-1];
  for (long z=0; z<nz; z++) {

    for (long x=0; x<sites_per_z; x++) {

      long i = wallcoord[z][x];
      if (z + offset_z < 0.5 * p->L[p->dim-1]) {
        // for small z:
        #ifdef HIGGS
        f->su2doublet[i][0] = 0.2 + 0.01*drand48();
        f->su2doublet[i][1] = 0.01*drand48();
        f->su2doublet[i][2] = 0.01*drand48();
        f->su2doublet[i][3] = 0.01*drand48();
        #endif
        #ifdef TRIPLET
        f->su2triplet[i][0] = 1.5 + 0.05*drand48();
        f->su2triplet[i][1] = 0.05*drand48();
        f->su2triplet[i][2] = 0.05*drand48();
        #endif
        // also set gauge links to (hopefully) help with thermalization
        for (int dir=0; dir<p->dim; dir++) {
          double u = 1.0 - 0.03*drand48();
    			f->su2link[i][dir][0] = u;
    			f->su2link[i][dir][1] = sqrt((double)(1.0 - u*u));
    			f->su2link[i][dir][2] = 0.0;
    			f->su2link[i][dir][3] = 0.0;
    		}
      } else {
        // for large z:
        #ifdef HIGGS
        f->su2doublet[i][0] = 1.5 + 0.05*drand48();
        f->su2doublet[i][1] = 0.05*drand48();
        f->su2doublet[i][2] = 0.05*drand48();
        f->su2doublet[i][3] = 0.05*drand48();
        #endif
        #ifdef TRIPLET
        f->su2triplet[i][0] = 0.2 + 0.01*drand48();
        f->su2triplet[i][1] = 0.01*drand48();
        f->su2triplet[i][2] = 0.01*drand48();
        #endif
        // gauge links:
        for (int dir=0; dir<p->dim; dir++) {
          double u = 1.0 - 0.4*drand48();
    			f->su2link[i][dir][0] = u;
    			f->su2link[i][dir][1] = sqrt((double)(1.0 - u*u));
    			f->su2link[i][dir][2] = 0.0;
    			f->su2link[i][dir][3] = 0.0;
    		}
      }
    }
  }
  // wall initialized, now just need to sync halo fields
  sync_halos(f, p, comlist);

  printf0(*p, "Wall profile initialized.\n");

}

/* Measure some quantity along the z-axis and print to a file.
* This routine assumes that ONLY the z-direction is split into MPI nodes!!
*/
void measure_wall(fields const* f, params const* p) {

  long nz = p->sliceL[p->dim-1];
  FILE* wallfile;

  barrier();

  // temporary arrays for storing the observables as functions of z
  double* f1 = make_singletfield(nz);
  double* f2 = make_singletfield(nz);

  if (!p->rank)
    wallfile = fopen("wallprofile", "w");

  for (long z=0; z<nz; z++) {
    double phi2 = 0.0;
    double Sigma2 = 0.0;

    for (long x=0; x<sites_per_z; x++) {
      long i = wallcoord[z][x];
      #ifdef HIGGS
      phi2 += doubletsq(f->su2doublet[i]);
      #endif
      #ifdef TRIPLET
      Sigma2 += tripletsq(f->su2triplet[i]);
      #endif
    } // end x loop
    #ifdef HIGGS
    f1[z] = phi2 / ((double) sites_per_z);
    #endif
    #ifdef TRIPLET
    f2[z] = Sigma2 / ((double) sites_per_z);
    #endif
  }

  #ifdef MPI
    // if not root, send
    if (p->rank != 0) {
      MPI_Send(f1, nz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(f2, nz, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
  #endif

  if (!p->rank) {
    // print to file from root node
    for (int rank=0; rank<p->size; rank++) {
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

  if (!p->rank) {
    fclose(wallfile);
    wall_count++;
  }
  barrier();
}

#endif

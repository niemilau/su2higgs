/** @file z_coord.c
*
* Routines for measuring stuff along a given direction, denoted "z direction".
* Here I take z = longest direction on the lattice, but can be changed in init_z_coord().
*
* It is assumed that all MPI nodes have the same size, specifically,
* p.sliceL[p.z_dir] and p.sites_per_z should be the same in all nodes.
* This is automatic in the current implementation of make_slices().
* Assumption could be lifted by probing the MPI message sizes before receiving.
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

  // how many sites for each z coordinate (on node and the total area of the "x,y plane")
  p->sites_per_z = 1;
  p->area = 1;
  if (p->dim == 1) {
    p->sites_per_z = p->sliceL[p->z_dir];
    p->area = p->L[p->z_dir];
  }
  else {

    for (int dir=0; dir<p->dim; dir++) {
      if (dir != p->z_dir) {
        p->sites_per_z *= p->sliceL[dir];
        p->area *= p->L[dir];
      }
    }

  }

  /* offset: each node normally loops over 0<= z < p.sliceL[p.z_dir].
  * Then p.offset_z + z should give the physical z coordinate */
  long xnode[p->dim]; // coordinates of the MPI node
  indexToCoords(p->dim, p->nslices, p->rank, xnode);

  p->offset_z = xnode[p->z_dir] * p->sliceL[p->z_dir];

  // site_at_z[z] is a list of all sites with z coordinate offset_z + z
  p->site_at_z = alloc_latticetable(p->sliceL[p->z_dir], p->sites_per_z);

  for (long nz=0; nz<p->sliceL[p->z_dir]; nz++) {
    long tot = 0;
    for (long i=0; i<p->sites; i++) {
      if (nz + p->offset_z == p->coords[i][p->z_dir]) {
        p->site_at_z[nz][tot] = i;
        tot++;
      }
    }
    if (tot != p->sites_per_z) {
      printf("Error counting sites in z_coord.c!\n");
      die(420);
    }
  }

  // test that all nodes have the same sites_per_z
  long test_z = p->sites_per_z;
  bcast_long(&test_z);
  if (test_z != p->sites_per_z) {
    printf("Error in z_coord.c! Node %d has sites_per_z = %ld, while root node has %ld\n", p->rank, p->sites_per_z, test_z);
    die(421);
  }

  // measurement labels and p.n_meas_z
  print_z_labels(p);

  printf0(*p, "Measuring profiles along direction %d every %d iterations\n", p->z_dir+1, p->meas_interval_z);

}

/* Labels for the measure_z file. Also sets p.n_meas_z.
*/
void print_z_labels(params* p) {
  int k = 1;

  if (!p->rank) {
  	FILE* f = fopen("labels_z", "w+");

    fprintf(f, "%d z\n", k); k++;
  	#ifdef HIGGS
  		fprintf(f, "%d phi^2\n", k); k++;
  	#endif
  	#ifdef TRIPLET
  		fprintf(f, "%d Sigma^2\n", k); k++;
  	#endif

  	fclose(f);
  }
  bcast_int(&k);
  p->n_meas_z = k-2; // k started from 1, not 0, and don't count the first column "z"
}

/* Measure stuff along the z axis.
* Each node works on arrays of size p.L[z] (full lattice, not slice!),
* and fills in only the elements in their own z range. These are then
* combined using reduce_sum(). We do p.n_meas_z different measurements
* so there are p.n_meas_z such arrays on each node.
* id is the identifier for the current measurement.
* TODO optimize communications, should be faster to send the whole arrays
* instead of reducing sums in each element separately.
*/
void measure_along_z(fields const* f, params const* p, long id) {
  long z_max = p->L[p->z_dir];
  long z_slice = p->sliceL[p->z_dir];
  FILE* file;
  int k; long z;

  if (p->n_meas_z <= 0) {
    return; // nothing to measure
  }
  /* arrays for storing the observables as functions of z
  * Can be understood as a field with z_max sites and n_meas_z DOFs */
  double** meas = make_field(z_max, p->n_meas_z);
  // initialize to 0. will remain 0 outside the range of z in my node
  for (z=0; z<z_max; z++) {
    for (k=0; k<p->n_meas_z; k++) {
      meas[z][k] = 0.0;
    }
  }

  // measure something at each z
  for (long z=0; z<z_slice; z++) {
    double phisq = 0.0;
    double Sigmasq = 0.0;

    for (long x=0; x<p->sites_per_z; x++) {
      long i = p->site_at_z[z][x];

      #ifdef HIGGS
        phisq += doubletsq(f->su2doublet[i]);
      #endif
      #ifdef TRIPLET
        Sigmasq += tripletsq(f->su2triplet[i]);
      #endif
    } // end area loop

    /* store using the physical, full z coordinate.
    * Use same ordering here as in labels_z() */
    k = 0;
    #ifdef HIGGS
      meas[z + p->offset_z][k] = phisq; k++;
    #endif
    #ifdef TRIPLET
      meas[z + p->offset_z][k] = Sigmasq; k++;
    #endif

  } // end z loop

  if (!p->rank) {
    file = fopen("measure_z", "a");
    // write header for the current set of measurements
    fprintf(file, "\nMeasurement %ld:\n", id);
  }

  // now reduce and write to file
  for (z=0; z<z_max; z++) {
    for (k=0; k<p->n_meas_z; k++) {
      meas[z][k] = reduce_sum(meas[z][k]);
    }
    k = 0;
    if (!p->rank) {
      fprintf(file, "%ld ", z);
      #ifdef HIGGS
        fprintf(file, "%g ", meas[z][k] / ((double)p->area) ); k++;
      #endif
      #ifdef TRIPLET
        fprintf(file, "%g ", meas[z][k] / ((double)p->area) ); k++;
      #endif
      fprintf(file, "\n");
  		fflush(file);
    }
  }

  // done
  if (!p->rank) {
    fclose(file);
  }

  free_field(meas);
}


/* Create an initial configuration where two phases coexist, separated by
* an interface at z = L[z_dir] / 2. The field values here need to be
* adjusted depending on the field content.
*/
void prepare_wall(fields* f, params const* p, comlist_struct* comlist) {

  long z_max = p->sliceL[p->z_dir];
  for (long z=0; z<z_max; z++) {

    for (long x=0; x<p->sites_per_z; x++) {

      long i = p->site_at_z[z][x];
      if (z + p->offset_z < 0.5 * p->L[p->z_dir]) {
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


#endif
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
void init_z_coord(lattice* l) {
  // which dir?
  int longest=l->dim-1;
  for (int dir=l->dim-1; dir>=0; dir--) {
    if (l->L[dir] > l->L[longest]) longest = dir;
  }
  l->z_dir = longest;

  // how many sites for each z coordinate (on node and the total area of the "x,y plane")
  l->sites_per_z = 1;
  l->area = 1;
  if (l->dim == 1) {
    l->sites_per_z = l->sliceL[l->z_dir];
    l->area = l->L[l->z_dir];
  }
  else {

    for (int dir=0; dir<l->dim; dir++) {
      if (dir != l->z_dir) {
        l->sites_per_z *= l->sliceL[dir];
        l->area *= l->L[dir];
      }
    }

  }

  /* offset: each node normally loops over 0<= z < p.sliceL[p.z_dir].
  * Then p.offset_z + z should give the physical z coordinate */

  l->offset_z = l->offset[l->z_dir];

  // site_at_z[z] is a list of all sites with z coordinate offset_z + z
  l->site_at_z = alloc_latticetable(l->sites_per_z, l->sliceL[l->z_dir]);

  for (long nz=0; nz<l->sliceL[l->z_dir]; nz++) {
    long tot = 0;
    for (long i=0; i<l->sites; i++) {
      if (nz + l->offset_z == l->coords[i][l->z_dir]) {
        l->site_at_z[nz][tot] = i;
        tot++;
      }
    }
    if (tot != l->sites_per_z) {
      printf("Error counting sites in z_coord.c!\n");
      die(420);
    }
  }

  // test that all nodes have the same sites_per_z
  long test_z = l->sites_per_z;
  bcast_long(&test_z, l->comm);
  if (test_z != l->sites_per_z) {
    printf("Error in z_coord.c! Node %d has sites_per_z = %ld, while root node has %ld\n", l->rank, l->sites_per_z, test_z);
    die(421);
  }

}

/* Labels for the measure_z file. Also sets l.n_meas_z.
*/
void print_z_labels(lattice const* l, params* p) {
  int k = 1;

  if (!l->rank) {
  	FILE* f = fopen("labels_z", "w+");

    fprintf(f, "%d z\n", k); k++;
  	#if (NHIGGS > 0)
  		fprintf(f, "%d phi^2\n", k); k++;
  	#endif
  	#ifdef TRIPLET
  		fprintf(f, "%d Sigma^2\n", k); k++;
  	#endif
    #ifdef SINGLET
  		fprintf(f, "%d S\n", k); k++;
      fprintf(f, "%d S^2\n", k); k++;
  	#endif


  	fclose(f);
  }
  bcast_int(&k, l->comm);
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
/*
void measure_along_z_old(lattice const* l, fields const* f, params const* p, long id) {
  long z_max = l->L[l->z_dir];
  long z_slice = l->sliceL[l->z_dir];
  FILE* file;
  int k; long z;

  if (p->n_meas_z <= 0) {
    return; // nothing to measure
  }

  // Arrays for storing the observables as functions of z
  // Can be understood as a field with z_max sites and n_meas_z DOFs
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

    for (long x=0; x<l->sites_per_z; x++) {
      long i = l->site_at_z[z][x];

      #ifdef HIGGS
        phisq += doubletsq(f->su2doublet[i]);
      #endif
      #ifdef TRIPLET
        Sigmasq += tripletsq(f->su2triplet[i]);
      #endif
    } // end area loop

    // store using the physical, full z coordinate.
    // Use same ordering here as in labels_z()
    k = 0;
    #ifdef HIGGS
      meas[z + l->offset_z][k] = phisq; k++;
    #endif
    #ifdef TRIPLET
      meas[z + l->offset_z][k] = Sigmasq; k++;
    #endif

  } // end z loop

  if (!l->rank) {
    file = fopen("measure_z", "a");
    // write header for the current set of measurements
    fprintf(file, "\nMeasurement %ld:\n", id);
  }

  // now combine results from all nodes
  for (z=0; z<z_max; z++) {
    for (k=0; k<p->n_meas_z; k++) {
      meas[z][k] = reduce_sum(meas[z][k], l->comm);
    }
  }

  // write to file
  for (z=0; z<z_max; z++) {
    k = 0;
    if (!l->rank) {
      fprintf(file, "%ld ", z);
      #ifdef HIGGS
        meas[z][k] = meas[z][k] / ((double)l->area);
        fprintf(file, "%g ", meas[z][k]);
        k++;
      #endif
      #ifdef TRIPLET
        fprintf(file, "%g ", meas[z][k] / ((double)l->area) ); k++;
      #endif
      fprintf(file, "\n");
  		fflush(file);
    }
  }

  // done
  if (!l->rank) {
    fclose(file);
  }

  free_field(meas);
}
*/

/* Measure stuff along the z axis.
* Each node works on arrays of size p.L[z] (full lattice, not slice!),
* and fills in only the elements in their own z range. These are then
* combined using reduce_sum(). We do p.n_meas_z different measurements
* so there are p.n_meas_z such arrays on each node. */
void measure_along_z(lattice const* l, fields const* f, params const* p, long id) {
  long z_max = l->L[l->z_dir];
  long z_slice = l->sliceL[l->z_dir];
  FILE* file;
  long z;

  // how many things to measure. This is set by print_z_labels()
  int n_meas_z = p->n_meas_z;


  /* arrays for storing the observables as functions of z
  * Can be understood as a field with z_max sites and n_meas_z DOFs */
  double** meas = make_field(z_max, n_meas_z);
  // initialize to 0. will remain 0 outside the range of z in my node
  for (z=0; z<z_max; z++) {
    for (int k=0; k<n_meas_z; k++) {
      meas[z][k] = 0.0;
    }
  }

  // measure stuff at each z
  for (long z=0; z<z_slice; z++) {

    for (long x=0; x<l->sites_per_z; x++) {

      // NB: in print_z_labels(), k=0 is the z coord, but that was not counted in n_meas_z so I start the measurements from k=0 here
      int k = 0; 

      long i = l->site_at_z[z][x];

      // Store measurements using the physical, full z coordinate. 
      // Also, use same ordering here as in print_z_labels().
      long z_phys = z + l->offset_z;

      #if (NHIGGS > 0)
        // phi^2
        meas[z_phys][k] += doubletsq(f->su2doublet[0][i]);
        k++;
      #endif

      #ifdef TRIPLET
        meas[z_phys][k] += tripletsq(f->su2triplet[i]);
        k++;
      #endif

      #ifdef SINGLET
        // S
        double s = f->singlet[i][0];
        meas[z_phys][k] += s;
        k++;

        // S^2
        meas[z_phys][k] += s*s;
        k++;
        
      #endif
    } // end area loop

  } // end z loop

  if (!l->rank) {
    file = fopen("measure_z", "a");
    // write header for the current set of measurements
    //fprintf(file, "\nMeasurement %ld:\n", id);
  }

  // now combine results from all nodes
  for (z=0; z<z_max; z++) {
    for (int k=0; k<n_meas_z; k++) {
      meas[z][k] = reduce_sum(meas[z][k], l->comm);
    }
  }

  // write plane averaged measurements to file
  for (z=0; z<z_max; z++) {
    if (!l->rank) {
      fprintf(file, "%ld ", z);

      for (int k=0; k<n_meas_z; k++) {
        fprintf(file, "%g ", meas[z][k] / ((double)l->area) );
      }

      fprintf(file, "\n");
    }
  }

  // done
  if (!l->rank) {
    fclose(file);
  }

  free_field(meas);
}


/* Create an initial configuration where two phases coexist, separated by
* an interface at z = L[z_dir] / 2. The field values here need to be
* adjusted depending on the field content. */
void prepare_wall(lattice* l, fields* f, params const* p) {

  long z_max = l->sliceL[l->z_dir];
  for (long z=0; z<z_max; z++) {

    for (long x=0; x<l->sites_per_z; x++) {

      long i = l->site_at_z[z][x];
      if (z + l->offset_z < 0.5 * l->L[l->z_dir]) {
        // for small z:
        #if (NHIGGS > 0)
          for (int db=0; db<NHIGGS; db++) {
            f->su2doublet[db][i][0] = 0.2 + 0.01*dran();
            f->su2doublet[db][i][1] = 0.01*dran();
            f->su2doublet[db][i][2] = 0.01*dran();
            f->su2doublet[db][i][3] = 0.01*dran();
          }
          
        #endif
        #ifdef TRIPLET
          f->su2triplet[i][0] = 1.5 + 0.05*dran();
          f->su2triplet[i][1] = 0.05*dran();
          f->su2triplet[i][2] = 0.05*dran();
        #endif
        #ifdef SINGLET
          f->singlet[i][0] = 0.8 + 0.01*dran();
        #endif
        // also set gauge links to (hopefully) help with thermalization
        for (int dir=0; dir<l->dim; dir++) {
          double u = 1.0 - 0.03*dran();
    			f->su2link[i][dir][0] = u;
    			f->su2link[i][dir][1] = sqrt((double)(1.0 - u*u));
    			f->su2link[i][dir][2] = 0.0;
    			f->su2link[i][dir][3] = 0.0;
    		}

      } else {
        // for large z:

        #if (NHIGGS > 0)
          for (int db=0; db<NHIGGS; db++) {
            f->su2doublet[db][i][0] = 1.2 + 0.05*dran();
            f->su2doublet[db][i][1] = 0.05*dran();
            f->su2doublet[db][i][2] = 0.05*dran();
            f->su2doublet[db][i][3] = 0.05*dran();
          }
          
        #endif
        #ifdef TRIPLET
          f->su2triplet[i][0] = 0.2 + 0.01*dran();
          f->su2triplet[i][1] = 0.01*dran();
          f->su2triplet[i][2] = 0.01*dran();
        #endif
        #ifdef SINGLET
          f->singlet[i][0] = 0.1 + 0.01*dran();
        #endif
        // gauge links:
        for (int dir=0; dir<l->dim; dir++) {
          double u = 1.0 - 0.4*dran();
    			f->su2link[i][dir][0] = u;
    			f->su2link[i][dir][1] = sqrt((double)(1.0 - u*u));
    			f->su2link[i][dir][2] = 0.0;
    			f->su2link[i][dir][3] = 0.0;
    		}
      }
    }
  }
  // wall initialized, now just need to sync halo fields
  sync_halos(l, f);

  printf0("Wall profile initialized.\n");

}

void prepare_wall_2hdm(lattice* l, fields* f, params const* p) {

  #if (NHIGGS>=2)

  long z_max = l->sliceL[l->z_dir];
  for (long z=0; z<z_max; z++) {

    for (long x=0; x<l->sites_per_z; x++) {

      // only sets wall for phi_2

      long i = l->site_at_z[z][x];
      if (z + l->offset_z < 0.5 * l->L[l->z_dir]) {
        // for small z:

        f->su2doublet[1][i][0] = 0.2 + 0.01*drand48();
        for (int a=1; a<SU2DB; a++) {
          f->su2doublet[1][i][a] = 0.01*drand48();
        }

      } else {
        // for large z:
        f->su2doublet[1][i][0] = 0.7 + 0.01*drand48();
        for (int a=1; a<SU2DB; a++) {
          f->su2doublet[1][i][a] = 0.05*drand48();
        }

      }
    }
  }
  // wall initialized, now just need to sync halo fields
  sync_halos(l, f);

  printf0(*l, "2HDM wall profile initialized.\n");

  #endif

}



#endif

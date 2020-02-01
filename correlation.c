/** @file correlation.c
*
* Routines for calculating correlation lengths.
* We use the observables described in hep-lat/9510020.
* In hep-lat/9612006, an optimized "blocking" method is described;
* this is not implemented here.
*
* Uses the same z-coordinate methods as wallprofile.c,
* so need to include wallprofile.h. Todo: better implementation
*
*/

#include "su2.h"
#include "wallprofile.h"

/*
* For Higgs, measure h(l) = 1/V * \sum_z R^2(z) R^2(z+l)
* where R^2 = 0.5 \sum_{x,y} Tr\he\Phi\Phi at fixed z.
* Higgs mass is found from the exponential fall-off off h(l), i.e.
* plot h(l) and fit h ~ e^{-m l}
*/

void measure_correlation(fields* f, params p) {

  long nz = p.sliceL[p.dim-1];
  FILE* corrfile;

  barrier();

  // temporary arrays for storing the observables as functions of z
  double* f1 = make_singletfield(nz);
  double* f2 = make_singletfield(nz);

  if (!p.rank)
    corrfile = fopen("correlation", "a");

  for (long z=0; z<nz; z++) {
    double hl = 0.0;
    /* for triplet, measure
    * s(l) = 1/V * \sum_z R^2(z) R^2(z+l)
    * with R^2(z) = \sum_{x,y} Tr \Sigma^2  */
    double sl = 0.0;

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
    fclose(corrfile);
  }
  barrier();
}

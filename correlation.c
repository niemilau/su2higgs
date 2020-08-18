/** @file correlation.c
*
* Routines for calculating correlation lengths.
*
* //TODO
*/

#ifdef CORRELATORS

#include "su2.h"

/* Calculate the total of some quantity over the all sites (in my node!) at location x_dir.
* The first argument is a pointer to the function that is called at each site whose physical coordinate
* equals x in the direction dir. The function operates on field[site].
* Note: does not divide by the number of sites */
double plane_sum(double (*funct)(double*), double** field, lattice* l, long x, int dir) {

  double res = 0.0;
  long x_node = x - l->offset[dir]; // coordinate on my node
  int skip = 0;

  if (x_node < 0 || x_node >= l->sliceL[dir]) {
    skip = 1; // out of bounds, so my node does not contribute
  }
  if (!skip) {
    for (long i=0; i<l->sites_per_coord[dir]; i++) {
      long site = l->sites_at_coord[dir][x_node][i];
      res += (*funct)(field[site]);
    }
  }

  res = allreduce(res, l->comm);
  return res;
}

#if (NHIGGS > 0)
/* Calculate H(l) = (1/V) \sum_z h(z) h(z+l),
* where h(z) = \sum_{x,y} Tr\Phi(z)\Phi(z).
* Below the argument d = l and dir specifies which direction z lives in. */
double higgs_correlator(lattice* l, fields const* f, int d, int dir, int higgs_id) {

  double res = 0.0;
  for (int z=0; z<l->L[dir]; z++) {

    /* periodicity: need to find site at coordinate z + d.
    * if x + d >= lattice length, we measure correlation in the "negative"
    * direction instead. */
    int zd = (z + d) % l->L[dir];

    double h1 = plane_sum(doubletsq, f->su2doublet[higgs_id], l, z, dir); // h1 <- h(z)
    double h2 = plane_sum(doubletsq, f->su2doublet[higgs_id], l, zd, dir); // h2 <- h(z+l)

    res += h1 * h2;
  }

  // also: doubletsq measures 0.5 Tr Phi^+ Phi, so multiply by 4 to get H(l)
  return 4.0*res / ((double) l->vol);
}

#endif // if NHIGGS > 0


#ifdef TRIPLET

/* Calculate (1/V) \sum_z sig(z) sig(z+l),
* where sig(z) = \sum_{x,y} Tr\Sigma(z)^2.
* Below the argument d = l and dir specifies which direction z lives in. */
double triplet_correlator(lattice* l, fields const* f, int d, int dir) {

  double res = 0.0;

  for (int z=0; z<l->L[dir]; z++) {

    // coordinate at distance d, accounting for periodicity
    int zd = (z + d) % l->L[dir];

    double s1 = plane_sum(tripletsq, f->su2triplet, l, z, dir);
    double s2 = plane_sum(tripletsq, f->su2triplet, l, zd, dir);

    res += s1 * s2;
  }

  return res / ((double) l->vol);

}

/* Calculate  sum_{x1,x2} Tr Sigma * F_{12}, which goes over
* to 2 / sqrt(beta) * \int d^2x Tr Sigma F_{12}
* in the cont. limit. Here F_{12} is calculated from the clover-improved plaquette.
* Should give better continuum limit. */
complex projected_photon_operator(lattice* l, fields const* f, int z, int dir, int* mom) {

  // momentum is quantized in all directions, (k1, k2, k3) = 2pi(n1/L1, n2/L2, n3/L3)
  // with n_i being integers. the input mom[i] = n_i

  complex res;
  res.re = 0.0; res.im = 0.0;

  // pre-calculate the sines and cosines??

  int dir1, dir2;
  if (dir == 0) {
    dir1 = 1; dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0; dir2 = 2;
  } else if (dir == 2) {
    dir1 = 0; dir2 = 1;
  }

  long z_node = z - l->offset[dir]; // coordinate on my node
  int skip = 0;

  if (z_node < 0 || z_node >= l->sliceL[dir]) {
    skip = 1; // out of bounds, so my node does not contribute
  }
  if (!skip) {
    for (long xy=0; xy<l->sites_per_coord[dir]; xy++) {
      long site = l->sites_at_coord[dir][z_node][xy];

      // calculate plaquette clover in dir1, dir2
      double clov[SU2LINK];
      clover_su2(l, f, site, dir1, dir2, clov);

      /* calculate field strength from clover O_munu:
      * g F_munu(x) = -i/8 [(O_munu - O^+_munu) - 1/N Tr(O_munu - O^+_munu) ],
      * but the trace vanishes for SU(2). */

      clov[0] = 0;
      for (int k=1; k<SU2LINK; k++) {
        clov[k] /= 4.0; // now clov = i g F_munu
      }

      double* a = f->su2triplet[site];
      // tr <- g Tr Sigma * F_munu, which is real:
      double tr = a[0]*clov[1] + a[1]*clov[2] + a[2]*clov[3];

      // multiply by exp(ip.k):
      double px = 0.0;
      for (int momdir=0; momdir<l->dim; momdir++) {
        if (mom[momdir] != 0) {
          px += ((double) (mom[momdir] * l->coords[site][momdir])) / ( (double)(l->L[momdir]) );
        }
      }
      /* note that blocked lattices use their own L[dir] and coords here.
      * For blocking level k, L_k = L / 2^k; x_k = x / 2^k,
      * so the ratio x / L is the same for all k. */

      res.re += tr * cos(2.0*M_PI * px);
      res.im += tr * sin(2.0*M_PI * px);
    }
  } // end else

  // normalize by area
  // double area = (double) (l->vol / l->L[dir]);

  res.re = allreduce(res.re, l->comm);
  res.im = allreduce(res.im, l->comm);
  return res;
}

/* Calculate correlation of the projected photon operator above, at distance d. */
complex projected_photon_correlator(lattice* l, fields const* f, int d, int dir) {

  complex res; res.re = 0.0; res.im = 0.0;

  // dunno how to generalize this to n. dimensions
  if (l->dim != 3 && l->rank) {
    printf("Warning: photon correlation function only sensible in 3 dimensions!!\n");
    return res;
  }

  // use first nonzero transverse momentum channel
  int mom[l->dim];
  if (dir == 0) {
    mom[1] = 1; mom[0] = 0; mom[2] = 0;
  } else if (dir == 1 || dir == 2) {
    mom[0] = 1; mom[1] = 0; mom[2] = 0;
  }

  for (int z=0; z<l->L[dir]; z++) {

    // coordinate at distance d, accounting for periodicity
    int zd = (z + d) % l->L[dir];

    complex op1 = projected_photon_operator(l, f, z, dir, mom);
    complex op2 = projected_photon_operator(l, f, zd, dir, mom);

    // O1.O2*
    op2.im = -1.0*op2.im;
    complex prod = cmult(op1, op2);
    res.re += prod.re;
    res.im += prod.im;
  }

  res.re /= (double) (l->vol);
  res.im /= (double) (l->vol);
  // imag part should be zero because the correlator is symmetric wrt. z <-> z+d

  return res;
}


#endif // end if TRIPLET


void print_labels_correlators() {

  int k = 1;

	FILE* f = fopen("labels_correl", "w");
  fprintf(f, "%d distance\n", k); k++;
  #if (NHIGGS > 0)
    fprintf(f, "%d H(l) = sum_z h(z) h(z+l) / V\n", k); k++; // Higgs correlator
  #endif
  #ifdef TRIPLET
    fprintf(f, "%d Sigma^2 correlator\n", k); k++;
    fprintf(f, "%d projected photon correlator RE\n", k); k++;
    fprintf(f, "%d projected photon correlator IM\n", k); k++;
  #endif

  fclose(f);

}

void measure_correlators(char* fname, lattice* l, fields const* f, params const* p, int dir, int meas_id) {

  FILE* file;

  if (!l->rank) {
		file = fopen(fname, "a");
    // write header for the current set of measurements
    fprintf(file, "\n =========== Measurement id: %d ===========\n", meas_id);
	}

  // no point measuring along the full length on a periodic lattice
  int max_distance = l->L[dir] / 2;

  for (int d=0; d<max_distance; d++) {
    #if (NHIGGS > 0)
      double hcorr = higgs_correlator(l, f, d, dir, 0);
    #endif
    #ifdef TRIPLET
      double tripcorr = triplet_correlator(l, f, d, dir);
      complex photoncorr = projected_photon_correlator(l, f, d, dir);
    #endif

    // write
    if (!l->rank) {
      fprintf(file, "%d ", d);
      #if (NHIGGS > 0)
        fprintf(file, "%g ", hcorr);
      #endif
      #ifdef TRIPLET
        fprintf(file, "%g ", tripcorr);
        fprintf(file, "%g %g ", photoncorr.re, photoncorr.im);
      #endif

      fprintf(file, "\n");
    }
  } // end d

  if (!l->rank) {
    fclose(file);
  }

}


#ifdef BLOCKING
/* Smear fields and calculate the correlators on a blocked lattice 'b' */
void measure_blocked_correlators(lattice* l, lattice* b, fields const* f, fields* f_b, params const* p,
      int const* block_dir, int dir, int id) {

  if (block_dir[dir] && !l->rank) {
    printf("Warning: calculating correlation lengths in a blocked direction\n");
  }
  fields f_smear;
  alloc_fields(l, &f_smear);

  #if (NHIGGS > 0)
    printf("Error in correlation.c: Higgs smearing not yet implemented!!!\n");
    free_fields(l, &f_smear);
    return;
  #endif

  smear_fields(l, f, &f_smear, block_dir);
  // transfer the smeared fields on the blocked lattice:
  make_blocked_fields(l, b, &f_smear, f_b);
  // then measure correlators along dir
  char fname[100];
  sprintf(fname, "correl%d", b->blocking_level);
  if (!b->standby) {
    measure_correlators(fname, b, f_b, p, dir, id);
  }

  free_fields(l, &f_smear);

}
#endif

#endif

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
* Note: does not divide by the number of sites, and no comms here. */
double plane_sum(double (*funct)(double*), double** field, lattice* l, long x, int dir) {

  double res = 0.0;
  long x_node = x - l->offset[dir]; // coordinate on my node

  if (x_node < 0 || x_node >= l->sliceL[dir]) {
    return 0.0; // out of bounds, so my node does not contribute
  } else {
    for (long i=0; i<l->sites_per_coord[dir]; i++) {
      long site = l->sites_at_coord[dir][x_node][i];
      res += (*funct)(field[site]);
    }
  }

  return res;
}

/* Calculate H(l) = (1/V) \sum_z h(z) h(z+l),
* where h(z) = \sum_{x,y} Tr\Phi(z)\Phi(z).
* Below the argument d = l and dir specifies which direction z lives in. */
double higgs_correlator(lattice* l, fields const* f, int d, int dir) {

  double res = 0.0;
  for (int z=0; z<l->L[dir]; z++) {

    /* periodicity: need to find site at coordinate z + d.
    * if x + d >= lattice length, we measure correlation in the "negative"
    * direction instead. */
    int zd = (z + d) % l->L[dir];

    double h1 = plane_sum(doubletsq, f->su2doublet, l, z, dir); // h1 <- h(z)
    double h2 = plane_sum(doubletsq, f->su2doublet, l, zd, dir); // h2 <- h(z+l)

    h1 = allreduce(h1, l->comm);
    h2 = allreduce(h2, l->comm);

    res += h1 * h2 / l->vol;
  }

  // also: doubletsq measures 0.5 Tr Phi^+ Phi, so multiply by 4 to get H(l)
  return 4.0*res;
}


#ifdef TRIPLET

/* Calculate (1/V) \sum_z sig(z) sig(z+l),
* where sig(z) = \sum_{x,y} Tr\Sigma(z)^2.
* Below the argument d = l and dir specifies which direction z lives in. */
double triplet_correlator(lattice* l, fields const* f, int d, int dir) {

  double res = 0.0;
  for (int z=0; z<l->L[dir]; z++) {

    // periodicity
    int zd = (z + d) % l->L[dir];

    double s1 = plane_sum(tripletsq, f->su2triplet, l, z, dir);
    double s2 = plane_sum(tripletsq, f->su2triplet, l, zd, dir);

    s1 = allreduce(s1, l->comm);
    s2 = allreduce(s2, l->comm);

    res += s1 * s2 / l->vol;
  }

  return res;
}

/* Calculate a "projected photon" operator in SU(2) + adjoint Higgs theory
* that is to be used for correlation measurements.  First calculates the projected
* field strength alpha_ij as specified in magfield.c, then O(z) = sum_{x} Im u_{1,2},
* where u_{1,2} is the U(1) plaquette and the x-sum is over the plane orthogonal to z-dir */
double projected_photon_operator_old(lattice* l, fields const* f, params const* p, int z, int dir,
      int* mom, double* res_re, double* res_im) {

  // momentum is quantized in all directions, (k1, k2, k3) = 2pi(n1/L1, n2/L2, n3/L3)
  // with n_i being integers. the input mom[i] = n_i

  // res[0] = real part, res[1] = imag part

  // dunno how to generalize this to n. dimensions
  if (l->dim != 3 && l->rank) {
    printf("Warning: photon correlation function only sensible in 3 dimensions!!\n");
    return 0.0;
  }

  double re = 0.0; double im = 0.0;

  // pre-calculate the sines and cosines??

  long z_node = z - l->offset[dir]; // coordinate on my node

  int dir1, dir2;
  if (dir == 0) {
    dir1 = 1; dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0; dir2 = 2;
  } else if (dir == 2) {
    dir1 = 0; dir2 = 1;
  }

  if (z_node < 0 || z_node >= l->sliceL[dir]) {
    return 0.0; // out of bounds, so my node does not contribute
  } else {
    for (long i=0; i<l->sites_per_coord[dir]; i++) {
      long site = l->sites_at_coord[dir][z_node][i];
      double alp = alpha_proj(l, f, p, site, dir1, dir2);

      /* alpha_ij = 2/g arg Tr p_ij, with p_ij = projected "plaquette" (2x2 matrix!)
      * For U(1) gauge theory we would have P_ij = e^{iga^2 F_ij}, so the above suggests
      * that the U(1) plaquette P_ij ~ Tr p_ij, so to get Im part of P_ij
      * we take sin(1/2 g alpha_ij) */
      double mod = sin(alp / sqrt(p->betasu2));

      // multiply by exp(ip.k):
      double px = 0.0;
      for (int momdir=0; momdir<l->dim; momdir++) {
        if (mom[momdir] != 0) {
          px += mom[momdir] * l->coords[site][momdir] / l->L[momdir];
        }
      }

      // real part
      re += mod * cos(2*M_PI * px);
      // imag part
      im += mod * sin(2*M_PI * px);
    }
  }

  *res_re = re;
  *res_im = im;
  // done
}

/* Calculate correlation of the projected photon operator above, at distance d. */
void projected_photon_correlator_old(lattice* l, fields const* f, params const* p, int d,
  int dir, double* res_re, double* res_im) {

  double re = 0.0; double im = 0.0;

  // use first nonzero transverse momentum channel
  int mom[l->dim];
  if (dir == 0) {
    mom[1] = 1; mom[0] = 0; mom[2] = 0;
  } else if (dir == 1 || dir == 2) {
    mom[0] = 1; mom[1] = 0; mom[2] = 0;
  }

  for (int z=0; z<l->L[dir]; z++) {

    double re1, im1;
    double re2, im2;

    // periodicity
    int zd = (z + d) % l->L[dir];

    projected_photon_operator_old(l, f, p, z, dir, mom, &re1, &im1);
    projected_photon_operator_old(l, f, p, zd, dir, mom, &re2, &im2);

    // combine from all nodes and normalize each operator by the "area"
    double area = (double) l->vol / l->L[dir];

    re1 = allreduce(re1, l->comm) / area;
    im1 = allreduce(im1, l->comm) / area;
    re2 = allreduce(re2, l->comm) / area;
    im2 = allreduce(im2, l->comm) / area;

    // O1 * Conjugate[O2], real and imag parts
    re += (re1*re2 + im1*im2);
    im += *res_im + im1*re2 - re1*im2;
  }

  *res_re = re;
  *res_im = im;
}

/* Calculate  -i * sum_{x1,x2} Tr Sigma * P_{12}, which goes over to 2 / sqrt(beta) * \int d^2x Tr Sigma F_{12}
* in the cont. limit */
double projected_photon_operator(lattice* l, fields const* f, int z, int dir,
      int* mom, double* res_re, double* res_im) {

  // momentum is quantized in all directions, (k1, k2, k3) = 2pi(n1/L1, n2/L2, n3/L3)
  // with n_i being integers. the input mom[i] = n_i

  // dunno how to generalize this to n. dimensions
  if (l->dim != 3 && l->rank) {
    printf("Warning: photon correlation function only sensible in 3 dimensions!!\n");
    return 0.0;
  }

  double re = 0.0; double im = 0.0;

  // pre-calculate the sines and cosines??

  long z_node = z - l->offset[dir]; // coordinate on my node

  int dir1, dir2;
  if (dir == 0) {
    dir1 = 1; dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0; dir2 = 2;
  } else if (dir == 2) {
    dir1 = 0; dir2 = 1;
  }

  if (z_node < 0 || z_node >= l->sliceL[dir]) {
    return 0.0; // out of bounds, so my node does not contribute

  } else {
    for (long i=0; i<l->sites_per_coord[dir]; i++) {
      long site = l->sites_at_coord[dir][z_node][i];

      // calculate plaquette in dir1, dir2
      double plaq[SU2LINK];
      su2plaquette(l, f, site, dir1, dir2, plaq);

      double* a = f->su2triplet[i];
      // tr <- -i Tr Sigma * plaq, which is real:
      double tr = a[0]*plaq[1] + a[1]*plaq[2] + a[2]*plaq[3];

      // multiply by exp(ip.k):
      double px = 0.0;
      for (int momdir=0; momdir<l->dim; momdir++) {
        if (mom[momdir] != 0) {
          px += mom[momdir] * l->coords[site][momdir] / l->L[momdir];
        }
      }
      /* note that blocked lattices use their own L[dir] and coords here.
      * For blocking level k, L_k = L / 2^k; x_k = x / 2^k,
      * so the ratio x / L is the same for all k. */

      // real part
      re += tr * cos(2*M_PI * px);
      // imag part
      im += tr * sin(2*M_PI * px); // for small momenta this is almost 0
    }
  } // end else

  *res_re = re;
  *res_im = im;
  // done
}

/* Calculate correlation of the projected photon operator above, at distance d. */
void projected_photon_correlator(lattice* l, fields const* f, int d, int dir, double* res_re, double* res_im) {

  double re = 0.0;
  double im = 0.0;
  // use first nonzero transverse momentum channel
  int mom[l->dim];
  if (dir == 0) {
    mom[1] = 1; mom[0] = 0; mom[2] = 0;
  } else if (dir == 1 || dir == 2) {
    mom[0] = 1; mom[1] = 0; mom[2] = 0;
  }

  for (int z=0; z<l->L[dir]; z++) {

    double re1, im1;
    double re2, im2;

    // periodicity
    int zd = (z + d) % l->L[dir];

    projected_photon_operator(l, f, z, dir, mom, &re1, &im1);
    projected_photon_operator(l, f, zd, dir, mom, &re2, &im2);

    // combine from all nodes and normalize each operator by the "area"
    double area = (double) l->vol / l->L[dir];

    re1 = allreduce(re1, l->comm) / area;
    im1 = allreduce(im1, l->comm) / area;
    re2 = allreduce(re2, l->comm) / area;
    im2 = allreduce(im2, l->comm) / area;



    // O1 * Conjugate[O2], re and imag parts
    re += re1*re2 + im1*im2;
    im += im1*re2 - re1*im2;
  }

  *res_re = re;
  *res_im = im;

}


#endif // end if TRIPLET


void print_labels_correlators() {

  int k = 1;

	FILE* f = fopen("labels_correlator", "w");
  fprintf(f, "%d distance\n", k); k++;
  #ifdef HIGGS
    fprintf(f, "%d H(l) = sum_z h(z) h(z+l) / V\n", k); k++; // Higgs correlator
  #endif
  #ifdef TRIPLET
    fprintf(f, "%d Sigma^2 correlator\n", k); k++;
    fprintf(f, "%d projected photon correlator RE\n", k); k++;
    fprintf(f, "%d projected photon correlator IM\n", k); k++;
    fprintf(f, "%d photon corr other method RE\n", k); k++;
    fprintf(f, "%d photon corr other method IM\n", k); k++;
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
    #ifdef HIGGS
      double hcorr = higgs_correlator(l, f, d, dir);
    #endif
    #ifdef TRIPLET
      double tripcorr = triplet_correlator(l, f, d, dir);
      double photon_re = 0.0; double photon_im = 0.0;
      projected_photon_correlator(l, f, d, dir, &photon_re, &photon_im);
      double old_re = 0.0; double old_im = 0.0;
      projected_photon_correlator_old(l, f, p, d, dir, &old_re, &old_im);
    #endif

    // write
    if (!l->rank) {
      fprintf(file, "%d ", d);
      #ifdef HIGGS
        fprintf(file, "%g ", hcorr);
      #endif
      #ifdef TRIPLET
        fprintf(file, "%g ", tripcorr);
        fprintf(file, "%g %g ", photon_re, photon_im);
        fprintf(file, "%g %g ", old_re, old_im);
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
  smear_fields(l, f, &f_smear, block_dir);
  // transfer the smeared fields on the blocked lattice:
  make_blocked_fields(l, b, &f_smear, f_b);
  // then measure correlators along dir
  char fname[100];
  sprintf(fname, "correlators_%d", b->blocking_level);
  if (!b->standby) {
    measure_correlators(fname, b, f_b, p, dir, id);
  }

  free_fields(l, &f_smear);

}
#endif

#endif

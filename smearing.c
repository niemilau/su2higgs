/** @file smearing.c
*
* Routines for performing smearing transformations on the fields.
* Used in conjuction with blocking routines in blocking.c.
*/

#include "su2.h"


/* Smearing routines. These construct a smeared field at each site
* by averaging over the field and its covariant connection with nearest neighbors.
* Naturally used with blocking routines.
* In all routines here the input
* smear_dir[j] = 1 if the direction j is to be smeared, 0 otherwise. */

#ifdef BLOCKING
/* NB! the BLOCKING flag here protects smear_link() in particular, which can fail
* in parallel implementation if the halo width is just 1 because it needs the
* staple calculated also at a neighbor site, which may be in halo. in su2.h
* HALOWIDTH is set to 2 if the BLOCKING flag is defined in makefile */


/* Smear the SU(2) link at site i and store in res (link components).
* This is done by calculating extended staples in the smearing directions,
* including only pure gauge contributions.
* See M. Teper, Phys.Lett.B 183 (1987) and hep-lat/9602006.
* However Kari seems to use slightly different blocking where he
* separately calculates two staples, and I don't think this is analytically
* the same as Teper's method. I have routines for both methods. */

/* Calculate smeared SU(2) link according to Eq. (15) in hep-lat/9602006,
* but normalize by the number of paths involved. Store in 'res', in SU(2) link parametrization.
* Note that this requires U_i(x+2i), so for MPI we need halo that extends two lattice spacings. */
void smear_link(lattice const* l, fields const* f, int const* smear_dir, double* res, long i, int dir) {

	if (!smear_dir[dir]) {
		printf("!!! WARNING: smearing gauge link without blocking the lattice dimension (in su2u1.c)\n");
	}

	int paths = 1;
	long next = l->next[i][dir];

	long x = i;
	memcpy(res, f->su2link[x][dir], SU2LINK * sizeof(*f->su2link[x][dir]));
	su2rot(res, f->su2link[next][dir], 0); // res <- Ui(x).Ui(x+i); now i = dir

	// Staple paths
	for (int j=0; j<l->dim; j++) {
		if (j != dir && smear_dir[j]) {
			double u[SU2LINK];
			// u <- Uj(x).Ui(x+j).Ui(x+i+j).Uj^+(x+2i)
			memcpy(u, f->su2link[x][j], SU2LINK * sizeof(*f->su2link[x][j]));
			next = l->next[x][j];
			su2rot(u, f->su2link[next][dir], 0);
			su2rot(u, f->su2link[ l->next[next][dir] ][dir], 0);
			next = l->next[ l->next[x][dir] ][dir];
			su2rot(u, f->su2link[next][j], 1);

			// upper staple done.
			for (int k=0; k<SU2LINK; k++) res[k] += u[k];

			// then the lower staple: u <- Uj^+(x-j).Ui(x-j).Ui(x+i-j).Uj(x+2i-j)
			long prev = l->prev[x][j];
			memcpy(u, f->su2link[prev][j], SU2LINK * sizeof(*f->su2link[prev][j]));
			// conjugate it:
			for (int k=1; k<SU2LINK; k++) u[k] *= -1.0;
			// then multiply with the other links
			su2rot(u, f->su2link[prev][dir], 0);
			next = l->next[prev][dir];
			su2rot(u, f->su2link[next][dir], 0);
			su2rot(u, f->su2link[ l->next[next][dir] ][j], 0);

			// lower staple done.
			for (int k=0; k<SU2LINK; k++) res[k] += u[k];

			paths += 2;
		}
	}

	// normalize and done.
	for (int k=0; k<SU2LINK; k++) res[k] /= paths;

}


/* Kari's method: eq. (6.5)-(6.6) in hep-lat/9612006.
* The smeared link is
* U_i(y) = V_i(x) V_i(x+i), where
*		V_i(x) = U_i(x) + \sum_j U_j(x) U_i(x+j) U_j^+(x+i)
* with the sum running over blocked directions, including backwards dirs, and j != i.
* Kari has normalization factor of 1/3 in V_i because for i=1,2, the original link
* plus a backwards and forward staple contribute in 3D. */
void smear_link_kari(lattice const* l, fields const* f, int const* smear_dir, double* res, long i, int dir) {

	if (!smear_dir[dir]) {
		printf("WARNING: smearing gauge link without blocking the lattice dimension (in su2u1.c)\n");
	}

	double stap[SU2LINK] = { 0.0 };
	memcpy(res, f->su2link[i][dir], SU2LINK * sizeof(*f->su2link[i][dir])); // will first construct res = V_dir(x)

	double v2[SU2LINK]; // v2 = V_dir(x+dir)
	memcpy(v2, f->su2link[ l->next[i][dir] ][dir], SU2LINK * sizeof(*f->su2link[i][dir]));

	int paths = 1; // how many "paths" are involved in smearing
	for (int j=0; j<l->dim; j++) {
		if (j != dir && smear_dir[j]) {
			// staples, but with the usual direction reversed (so take Hermitian conjugate)
			su2staple_wilson_onedir(l, f, i, dir, j, 1, stap);
			for (int k=0; k<SU2LINK; k++) {
				res[k] += stap[k];
			}

			// repeat for site x+dir
			su2staple_wilson_onedir(l, f, l->next[i][dir], dir, j, 1, stap);
			for (int k=0; k<SU2LINK; k++) {
				v2[k] += stap[k];
			}

			paths += 2; // counterwise and clockwise staples
		}
	}

	// normalize
	for (int k=0; k<SU2LINK; k++) {
		res[k] /= ((double) paths);
		v2[k] /= ((double) paths);
	}

	su2rot(res, v2, 0); // res <- V1.V2 = V_dir(x) V_dir(x+dir)
	// in general this is unitary but not necessarily in SU(2). so normalize again:
	double det = su2sqr(res);
	det = sqrt(det);
	for (int k=0; k<SU2LINK; k++) {
		res[k] /= det;
	}

}


// TODO doublet
void smear_doublet(lattice const* l, fields const* f, int const* smear_dir, double* res, long i) {
}

#ifdef SINGLET

/* Smear a gauge singlet field at site i and store in res[0].
* The smearing transformation is
* 	S_smear(x) = 1/N (S(x) + sum_j S(x+i))
* where N = number of sites involved, and j runs over positive and negative directions. */
double smear_singlet(lattice const* l, fields const* f, int const* smear_dir, double* res, long i) {

	int sites = 1; // how many sites are involved in smearing
	res[0] = f->singlet[i][0];

	for (int dir=0; dir<l->dim; dir++) {
		if (smear_dir[dir]) {

			long next = l->next[i][dir];
			long prev = l->prev[i][dir];

			res[0] += f->singlet[next][0];
			res[0] += f->singlet[prev][0];
			sites += 2; // involves 2 nearest neighbors
		}
	}
	res[0] /= sites;

}

#endif

/* Smear the triplet field at site i and store in res (triplet components).
* This is done by calculating
* 	Sigma(x) + sum_j U_j(x) Sigma(x+j) U^+_j (x)
* and normalizing with the number of sites involved. The sum is over
* directions specified in smear_dir, and contains both forward and backward terms.
* Here U_{-j}(x) = U^+_j(x-j). */
void smear_triplet(lattice const* l, fields const* f, int const* smear_dir, double* res, long i) {

	int sites = 1; // how many sites are involved in smearing
	double cov[SU2TRIP] = {0.0};

	for (int dir=0; dir<l->dim; dir++) {
		if (smear_dir[dir]) {
			// forward connection Ui(x) Sigma(x+i) Ui^+(x)
			double* u = f->su2link[i][dir];
			long next = l->next[i][dir];
			double* b = f->su2triplet[next];
			cov[0] += b[0]*(u[0]*u[0]) + b[0]*(u[1]*u[1]) - 2*b[2]*u[0]*u[2]
							+ 2*b[1]*u[1]*u[2] - b[0]*(u[2]*u[2]) + 2*b[1]*u[0]*u[3]
							+ 2*b[2]*u[1]*u[3] - b[0]*(u[3]*u[3]);

			cov[1] += b[1]*(u[0]*u[0]) + 2*b[2]*u[0]*u[1] - b[1]*(u[1]*u[1])
							+ 2*b[0]*u[1]*u[2] + b[1]*(u[2]*u[2]) - 2*b[0]*u[0]*u[3]
							+ 2*b[2]*u[2]*u[3] - b[1]*(u[3]*u[3]);

			cov[2] += b[2]*(u[0]*u[0]) - 2*b[1]*u[0]*u[1] - b[2]*(u[1]*u[1])
							+ 2*b[0]*u[0]*u[2] - b[2]*(u[2]*u[2]) + 2*b[0]*u[1]*u[3]
							+ 2*b[1]*u[2]*u[3] + b[2]*(u[3]*u[3]);

			// backward connection Ui^+(x-i) Sigma(x-i) Ui(x-i)
			long prev = l->prev[i][dir];
			u = f->su2link[prev][dir];
			b = f->su2triplet[prev];
			cov[0] += b[0]*(u[0]*u[0]) + b[0]*(u[1]*u[1]) + 2*b[2]*u[0]*u[2]
							+ 2*b[1]*u[1]*u[2] - b[0]*(u[2]*u[2]) - 2*b[1]*u[0]*u[3]
							+ 2*b[2]*u[1]*u[3] - b[0]*(u[3]*u[3]);

			cov[1] += b[1]*(u[0]*u[0]) - 2*b[2]*u[0]*u[1] - b[1]*(u[1]*u[1])
							+ 2*b[0]*u[1]*u[2] + b[1]*(u[2]*u[2]) + 2*b[0]*u[0]*u[3]
							+ 2*b[2]*u[2]*u[3] - b[1]*(u[3]*u[3]);

			cov[2] += b[2]*(u[0]*u[0]) + 2*b[1]*u[0]*u[1] - b[2]*(u[1]*u[1])
							- 2*b[0]*u[0]*u[2] - b[2]*(u[2]*u[2]) + 2*b[0]*u[1]*u[3]
							+ 2*b[1]*u[2]*u[3] + b[2]*(u[3]*u[3]);

			sites += 2; // involves 2 nearest neighbors
		}
	}

	for (int k=0; k<SU2TRIP; k++) {
		res[k] = f->su2triplet[i][k] + cov[k];
		res[k] = res[k] / ((double) sites);
	}
}


/* Smear all fields and store in f_b. Does not smear fields at odd sites in smearing
* directions, because those are not needed for blocked lattices. */
void smear_fields(lattice const* l, fields const* f, fields* f_b, int const* block_dir) {

  // no halos
  for (long i=0; i<l->sites; i++) {
    int skip = 0;

    // smear fields only on the sites that end up on the blocked lattice
    for (int dir=0; dir<l->dim; dir++) {
      if (block_dir[dir] && l->coords[i][dir] % 2 != 0) {
        // no need to smear these
        skip = 1;
        break;
      }
    }

    if (!skip) {
      // create blocked fields by smearing in the specified dirs.

      // gauge links. links pointing in non-smeared directions remain unchanged
      for (int dir=0; dir<l->dim; dir++) {
        if (block_dir[dir]) {
          smear_link(l, f, block_dir, f_b->su2link[i][dir], i, dir);
          #ifdef U1

          #endif
        } else {
          memcpy(f_b->su2link[i][dir], f->su2link[i][dir], SU2LINK * sizeof(*f->su2link[i][dir]));
          #ifdef U1
            f_b->u1link[i][dir] = f->u1link[i][dir];
          #endif
        }
      } // end dir

      // scalars
      #if (NHIGGS > 0)

      #endif
			#ifdef SINGLET
				smear_singlet(l, f, block_dir, f_b->singlet[i], i);
			#endif
      #ifdef TRIPLET
        smear_triplet(l, f, block_dir, f_b->su2triplet[i], i);
      #endif

    }

  } // end i
}

#endif // ifdef BLOCKING

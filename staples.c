/** @file staples.c
*
* Routines for calculating local "staples" for different fields.
* By staple I mean the part of the action that depends linearly on a local field,
* i.e.
* 	S_x = f_a(x) * s_a(x) + O(f^2) + f-independent
* where a labels the components of the field at site x. s_a is the staple.
*/

#include "su2.h"

/* ----- SU(2) gauge links ----- */

/* Calculate untraced staple in counterclockwise direction ("usual" direction),
*  and store it in V.
* Specifically, calculates:
* 	U1 U2^+ U3^+ */
void su2staple_counterwise(double* V, double* u1, double* u2, double* u3) {

	V[0] = u1[0]*u2[0]*u3[0] + u1[1]*u2[1]*u3[0] + u1[2]*u2[2]*u3[0] + u1[3]*u2[3]*u3[0] +
   u1[1]*u2[0]*u3[1] - u1[0]*u2[1]*u3[1] - u1[3]*u2[2]*u3[1] +
   u1[2]*u2[3]*u3[1] + u1[2]*u2[0]*u3[2] + u1[3]*u2[1]*u3[2] -
   u1[0]*u2[2]*u3[2] - u1[1]*u2[3]*u3[2] + u1[3]*u2[0]*u3[3] -
   u1[2]*u2[1]*u3[3] + u1[1]*u2[2]*u3[3] - u1[0]*u2[3]*u3[3];

	 V[1] = u1[1]*u2[0]*u3[0] - u1[0]*u2[1]*u3[0] - u1[3]*u2[2]*u3[0] + u1[2]*u2[3]*u3[0] -
   u1[0]*u2[0]*u3[1] - u1[1]*u2[1]*u3[1] - u1[2]*u2[2]*u3[1] -
   u1[3]*u2[3]*u3[1] - u1[3]*u2[0]*u3[2] + u1[2]*u2[1]*u3[2] -
   u1[1]*u2[2]*u3[2] + u1[0]*u2[3]*u3[2] + u1[2]*u2[0]*u3[3] +
   u1[3]*u2[1]*u3[3] - u1[0]*u2[2]*u3[3] - u1[1]*u2[3]*u3[3];

	 V[2] = u1[2]*u2[0]*u3[0] + u1[3]*u2[1]*u3[0] - u1[0]*u2[2]*u3[0] - u1[1]*u2[3]*u3[0] +
   u1[3]*u2[0]*u3[1] - u1[2]*u2[1]*u3[1] + u1[1]*u2[2]*u3[1] -
   u1[0]*u2[3]*u3[1] - u1[0]*u2[0]*u3[2] - u1[1]*u2[1]*u3[2] -
   u1[2]*u2[2]*u3[2] - u1[3]*u2[3]*u3[2] - u1[1]*u2[0]*u3[3] +
   u1[0]*u2[1]*u3[3] + u1[3]*u2[2]*u3[3] - u1[2]*u2[3]*u3[3];

	 V[3] = u1[3]*u2[0]*u3[0] - u1[2]*u2[1]*u3[0] + u1[1]*u2[2]*u3[0] - u1[0]*u2[3]*u3[0] -
   u1[2]*u2[0]*u3[1] - u1[3]*u2[1]*u3[1] + u1[0]*u2[2]*u3[1] +
   u1[1]*u2[3]*u3[1] + u1[1]*u2[0]*u3[2] - u1[0]*u2[1]*u3[2] -
   u1[3]*u2[2]*u3[2] + u1[2]*u2[3]*u3[2] - u1[0]*u2[0]*u3[3] -
   u1[1]*u2[1]*u3[3] - u1[2]*u2[2]*u3[3] - u1[3]*u2[3]*u3[3];
}

/* Calculate untraced staple in clockwise direction ("inverse" direction),
*  and store it in V.
* Specifically, calculates:
* 	U1^+ U2^+ U3
*/
void su2staple_clockwise(double* V, double* u1, double* u2, double* u3) {

	V[0] = u1[0]*u2[0]*u3[0] - u1[1]*u2[1]*u3[0] - u1[2]*u2[2]*u3[0] - u1[3]*u2[3]*u3[0] +
   u1[1]*u2[0]*u3[1] + u1[0]*u2[1]*u3[1] - u1[3]*u2[2]*u3[1] +
   u1[2]*u2[3]*u3[1] + u1[2]*u2[0]*u3[2] + u1[3]*u2[1]*u3[2] +
   u1[0]*u2[2]*u3[2] - u1[1]*u2[3]*u3[2] + u1[3]*u2[0]*u3[3] -
   u1[2]*u2[1]*u3[3] + u1[1]*u2[2]*u3[3] + u1[0]*u2[3]*u3[3];

	 V[1] = -u1[1]*u2[0]*u3[0] - u1[0]*u2[1]*u3[0] + u1[3]*u2[2]*u3[0] -
   u1[2]*u2[3]*u3[0] + u1[0]*u2[0]*u3[1] - u1[1]*u2[1]*u3[1] -
   u1[2]*u2[2]*u3[1] - u1[3]*u2[3]*u3[1] - u1[3]*u2[0]*u3[2] +
   u1[2]*u2[1]*u3[2] - u1[1]*u2[2]*u3[2] - u1[0]*u2[3]*u3[2] +
   u1[2]*u2[0]*u3[3] + u1[3]*u2[1]*u3[3] + u1[0]*u2[2]*u3[3] - u1[1]*u2[3]*u3[3];

	 V[2] = -(u1[2]*u2[0]*u3[0]) - u1[3]*u2[1]*u3[0] - u1[0]*u2[2]*u3[0] +
   u1[1]*u2[3]*u3[0] + u1[3]*u2[0]*u3[1] - u1[2]*u2[1]*u3[1] +
   u1[1]*u2[2]*u3[1] + u1[0]*u2[3]*u3[1] + u1[0]*u2[0]*u3[2] -
   u1[1]*u2[1]*u3[2] - u1[2]*u2[2]*u3[2] - u1[3]*u2[3]*u3[2] -
   u1[1]*u2[0]*u3[3] - u1[0]*u2[1]*u3[3] + u1[3]*u2[2]*u3[3] - u1[2]*u2[3]*u3[3];

	 V[3] = -(u1[3]*u2[0]*u3[0]) + u1[2]*u2[1]*u3[0] - u1[1]*u2[2]*u3[0] -
   u1[0]*u2[3]*u3[0] - u1[2]*u2[0]*u3[1] - u1[3]*u2[1]*u3[1] -
   u1[0]*u2[2]*u3[1] + u1[1]*u2[3]*u3[1] + u1[1]*u2[0]*u3[2] +
   u1[0]*u2[1]*u3[2] - u1[3]*u2[2]*u3[2] + u1[2]*u2[3]*u3[2] +
   u1[0]*u2[0]*u3[3] - u1[1]*u2[1]*u3[3] - u1[2]*u2[2]*u3[3] - u1[3]*u2[3]*u3[3];
}


/* Calculate full untraced staple for an SU(2) link and stores it in V
* Only the contribution from Wilson action is included here, modulo the beta prefactor(s).
* Specifically, calculates:
* 	\sum_{nu != mu} (U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
 																	+ U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu) )
* where mu = dir. */
void su2staple_wilson(lattice const* l, fields const* f, long i, int dir, double* V) {
	double tot[SU2LINK] = { 0.0 };
	double* u1 = NULL;
	double* u2 = NULL;
	double* u3 = NULL;

	for (int j=0; j<l->dim; j++) {
		if (j != dir) {
			// "upper" staple
			u1 = f->su2link[ l->next[i][dir] ][j];
			u2 = f->su2link[ l->next[i][j] ][dir];
			u3 = f->su2link[i][j];
			su2staple_counterwise(V, u1, u2, u3);
			for(int k=0; k<SU2LINK; k++){
				tot[k] += V[k];
			}
			// "lower" staple
			u1 = f->su2link[ l->prev[(l->next[i][dir])][j] ][j];
			u2 = f->su2link[(l->prev[i][j])][dir];
			u3 = f->su2link[(l->prev[i][j])][j];
			su2staple_clockwise(V, u1, u2, u3);;
			for(int k=0; k<SU2LINK; k++){
				tot[k] += V[k];
			}
		}
	}

	for(int k=0; k<SU2LINK; k++){
		V[k] = tot[k];
	}

}


/* Calculate total untraced staple for a SU(2) link and store in V.
* Staple S here refers to the matrix multiplying link U in the trace, action ~ Tr US.
*	Specifically, calculates:
* 	-\beta * 0.5 * su2staple_wilson()
* plus whatever comes from other terms in the action.
* Wilson staple is
*    -----
*   |     |
*   V     |
*   ------->
*   ^     |
*   |     |
*    -----
*
* Triplet is not included here, because its hopping term is quadratic in the link
* and cannot be expressed as a simple staple. */
void su2link_staple(lattice const* l, fields const* f, params const* p, long i, int dir, double* V) {

	su2staple_wilson(l, f, i, dir, V);
	for (int k=0; k<4; k++) {
		V[k] *= -0.5 * p->betasu2;
	}

	#if (NHIGGS > 0)

		// Higgs doublet hopping terms: -Tr U_j Phi(x+j) exp(-i a_j(x) sigma_3) Phi(x)^+,
		// a_j(x) = 0 if hypercharge is neglected.
		double** phi;
		long nextsite = l->next[i][dir];
		double nextphi[SU2DB];
		double currentphi[SU2DB];

		for (int db=0; db<NHIGGS; db++) {

			phi = f->su2doublet[db];
			// we want Hermitian conjugate of Phi(x):
			for (int d=0; d<SU2DB; d++) {
				currentphi[d] = phi[i][d];
				if (d > 0) currentphi[d] *= -1.0;
			}

			memcpy(nextphi, phi[nextsite], SU2DB * sizeof(*nextphi));

			#ifdef U1
				// the U(1) contribution can be written as
				// I cos(a) - i sin(a) sigma_3, so in our notation it can be
				// treated as a doublet field with components
				// a[0] = sqrt(2) cos(a), a[1] = 0, a[2] = 0, a[3] = -sqrt(2) sin(a).
				// I assume Higgs hypercharge Y=1.
				double s = sin(f->u1link[i][dir]);
				double c = cos(f->u1link[i][dir]);
				double b[4];
				for (int k=0; k<SU2DB; k++) b[k] = nextphi[k];

				// calculate Phi(x+j) times the hypercharge bit
				nextphi[0] = c*b[0] + s*b[3];
				nextphi[1] = c*b[1] + s*b[2];
				nextphi[2] = -s*b[1] + c*b[2];
				nextphi[3] = -s*b[0] + c*b[3];
			#endif

			/* Here we use su2rot to do the multiplication.
			 However, we need a factor 1/sqrt(2) becase the doublet components
			 are normalized differently from SU(2) links. For the same reason
			 there is another 1/sqrt(2) when we add the scalar staple to the
			 link staple, so overall we add to the Wilson staple
			 -0.5 times what su2rot() of two doublets gives us. */
			su2rot(nextphi, currentphi, 0);
			for (int k=0; k<SU2LINK; k++) {
				V[k] -= 0.5 * nextphi[k];
			}
		} // doublet hoppings done

	#endif

}

/* Same as su2staple_wilson(), but only does the staple for U_mu(x) in one direction = nu.
* If dagger = 1, also takes the Hermitian conjugate of both forward and backward staples */
void su2staple_wilson_onedir(lattice const* l, fields const* f, long i, int mu, int nu, int dagger, double* res) {
	if (mu == nu) {
		res[0] = 1.0;
		for(int k=1; k<SU2LINK; k++) {
			res[k] = 0.0;
		}
		return;
	}

	double tot[SU2LINK] = { 0.0 };
	double* u1 = NULL;
	double* u2 = NULL;
	double* u3 = NULL;

	// "upper" staple U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
	u1 = f->su2link[ l->next[i][mu] ][nu];
	u2 = f->su2link[ l->next[i][nu] ][mu];
	u3 = f->su2link[i][nu];
	su2staple_counterwise(tot, u1, u2, u3);
	for(int k=0; k<SU2LINK; k++) {
		// take conjugate if needed
		if (dagger && k != 0) tot[k] = -1.0 * tot[k];
		res[k] = tot[k];
	}

	// "lower" staple U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
	long site = l->next[i][mu];
	site = l->prev[site][nu];
	u1 = f->su2link[site][nu];
	u2 = f->su2link[ l->prev[i][nu] ][mu];
	u3 = f->su2link[ l->prev[i][nu] ][nu];
	su2staple_clockwise(tot, u1, u2, u3);;
	for(int k=0; k<SU2LINK; k++) {
		// take conjugate if needed
		if (dagger && k != 0) tot[k] = -1.0 * tot[k];
		res[k] += tot[k];
	}
}

#ifdef U1
/* ----- U(1) links ----- */

/* Same as su2staple_wilson_onedir(), but for U(1) links
* Remember: my u1link[dir] is actually the phase. */
complex u1staple_wilson_onedir(lattice const* l, fields const* f, long i, int mu, int nu, int dagger) {

	complex res;
	res.re = 1.0; res.im = 0.0;

	if (mu == nu) {
		return res;
	}

	// "upper" staple U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
	double a1 = f->u1link[ l->next[i][mu] ][nu];
	double a2 = f->u1link[ l->next[i][nu] ][mu];
	double a3 = f->u1link[i][nu];

	double staple_phase = a1 - a2 - a3;
	if (dagger) staple_phase *= -1.0;

	res.re = cos(staple_phase);
	res.im = sin(staple_phase);

	// "lower" staple U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
	long site = l->next[i][mu];
	site = l->prev[site][nu];
	a1 = f->u1link[site][nu];
	a2 = f->u1link[ l->prev[i][nu] ][mu];
	a3 = f->u1link[ l->prev[i][nu] ][nu];

	staple_phase = -1.0*a1 - a2 + a3;
	if (dagger) staple_phase *= -1.0;

	res.re += cos(staple_phase);
	res.im += sin(staple_phase);

	return res;
}

#endif

/* ----- SU(2) doublets ----- */

#if (NHIGGS > 0)
/* Phi(x) = 1/sqrt(2) (f[0] I + i f[i]*sig[i]),
* --> calculate s[a] so that S[Phi(x)] = f[a] s[a] + O(f^2)
* but include only hopping terms (there can be other, local additions from the potential).
* Hopping is -\sum_i Re Tr Phi^+ U_i Phi(x+i) e^{-i Y alpha_i(x) sigma_3}, and minus sign is included here. */
void staple_doublet(double* res, lattice const* l, fields const* f, params const* p, long i, int higgs_id) {

	for (int k=0; k<SU2DB; k++) res[k] = 0.0;

	double** phi = f->su2doublet[higgs_id];
	double* u, *b;

	// hopping terms
	for (int dir=0; dir<l->dim; dir++) {

		u = f->su2link[i][dir];
		long next = l->next[i][dir];
		b = phi[next];

		// Forward direction (these were obtained in Mathematica):
		#ifndef U1
			res[0] += -(b[0]*u[0]) + b[1]*u[1] + b[2]*u[2] + b[3]*u[3];
			res[1] += -(b[1]*u[0]) - b[0]*u[1] + b[3]*u[2] - b[2]*u[3];
			res[2] += -(b[2]*u[0]) - b[3]*u[1] - b[0]*u[2] + b[1]*u[3];
			res[3] += -(b[3]*u[0]) + b[2]*u[1] - b[1]*u[2] - b[0]*u[3];
		#else
			// hypercharge Y = 1
			double ss = sin(f->u1link[i][dir]);
			double cc = cos(f->u1link[i][dir]);
			res[0] += -(cc*b[0]*u[0]) - ss*b[3]*u[0] + cc*b[1]*u[1] + ss*b[2]*u[1]
							- ss*b[1]*u[2] + cc*b[2]*u[2] - ss*b[0]*u[3] + cc*b[3]*u[3];
			res[1] += -(cc*b[1]*u[0]) - ss*b[2]*u[0] - cc*b[0]*u[1] - ss*b[3]*u[1]
							- ss*b[0]*u[2] + cc*b[3]*u[2] + ss*b[1]*u[3] - cc*b[2]*u[3];
			res[2] += ss*b[1]*u[0] - cc*b[2]*u[0] + ss*b[0]*u[1] - cc*b[3]*u[1]
							- cc*b[0]*u[2] - ss*b[3]*u[2] + cc*b[1]*u[3] + ss*b[2]*u[3];
			res[3] += ss*b[0]*u[0] - cc*b[3]*u[0] - ss*b[1]*u[1] + cc*b[2]*u[1]
							- cc*b[1]*u[2] - ss*b[2]*u[2] - cc*b[0]*u[3] - ss*b[3]*u[3];
		#endif

		// same for backwards directions
		long prev = l->prev[i][dir];
		u = f->su2link[prev][dir];
		b = phi[prev];

		#ifndef U1
			res[0] += -(b[0]*u[0]) - b[1]*u[1] - b[2]*u[2] - b[3]*u[3];
			res[1] += -(b[1]*u[0]) + b[0]*u[1] - b[3]*u[2] + b[2]*u[3];
			res[2] += -(b[2]*u[0]) + b[3]*u[1] + b[0]*u[2] - b[1]*u[3];
			res[3] += -(b[3]*u[0]) - b[2]*u[1] + b[1]*u[2] + b[0]*u[3];
		#else
			ss = sin(f->u1link[prev][dir]);
			cc = cos(f->u1link[prev][dir]);
			res[0] += -(cc*b[0]*u[0]) + ss*b[3]*u[0] - cc*b[1]*u[1] + ss*b[2]*u[1]
							- ss*b[1]*u[2] - cc*b[2]*u[2] - ss*b[0]*u[3] - cc*b[3]*u[3];
			res[1] += -(cc*b[1]*u[0]) + ss*b[2]*u[0] + cc*b[0]*u[1] - ss*b[3]*u[1]
							- ss*b[0]*u[2] - cc*b[3]*u[2] + ss*b[1]*u[3] + cc*b[2]*u[3];
			res[2] += -(ss*b[1]*u[0]) - cc*b[2]*u[0] + ss*b[0]*u[1] + cc*b[3]*u[1]
							+ cc*b[0]*u[2] - ss*b[3]*u[2] - cc*b[1]*u[3] + ss*b[2]*u[3];
			res[3] += -(ss*b[0]*u[0]) - cc*b[3]*u[0] - ss*b[1]*u[1] - cc*b[2]*u[1]
							+ cc*b[1]*u[2] - ss*b[2]*u[2] + cc*b[0]*u[3] - ss*b[3]*u[3];
		#endif

	} // hoppings done

}

#endif // if NHIGGS > 0

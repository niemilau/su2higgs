/** @file su2u1.c
*
* Routines for operating on SU(2) and U(1) fields and calculating the action.
* U(1) sector can be turned off by leaving out the U1 flag in makefile.
* Overall conventions are taken from hep-lat/9504001. Specifically:
*
*	Everywhere in this code, the SU(2) links are parameterized with real numbers as
* 	U = u0 I + i ( u1 sig1 + u2 sig2 + u3 sig3 ),
* where: i is the imaginary unit, I is the unit matrix, sig are the Pauli matrices.
*
*	The condition for U to be in SU(2) is
*		det U = u0^2 + u1^2 + u2^3 + u3^3 = 1.
*
* The component u_a is stored in fields.su2link[i][dir][a].
* Here i is the lattice site index (i = 0, ... vol-1) and dir is the direction.
*
*	In terms of the adjoint gauge fields A_i, the link is
*		U_i(x) = exp(0.5*i*g * sig_a * A^a_i(x)), with g = gauge coupling in lattice units.
*
*	SU(2) doublet scalars are parametrized as:
* 	\Phi = 1/sqrt(2) ( a0 I + i ( a1 sig1 + a2 sig2 + a3 sig3 ) ),
* We have then that Tr \Phi^+ \Phi = 2*(a[i]^2), i=0,1,2,3.
*	Higgs potential is \frac12 m^2 Tr \Phi^+ \Phi + \frac14 \lambda (Tr \Phi^+ \Phi)^2.
*
*	The normalizations of the ai components are now the same as usually in continuum.
* See also the Mathematica notebook su2.nb.
*
*
*	Real SU(2) triplets are parametrized as:
*		A = \frac12 A[a] sigma[a],  a = 1,2,3.
* This matches ref. hep-ph/9704416 but is different from that of hep-lat/9504001,
* which uses A = A[a] sigma[a] and has the action written in a slightly different,
* but equal, form.
*
*
*	My potential is  V = m^2 Tr A^2 + b4 (Tr A^2)^2 + 0.5 a2 Tr(\Phi^+ \Phi) Tr A^2,
*	where the normalization matches what is usually used in continuum for the A[a].
*
*
*	U(1) links are written simply as U_j(x) = exp(i a_j(x)), where a_j(x) is real
* and stored in u1link[x][j]. We assume "standard" Wilson action, i.e.
* 	gamma = 1, see doc.
//TODO
*
*/

#include "su2.h"


// calculate the determinant, or norm, of a matrix in our SU(2) parametrization
double su2sqr(double *u) {
	return u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
}


// rotate U1 from the right by U2 and store again in U1: U1 <- U1.U2.
void su2rot(double *u1, double *u2) {

	double new1 = u1[0]*u2[0] - u1[1]*u2[1] - u1[2]*u2[2] - u1[3]*u2[3];
	double new2 = u1[1]*u2[0] + u1[0]*u2[1] + u1[3]*u2[2] - u1[2]*u2[3];
	double new3 = u1[2]*u2[0] - u1[3]*u2[1] + u1[0]*u2[2] + u1[1]*u2[3];
	double new4 = u1[3]*u2[0] + u1[2]*u2[1] - u1[1]*u2[2] + u1[0]*u2[3];
	u1[0] = new1;
	u1[1] = new2;
	u1[2] = new3;
	u1[3] = new4;
}


/* Calculate trace of four SU(2) matrices as used in Wilson action. Returns:
*  	Re Tr U1.U2.U3^+.U4^+
*	Note that for SU(2), the trace is always real.
*/
double su2trace4(double *u1, double *u2, double *u3, double *u4) {

	return 2.0 * (u1[0]*u2[0]*u3[0]*u4[0] - u1[1]*u2[1]*u3[0]*u4[0] - u1[2]*u2[2]*u3[0]*u4[0] -
   u1[3]*u2[3]*u3[0]*u4[0] + u1[1]*u2[0]*u3[1]*u4[0] + u1[0]*u2[1]*u3[1]*u4[0] +
   u1[3]*u2[2]*u3[1]*u4[0] - u1[2]*u2[3]*u3[1]*u4[0] + u1[2]*u2[0]*u3[2]*u4[0] -
   u1[3]*u2[1]*u3[2]*u4[0] + u1[0]*u2[2]*u3[2]*u4[0] + u1[1]*u2[3]*u3[2]*u4[0] +
   u1[3]*u2[0]*u3[3]*u4[0] + u1[2]*u2[1]*u3[3]*u4[0] - u1[1]*u2[2]*u3[3]*u4[0] +
   u1[0]*u2[3]*u3[3]*u4[0] + u1[1]*u2[0]*u3[0]*u4[1] + u1[0]*u2[1]*u3[0]*u4[1] +
   u1[3]*u2[2]*u3[0]*u4[1] - u1[2]*u2[3]*u3[0]*u4[1] - u1[0]*u2[0]*u3[1]*u4[1] +
   u1[1]*u2[1]*u3[1]*u4[1] + u1[2]*u2[2]*u3[1]*u4[1] + u1[3]*u2[3]*u3[1]*u4[1] -
   u1[3]*u2[0]*u3[2]*u4[1] - u1[2]*u2[1]*u3[2]*u4[1] + u1[1]*u2[2]*u3[2]*u4[1] -
   u1[0]*u2[3]*u3[2]*u4[1] + u1[2]*u2[0]*u3[3]*u4[1] - u1[3]*u2[1]*u3[3]*u4[1] +
   u1[0]*u2[2]*u3[3]*u4[1] + u1[1]*u2[3]*u3[3]*u4[1] + u1[2]*u2[0]*u3[0]*u4[2] -
   u1[3]*u2[1]*u3[0]*u4[2] + u1[0]*u2[2]*u3[0]*u4[2] + u1[1]*u2[3]*u3[0]*u4[2] +
   u1[3]*u2[0]*u3[1]*u4[2] + u1[2]*u2[1]*u3[1]*u4[2] - u1[1]*u2[2]*u3[1]*u4[2] +
   u1[0]*u2[3]*u3[1]*u4[2] - u1[0]*u2[0]*u3[2]*u4[2] + u1[1]*u2[1]*u3[2]*u4[2] +
   u1[2]*u2[2]*u3[2]*u4[2] + u1[3]*u2[3]*u3[2]*u4[2] - u1[1]*u2[0]*u3[3]*u4[2] -
   u1[0]*u2[1]*u3[3]*u4[2] - u1[3]*u2[2]*u3[3]*u4[2] + u1[2]*u2[3]*u3[3]*u4[2] +
   u1[3]*u2[0]*u3[0]*u4[3] + u1[2]*u2[1]*u3[0]*u4[3] - u1[1]*u2[2]*u3[0]*u4[3] +
   u1[0]*u2[3]*u3[0]*u4[3] - u1[2]*u2[0]*u3[1]*u4[3] + u1[3]*u2[1]*u3[1]*u4[3] -
   u1[0]*u2[2]*u3[1]*u4[3] - u1[1]*u2[3]*u3[1]*u4[3] + u1[1]*u2[0]*u3[2]*u4[3] +
   u1[0]*u2[1]*u3[2]*u4[3] + u1[3]*u2[2]*u3[2]*u4[3] - u1[2]*u2[3]*u3[2]*u4[3] -
   u1[0]*u2[0]*u3[3]*u4[3] + u1[1]*u2[1]*u3[3]*u4[3] + u1[2]*u2[2]*u3[3]*u4[3] +
   u1[3]*u2[3]*u3[3]*u4[3]);

}

/* Calculate plaquette trace in the (dir1, dir2) plane
*	of SU(2) link at lattice site i. Returns:
* 	Re Tr U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
*	where mu = dir1, nu = dir2, x = site at index i */
double su2ptrace(lattice const* l, fields const* f, long i, int dir1, int dir2) {

	double *u1 = f->su2link[i][dir1];
	double *u2 = f->su2link[ l->next[i][dir1] ][dir2];
	double *u3 = f->su2link[ l->next[i][dir2] ][dir1];
	double *u4 = f->su2link[i][dir2];

	return su2trace4(u1, u2, u3, u4);
}

/* Calculate the SU(2) plaquette and store in u1 */
void su2plaquette(lattice const* l, fields const* f, long i, int dir1, int dir2, double* u1) {

	double u2[SU2LINK], u3[SU2LINK], u4[SU2LINK];
	memcpy(u1, f->su2link[i][dir1], SU2LINK * sizeof(*u1));
	memcpy(u2, f->su2link[ l->next[i][dir1] ][dir2], SU2LINK * sizeof(*u2));
	memcpy(u3, f->su2link[ l->next[i][dir2] ][dir1], SU2LINK * sizeof(*u3));
	memcpy(u4, f->su2link[i][dir2], SU2LINK * sizeof(*u4));

	// conjugate u3 and u4
	for (int k=1; k<SU2LINK; k++) {
		u3[k] = -1.0 * u3[k];
		u4[k] = -1.0 * u4[k];
	}

	su2rot(u3, u4); // u3 <- u3.u4
	su2rot(u2, u3); // u2 <- u2.u3
	su2rot(u1, u2); // u1 <- u1.u2
}


/* Calculate untraced staple in counterclockwise direction ("usual" direction),
*  and store it in V.
* Specifically, calculates:
* 	U1 U2^+ U3^+
*/
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
* where mu = dir.
*/
void su2staple_wilson(lattice const* l, fields const* f, params const* p, long i, int dir, double* V) {
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
* Triplet is not included here, because its' hopping term is quadratic in the link
* and cannot be expressed as a simple staple.
*/
void su2link_staple(lattice const* l, fields const* f, params const* p, long i, int dir, double* V) {

	su2staple_wilson(l, f, p, i, dir, V);
	for (int k=0; k<4; k++) {
		V[k] *= -0.5 * p->betasu2;
	}
	#ifdef HIGGS
		// add hopping term: -Tr U_j Phi(x+j) exp(-i a_j(x) sigma_3) Phi(x)^+,
		// a_j(x) = 0 if hypercharge is neglected.
		double nextphi[4];
		double currentphi[4];
		long nextsite = l->next[i][dir];
		// we want Hermitian conjugate of Phi(x):

		currentphi[0] = f->su2doublet[i][0];
		currentphi[1] = -1.0*f->su2doublet[i][1];
		currentphi[2] = -1.0*f->su2doublet[i][2];
		currentphi[3] = -1.0*f->su2doublet[i][3];
		nextphi[0] = f->su2doublet[nextsite][0];
		nextphi[1] = f->su2doublet[nextsite][1];
		nextphi[2] = f->su2doublet[nextsite][2];
		nextphi[3] = f->su2doublet[nextsite][3];

		#ifdef U1
			// the U(1) contribution can be written as
			// I cos(a) - i sin(a) sigma_3, so in our notation it can be
			// treated as a doublet field with components
			// a[0] = sqrt(2) cos(a), a[1] = 0, a[2] = 0, a[3] = -sqrt(2) sin(a)
			double s = sin(f->u1link[i][dir]);
			double c = cos(f->u1link[i][dir]);
			double b[4];
			for (int k=0; k<4; k++) {
				b[k] = nextphi[k];
			}
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
		 -0.5 times what su2rot() of two doublets gives us.
		*/
		su2rot(nextphi, currentphi);
		for (int k=0; k<4; k++) {
			V[k] -= 0.5 * nextphi[k];
		}
	#endif
}


/* Calculate Wilson action for a single SU(2) link.  Explicitly, calculates:
* beta * \Sum_{i < j} ( 1 - 0.5 Re Tr U_i (x) U_j (x+i) U_i(x+j)^+ U_j(x)^+ )
*/
long double local_su2wilson(lattice const* l, fields const* f, params const* p, long i) {

	long double res = 0.0;

	for (int dir1 = 0; dir1 < l->dim; dir1++) {
		for (int dir2 = 0; dir2 < dir1; dir2++ ) {
			res += (1.0 - 0.5 * su2ptrace(l, f, i, dir2, dir1) );
		}
	}

	res = p->betasu2 * res;

	return res;

}


/* Calculate the contribution to the total action due to a single SU(2) link.
* This requires calculating the plaquette trace in two plaquettes in each plane
* that contain the link at site i.
*
* The constant term in beta \sum (1 - 0.5 ptrace) is also included for convenience.
*/
double localact_su2link(lattice const* l, fields const* f, params const* p, long i, int dir) {

	double tot = 0.0;

	for (int dir2 = 0; dir2<l->dim; dir2++) {
		if (dir2 != dir) {
			tot += (1.0 - 0.5 * su2ptrace(l, f, i, dir, dir2));
			tot += (1.0 - 0.5 * su2ptrace(l, f, l->prev[i][dir2], dir, dir2));
		}
	}
	tot *= p->betasu2;

/* same using staples, for debugging. Both work.
	double staple[4];
	double V[4];
	memcpy(V, f.su2link[i][dir], SU2LINK*sizeof(double));
	su2staple_wilson(f, p, i, dir, staple);
	su2rot(V,staple);

	tot = p.betasu2 * (1.0 * p.dim * 0.5 - 0.5*2*V[0]);
	*/

	// hopping terms:
	#ifdef HIGGS
		tot += hopping_doublet_forward(l, f, p, i, dir);
	#endif
	#ifdef TRIPLET
		tot += hopping_triplet_forward(l, f, p, i, dir);
	#endif

	return tot;
}

/**********************************
* 	Routines for U(1) fields		*
***********************************/

/* Calculate plaquette trace in the (dir1, dir2) plane
*	of U(1) link at lattice site i, i.e. same as su2ptrace()
* but for U(1). Note that since we take the real part,
* the result is just a cosine.
*/
double u1ptrace(lattice const* l, fields const* f, long i, int dir1, int dir2) {

	double u1 = f->u1link[i][dir1];
	double u2 = f->u1link[ l->next[i][dir1] ][dir2];
	double u3 = f->u1link[ l->next[i][dir2] ][dir1];
	double u4 = f->u1link[i][dir2];

	return cos(u1 + u2 - u3 - u4);

}

/* Calculate Wilson action for a single U(1) link. Construction is so that
* a loop over sites gives the total Wilson action. This is NOT the full Contribution
* due to a single link, so do NOT use this in update algorithms.
*  Specifically, calculates:
* beta_U1 * \Sum_{i < j} [1 - cos(a_i(x) + a_j(x+i) - a_i(x+j) - a_j(x)]
*/
double local_u1wilson(lattice const* l, fields const* f, params const* p, long i) {

	double res = 0.0;

	for (int dir1 = 0; dir1 < l->dim; dir1++) {
		for (int dir2 = 0; dir2 < dir1; dir2++ ) {
			res += (1.0 - u1ptrace(l, f, i, dir2, dir1) );
		}
	}

	return p->betau1 * res;
}


/* Calculate the contribution to the total action due to a single U(1) link.
* This requires calculating the plaquette trace in two plaquettes in each plane
* that contain the link at site i. Used in metropolis updates.
*
* See localact_su2link() for SU(2) version.
*/
double localact_u1link(lattice const* l, fields const* f, params const* p, long i, int dir) {

	double tot = 0.0;

	for (int dir2 = 0; dir2<l->dim; dir2++) {
		if (dir2 != dir) {
			tot += (1.0 - u1ptrace(l, f, i, dir, dir2));
			tot += (1.0 - u1ptrace(l, f, l->prev[i][dir2], dir, dir2));
		}
	}
	tot *= p->betau1;

	// hopping terms:
	#ifdef HIGGS
		tot += hopping_doublet_forward(l, f, p, i, dir);
	#endif

	return tot;
}


/**********************************
* 	Routines for SU(2) doublets		*
***********************************/

/* Calculate 0.5 Tr \Phi^+ \Phi
*	 which corresponds to \phi^\dagger \phi in continuum notation.
*/
double doubletsq(double* a) {
	return 0.5 * (a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);
}


/* Calculate trace of two SU(2) doublets and a SU(2) link. Used in hopping terms.
* Specifically, calculates:
*		Tr \Phi_1^+ U \Phi_2.
* Note that this is always real.
*/
double hopping_trace(double* phi1, double* u, double* phi2) {

	return phi1[0]*phi2[0]*u[0] + phi1[1]*phi2[1]*u[0] + phi1[2]*phi2[2]*u[0] +
   phi1[3]*phi2[3]*u[0] + phi1[1]*phi2[0]*u[1] - phi1[0]*phi2[1]*u[1] -
   phi1[3]*phi2[2]*u[1] + phi1[2]*phi2[3]*u[1] + phi1[2]*phi2[0]*u[2] +
   phi1[3]*phi2[1]*u[2] - phi1[0]*phi2[2]*u[2] - phi1[1]*phi2[3]*u[2] +
   phi1[3]*phi2[0]*u[3] - phi1[2]*phi2[1]*u[3] + phi1[1]*phi2[2]*u[3] -
   phi1[0]*phi2[3]*u[3];
}

/* Same as hopping_trace(), but includes hypercharge.
* Specifically, calculates:
*		Re Tr \Phi_1^+ U \Phi_2 exp[-i a sigma_3].
*/
double hopping_trace_su2u1(double* phi1, double* u, double* phi2, double a) {
	double s = sin(a);
	double c = cos(a);

	return c*phi1[0]*phi2[0]*u[0] - s*phi1[3]*phi2[0]*u[0] + c*phi1[1]*phi2[1]*u[0] -
   s*phi1[2]*phi2[1]*u[0] + s*phi1[1]*phi2[2]*u[0] + c*phi1[2]*phi2[2]*u[0] +
   s*phi1[0]*phi2[3]*u[0] + c*phi1[3]*phi2[3]*u[0] + c*phi1[1]*phi2[0]*u[1] -
   s*phi1[2]*phi2[0]*u[1] - c*phi1[0]*phi2[1]*u[1] + s*phi1[3]*phi2[1]*u[1] -
   s*phi1[0]*phi2[2]*u[1] - c*phi1[3]*phi2[2]*u[1] + s*phi1[1]*phi2[3]*u[1] +
   c*phi1[2]*phi2[3]*u[1] + s*phi1[1]*phi2[0]*u[2] + c*phi1[2]*phi2[0]*u[2] +
   s*phi1[0]*phi2[1]*u[2] + c*phi1[3]*phi2[1]*u[2] - c*phi1[0]*phi2[2]*u[2] +
   s*phi1[3]*phi2[2]*u[2] - c*phi1[1]*phi2[3]*u[2] + s*phi1[2]*phi2[3]*u[2] +
   s*phi1[0]*phi2[0]*u[3] + c*phi1[3]*phi2[0]*u[3] - s*phi1[1]*phi2[1]*u[3] -
   c*phi1[2]*phi2[1]*u[3] + c*phi1[1]*phi2[2]*u[3] - s*phi1[2]*phi2[2]*u[3] -
   c*phi1[0]*phi2[3]*u[3] + s*phi1[3]*phi2[3]*u[3];

}


/* Calculate the hopping term for a SU(2) doublet at site i,
* in the "forward" direction. Specifically, calculates
*		-Tr \Phi(x)^+ U_j(x) \Phi(x+j) exp(-i \alpha_j(x) \sigma_3)
* for j = dir and \alpha_j(x) = 0 if hypercharge is neglected.
*/
double hopping_doublet_forward(lattice const* l, fields const* f, params const* p, long i, int dir) {
	double *phi1 = f->su2doublet[i];
	double *phi2 = NULL;
	double *U = NULL;
	double tot = 0.0;

	phi2 = f->su2doublet[(l->next[i][dir])];
	U = f->su2link[i][dir];

	#ifndef U1
	// no U(1)
	tot -= hopping_trace(phi1, U, phi2);
	#else
	// include U(1)
	tot -= hopping_trace_su2u1(phi1, U, phi2, f->u1link[i][dir]);
	#endif

	return tot;
}

/* Calculate the hopping term for a SU(2) doublet at site i,
* in the "backwards" direction. Specifically, calculates
*		- Tr \Phi(x-j)^+ U_j(x-j) \Phi(x) exp(-i \alpha_j(x-j) \sigma_3)
* for j = dir and \alpha_j(x-j) = 0 if hypercharge is neglected.
*/
double hopping_doublet_backward(lattice const* l, fields const* f, params const* p, long i, int dir) {
	double *phi1 = NULL;
	double *phi2 = f->su2doublet[i];
	double *U = NULL;
	double tot = 0.0;

	long previous = l->prev[i][dir];

	phi1 = f->su2doublet[previous];
	U = f->su2link[previous][dir];

	#ifndef U1
	// no U(1)
	tot -= hopping_trace(phi1, U, phi2);
	#else
	// include U(1)
	tot -= hopping_trace_su2u1(phi1, U, phi2, f->u1link[previous][dir]);
	#endif

	return tot;
}

/* Calculate the full covariant derivative for a SU(2) doublet at site i,
* in the "forward" directions. Specifically, calculates
*		\sum_j [ Tr\Phi(x)^+ \Phi(x) - Tr \Phi(x)^+ U_j(x) \Phi(x+j) exp(-i \alpha_j(x) \sigma_3) ]
*/
double covariant_doublet(lattice const* l, fields const* f, params const* p, long i) {
	double tot = 0.0;
	double mod = doubletsq(f->su2doublet[i]);
	for (int dir=0; dir<l->dim; dir++){
		// multiply by 2 here because doubletsq gives 0.5 Tr Phi^+ Phi
		tot += 2.0 * mod + hopping_doublet_forward(l, f, p, i, dir);
	}

	return tot;
}


/* Calculate the Higgs potential for a su2doublet field at site i.
* Specifically: V = 0.5 m^2 Tr(\Phi^+ \Phi) + 0.25 * \lambda (Tr(\Phi^+ \Phi))^2 in the SM.
* Contributions from other fields that couple to Higgs are also included.
* Used in localact_doublet(), which is used in metropolis update.
*
*/
double higgspotential(fields const* f, params const* p, long i) {

	double mod = doubletsq(f->su2doublet[i]);
	double pot = p->msq_phi * mod + p->lambda_phi * mod*mod;

	#ifdef TRIPLET
	// add term 0.5 Tr Phi^+ Phi Tr A^2
		pot += p->a2 * mod * tripletsq(f->su2triplet[i]);
	#endif

	return pot;
}

/* Calculate the action due to single su2doublet field at site i.
* This includes the potential, as well as hopping terms
* in "forward" and "backwards" directions.
* Used in metropolis update.
*/
double localact_doublet(lattice const* l, fields const* f, params const* p, long i) {

	double tot = 0.0;
	// Full covariant derivative with the local Tr Phi^+ Phi included,
	// and summed over directions:
	tot += covariant_doublet(l, f, p, i);

	for (int dir=0; dir<l->dim; dir++) {
		// contribution from backwards hopping terms:
		tot += hopping_doublet_backward(l, f, p, i, dir);
	}

	tot += higgspotential(f, p, i);

	return tot;
}

/**********************************
* 	Routines for SU(2) triplets		*
***********************************/
// We assume zero hypercharge for these

/* Calculate Tr A^2 for an adjoint field.
*	 This corresponds to 0.5 A^a A^a in continuum notation. */
double tripletsq(double* a) {
	return 0.5*(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

/* Calculate trace of two SU(2) doublets and two SU(2) links. Used in adjoint hopping terms.
* Specifically, calculates:
*		Tr A1 U A2 U^+.
* Note that this is always real.
*/
double hopping_trace_triplet(double* a1, double* u, double* a2) {

	return 0.5 * ( a1[0]*a2[0]*(u[0]*u[0]) + a1[1]*a2[1]*(u[0]*u[0]) + a1[2]*a2[2]*(u[0]*u[0]) -
   2.0*a1[2]*a2[1]*u[0]*u[1] + 2.0*a1[1]*a2[2]*u[0]*u[1] + a1[0]*a2[0]*(u[1]*u[1]) -
   a1[1]*a2[1]*(u[1]*u[1]) - a1[2]*a2[2]*(u[1]*u[1]) + 2.0*a1[2]*a2[0]*u[0]*u[2] -
   2.0*a1[0]*a2[2]*u[0]*u[2] + 2.0*a1[1]*a2[0]*u[1]*u[2] + 2*a1[0]*a2[1]*u[1]*u[2] -
   a1[0]*a2[0]*(u[2]*u[2]) + a1[1]*a2[1]*(u[2]*u[2]) - a1[2]*a2[2]*(u[2]*u[2]) -
   2.0*a1[1]*a2[0]*u[0]*u[3] + 2.0*a1[0]*a2[1]*u[0]*u[3] + 2.0*a1[2]*a2[0]*u[1]*u[3] +
   2.0*a1[0]*a2[2]*u[1]*u[3] + 2.0*a1[2]*a2[1]*u[2]*u[3] + 2.0*a1[1]*a2[2]*u[2]*u[3] -
   a1[0]*a2[0]*(u[3]*u[3]) - a1[1]*a2[1]*(u[3]*u[3]) + a1[2]*a2[2]*(u[3]*u[3]) );
}

/* Calculate the hopping term for a SU(2) triplet at site i,
* in the "forward" direction. Specifically, calculates
*		-2 Tr A(x) U_j(x) A(x+j) U_j(x)^+
* for j = dir.
*/
double hopping_triplet_forward(lattice const* l, fields const* f, params const* p, long i, int dir) {
	double *a1 = f->su2triplet[i];
	double *a2 = NULL;
	double *U = NULL;
	double tot = 0.0;

	a2 = f->su2triplet[l->next[i][dir]];
	U = f->su2link[i][dir];

	tot -= 2.0 * hopping_trace_triplet(a1, U, a2);

	return tot;
}

/* Calculate the hopping term for a SU(2) triplet at site i,
* in the "backwards" direction. Specifically, calculates
*		-2 Tr A(x-j) U_j(x-j) A(x) U_j(x-j)^+
* for j = dir.
*/
double hopping_triplet_backward(lattice const* l, fields const* f, params const* p, long i, int dir) {
	double *a1 = NULL;
	double *a2 = f->su2triplet[i];
	double *U = NULL;
	double tot = 0.0;

	long previous = l->prev[i][dir];

	a1 = f->su2triplet[previous];
	U = f->su2link[previous][dir];

	tot -= 2.0 * hopping_trace_triplet(a1, U, a2);

	return tot;
}

/* Calculate the full covariant derivative for a SU(2) triplet at site i,
* in the "forward" directions. Specifically, calculates
*		2 \sum_j [ Tr A^2 - Tr A(x) U_j(x) A(x+j) U_j(x)^+ ]
* using hopping_triplet_forward().
*/
double covariant_triplet(lattice const* l, fields const* f, params const* p, long i) {
	double tot = 0.0;
	double mod = tripletsq(f->su2triplet[i]);
	for (int dir=0; dir<l->dim; dir++) {
		tot += 2.0 * mod + hopping_triplet_forward(l, f, p, i, dir);
	}

	return tot;
}

/* Calculate the action due to single su2triplet field at site i.
* This includes the potential, as well as hopping terms
* in "forward" and "backwards" directions.
* Used in metropolis update.
*/
double localact_triplet(lattice const* l, fields const* f, params const* p, long i) {

	double tot = 0.0;
	// Full covariant derivative with the local 2*Tr A^2 included,
	// and summed over directions:
	tot += covariant_triplet(l, f, p, i);

	for (int dir=0; dir<l->dim; dir++) {
		// contribution from backwards hopping terms:
		tot += hopping_triplet_backward(l, f, p, i, dir);
	}

	// add potential
	double mod = tripletsq(f->su2triplet[i]);

	tot += p->msq_triplet * mod + p->b4 * mod * mod;
	#ifdef HIGGS
		// add term 0.5 Tr Phi^+ Phi Tr A^2
		tot += p->a2 * doubletsq(f->su2doublet[i]) * mod;
	#endif

	return tot;
}


/* Smearing routines. These construct a smeared field at each site
* by averaging over the field and its covariant connection with nearest neighbors.
* Naturally used with blocking routines.
* In all routines here the input
* smear_dir[j] = 1 if the direction j is to be smeared, 0 otherwise. */


/* Smear the SU(2) link at site i and store in res (link components).
* This is done by calculating extended staples in the smearing directions,
* including only pure gauge contributions.
* See M. Teper, Phys.Lett.B 183 (1987);
* however Kari seems to use slightly different blocking where he
* separately calculates two staples:
* 	U_i(y) = V_i(x) V_i(x+i), where
*		V_i(x) = U_i(x) + \sum_j U_j(x) U_i(x+j) U_j^+(x+i)
* with the sum running over blocked directions, including backwards dirs, and j != i.
* Kari has normalization factor of 1/3 in V_i because for i=1,2, the original link
* plus a backwards and forward staple contributes. I am using Kari's method here.
*/
void smear_link(lattice const* l, fields const* f, int const* smear_dir, double* res, long i, int dir) {

	if (!smear_dir[dir]) {
		printf("WARNING: smearing gauge link without smearing the lattice dimension (in su2u1.c)\n");
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
		res[k] /= paths;
		v2[k] /= paths;
	}

	su2rot(res, v2); // res <- V1.V2 = V_dir(x) V_dir(x+dir)
	// in general this is unitary but not necessarily in SU(2). so normalize again:
	double det = su2sqr(res);
	det = sqrt(det);
	for (int k=0; k<SU2LINK; k++) {
		res[k] /= det;
	}

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

	// "upper" staple (U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
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
	u3 = f->su2link[ l->prev[i][nu] ][mu];
	su2staple_clockwise(tot, u1, u2, u3);;
	for(int k=0; k<SU2LINK; k++) {
		// take conjugate if needed
		if (dagger && k != 0) tot[k] = -1.0 * tot[k];
		res[k] += tot[k];
	}

}


/* Smear the triplet field at site i and store in res (triplet components).
* This is done by calculating
* 	Sigma(x) + sum_j U_j(x) Sigma(x+j) U^+_j (x)
* and normalizing with the number of sites. The sum is over
* directions specified in smear_dir, and contains both forward and backward terms.
* Here U_{-j}(x) = U^+_j(x-j). */
void smear_triplet(lattice const* l, fields const* f, int const* smear_dir, double* res, long i) {

	int sites = 1; // how many sites are involved in smearing
	double cov[SU2TRIP] = {0.0};

	for (int dir=0; dir<l->dim; dir++) {
		if (smear_dir[dir]) {
			// forward connection
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

			// backward connection
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
		res[k] /= sites;
	}
	// done
}


/* Smear all fields and store in f_b. Does not smear fields at odd sites in smearing
* directions, because those are not needed for blocked lattices.
*/
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
      #ifdef HIGGS

      #endif
      #ifdef TRIPLET
        smear_triplet(l, f, block_dir, f_b->su2triplet[i], i);
      #endif

    }

  } // end i
}

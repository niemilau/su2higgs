/** @file su2.c
*
* Routines for operating on SU(2) fields and calculating the SU(2) action.
* Overall conventions are taken from hep-lat/9504001. Specifically:
*
*	Everywhere in this code, the SU(2) links are parameterized with real numbers as
* 	U = u0 I + i ( u1 sig1 + u2 sig2 + u3 sig3 ),
* where: i is the imaginary unit, I is the unit matrix, sig are the Pauli matrices.
*
*	In general such U is a U(2) matrix, and the condition for SU(2) is
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
*
*	My potential is  V = m^2 Tr A^2 + b4 (Tr A^2)^2 + 0.5 Tr(\Phi^+ \Phi) Tr A^2,
*	where the normalization matches what is usually used in continuum for the A[a].
* This matches ref. hep-ph/9704416 but is different from that of hep-lat/9504001.
*
//TODO
*
*/

#include "su2.h"


// calculate the determinant, or norm, of a matrix in our SU(2) parametrization
// mainly for debugging purposes
double su2sqr(double *u) {
	return u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
}


// rotate U1 from the right by U2 and store again in U1: U1 = U1.U2.
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
inline double su2trace4(double *u1, double *u2, double *u3, double *u4) {

	return 2 * (u1[0]*u2[0]*u3[0]*u4[0] - u1[1]*u2[1]*u3[0]*u4[0] - u1[2]*u2[2]*u3[0]*u4[0] -
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
*	where mu = dir1, nu = dir2, x = site at index i
*/
double su2ptrace(fields f, params p, long i, int dir1, int dir2) {

	// set pointers to point to the component arrays of the links
	double *u1 = &(*f.su2link[i][dir1]);
	double *u2 = &(*f.su2link[ p.next[i][dir1] ][dir2]);
	double *u3 = &(*f.su2link[ p.next[i][dir2] ][dir1]);
	double *u4 = &(*f.su2link[i][dir2]);

	return su2trace4(u1, u2, u3, u4);

}


/* Calculate untraced staple in counterclockwise direction ("usual" direction),
*  and store it in V.
* Specifically, calculates:
* 	U1 U2^+ U3^+
*/
inline void su2staple_counterwise(double* V, double* u1, double* u2, double* u3) {

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
inline void su2staple_clockwise(double* V, double* u1, double* u2, double* u3) {

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
void su2staple_wilson(fields f, params p, long i, int dir, double* V) {
	double tot[4] = {0.0, 0.0, 0.0, 0.0};
	double* u1 = NULL;
	double* u2 = NULL;
	double* u3 = NULL;

	for (int j=0; j<p.dim; j++) {
		if (j != dir) {
			// "upper" staple
			u1 = f.su2link[ p.next[i][dir] ][j];
			u2 = f.su2link[ p.next[i][j] ][dir];
			u3 = f.su2link[i][j];
			su2staple_counterwise(V, u1, u2, u3);
			for(int k=0; k<4; k++){
				tot[k] += V[k];
			}
			// "lower" staple
			u1 = f.su2link[ p.prev[(p.next[i][dir])][j] ][j];
			u2 = f.su2link[(p.prev[i][j])][dir];
			u3 = f.su2link[(p.prev[i][j])][j];
			su2staple_clockwise(V, u1, u2, u3);
			//su2staple_counterwise(V, u1, u2, u3);
			for(int k=0; k<4; k++){
				tot[k] += V[k];
			}
		}
	}

	for(int k=0; k<4; k++){
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
*/
void su2link_staple(fields f, params p, long i, int dir, double* V) {

	su2staple_wilson(f, p, i, dir, V);
	for (int k=0; k<4; k++) {
		V[k] *= -0.5 * p.betasu2;
	}
	#ifdef HIGGS
		// add hopping term: -Tr U_j \Phi(x+j)\Phi(x)^+
		double nextphi[4];
		double currentphi[4];
		long nextsite = p.next[i][dir];
		// we want Hermitian conjugate of Phi(x):

		currentphi[0] = f.su2doublet[i][0];
		currentphi[1] = -1.0*f.su2doublet[i][1];
		currentphi[2] = -1.0*f.su2doublet[i][2];
		currentphi[3] = -1.0*f.su2doublet[i][3];
		nextphi[0] = f.su2doublet[nextsite][0];
		nextphi[1] = f.su2doublet[nextsite][1];
		nextphi[2] = f.su2doublet[nextsite][2];
		nextphi[3] = f.su2doublet[nextsite][3];

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
long double local_su2wilson(fields f, params p, long i) {

	long double res = 0.0;

	for (int dir1 = 0; dir1 < p.dim; dir1++) {
		for (int dir2 = 0; dir2 < dir1; dir2++ ) {
			res += (1.0 - 0.5 * su2ptrace(f, p, i, dir2, dir1) );
		}
	}

	res = p.betasu2 * res;

	return res;

}


/* Calculate the contribution to the total action due to a single SU(2) link.
* This requires calculating the plaquette trace in two plaquettes in each plane
* that contain the link at site i.
*
* The constant term in beta \sum (1 - 0.5 ptrace) is also included for convenience.
// TODO remove the extra term for optimization?
*/
double localact_su2link(fields f, params p, long i, int dir) {

	double tot = 0.0;

	for (int dir2 = 0; dir2<p.dim; dir2++) {
		if (dir2 != dir) {
			tot += (1.0 - 0.5 * su2ptrace(f, p, i, dir, dir2));
			tot += (1.0 - 0.5 * su2ptrace(f, p, p.prev[i][dir2], dir, dir2));
		}
	}
	tot *= p.betasu2;

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
		tot += hopping_doublet_forward(f, p, i, dir);
	#endif
	#ifdef TRIPLET
		tot += hopping_triplet_forward(f, p, i, dir);
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
*		Tr \Phi_1^+ U_ \Phi_2.
* Note that this is always real.
*/
inline double hopping_trace(double* phi1, double* u, double* phi2) {

	return phi1[0]*phi2[0]*u[0] + phi1[1]*phi2[1]*u[0] + phi1[2]*phi2[2]*u[0] +
   phi1[3]*phi2[3]*u[0] + phi1[1]*phi2[0]*u[1] - phi1[0]*phi2[1]*u[1] -
   phi1[3]*phi2[2]*u[1] + phi1[2]*phi2[3]*u[1] + phi1[2]*phi2[0]*u[2] +
   phi1[3]*phi2[1]*u[2] - phi1[0]*phi2[2]*u[2] - phi1[1]*phi2[3]*u[2] +
   phi1[3]*phi2[0]*u[3] - phi1[2]*phi2[1]*u[3] + phi1[1]*phi2[2]*u[3] -
   phi1[0]*phi2[3]*u[3];
}


/* Calculate the hopping term for a SU(2) doublet at site i,
* in the "forward" direction. Specifically, calculates
*		-Tr \Phi(x)^+ U_j(x) \Phi(x+j)
* for j = dir.
*/
double hopping_doublet_forward(fields f, params p, long i, int dir) {
	double *phi1 = f.su2doublet[i];
	double *phi2 = NULL;
	double *U = NULL;
	double tot = 0.0;

	phi2 = f.su2doublet[(p.next[i][dir])];
	U = f.su2link[i][dir];

	tot -= hopping_trace(phi1, U, phi2);

	return tot;
}

/* Calculate the hopping term for a SU(2) doublet at site i,
* in the "backwards" direction. Specifically, calculates
*		- Tr \Phi(x-j)^+ U_j(x-j) \Phi(x)
* for j = dir.
*/
double hopping_doublet_backward(fields f, params p, long i, int dir) {
	double *phi1 = NULL;
	double *phi2 = f.su2doublet[i];
	double *U = NULL;
	double tot = 0.0;

	long previous = p.prev[i][dir];

	phi1 = f.su2doublet[previous];
	U = f.su2link[previous][dir];

	tot -= hopping_trace(phi1, U, phi2);

	return tot;
}

/* Calculate the covariant derivative for a SU(2) doublet at site i,
* in the "forward" directions. Specifically, calculates
*		\sum_j [ Tr\Phi(x)^+ \Phi(x) - Tr \Phi(x)^+ U_j(x) \Phi(x+j) ]
*/
double covariant_doublet(fields f, params p, long i) {
	double tot = 0.0;
	double mod = doubletsq(f.su2doublet[i]);
	for (int dir=0; dir<p.dim; dir++){
		// multiply by 2 here because doubletsq gives 0.5 Tr Phi^+ Phi
		tot += 2.0 * mod + hopping_doublet_forward(f, p, i, dir);
	}

	return tot;
}


/* Calculate the Higgs potential for a su2doublet field at site i.
* Specifically: V = 0.5 m^2 Tr(\Phi^+ \Phi) + 0.25 * \lambda (Tr(\Phi^+ \Phi))^2 in the SM.
* Contributions from other fields that couple to Higgs are also included.
* Used in localact_doublet(), which is used in metropolis update.
*
*/
double higgspotential(fields f, params p, long i) {

	double mod = doubletsq(f.su2doublet[i]);
	double pot = p.msq_phi * mod + p.lambda_phi * mod*mod;

	#ifdef TRIPLET
	// add term 0.5 Tr \Phi^+ \Phi Tr A^2
		pot += p.a2 * mod * tripletsq(f.su2triplet[i]);
	#endif

	return pot;
}

/* Calculate the action due to single su2doublet field at site i.
* This includes the potential, as well as hopping terms
* in "forward" and "backwards" directions.
* Used in metropolis update.
*/
double localact_doublet(fields f, params p, long i) {

	double tot = 0.0;
	// Full covariant derivative with the local Tr Phi^+ Phi included,
	// and summed over directions:
	tot += covariant_doublet(f, p, i);

	for (int dir=0; dir<p.dim; dir++) {
		// contribution from backwards hopping terms:
		tot += hopping_doublet_backward(f, p, i, dir);
	}

	tot += higgspotential(f, p, i);

	return tot;
}

/**********************************
* 	Routines for SU(2) triplets		*
***********************************/

/* Calculate Tr A^2 for an adjoint field.
*	 This corresponds to 0.5 A^a A^a in continuum notation.
*/
double tripletsq(double* a) {
	return 0.5*(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

/* Calculate trace of two SU(2) doublets and two SU(2) links. Used in adjoint hopping terms.
* Specifically, calculates:
*		Tr A1 U A2 U^+.
* Note that this is always real.
*/
inline double hopping_trace_triplet(double* a1, double* u, double* a2) {

	return 0.5 * ( a1[0]*a2[0]*(u[0]*u[0]) + a1[1]*a2[1]*(u[0]*u[0]) + a1[2]*a2[2]*(u[0]*u[0]) -
   2*a1[2]*a2[1]*u[0]*u[1] + 2*a1[1]*a2[2]*u[0]*u[1] + a1[0]*a2[0]*(u[1]*u[1]) -
   a1[1]*a2[1]*(u[1]*u[1]) - a1[2]*a2[2]*(u[1]*u[1]) + 2*a1[2]*a2[0]*u[0]*u[2] -
   2*a1[0]*a2[2]*u[0]*u[2] + 2*a1[1]*a2[0]*u[1]*u[2] + 2*a1[0]*a2[1]*u[1]*u[2] -
   a1[0]*a2[0]*(u[2]*u[2]) + a1[1]*a2[1]*(u[2]*u[2]) - a1[2]*a2[2]*(u[2]*u[2]) -
   2*a1[1]*a2[0]*u[0]*u[3] + 2*a1[0]*a2[1]*u[0]*u[3] + 2*a1[2]*a2[0]*u[1]*u[3] +
   2*a1[0]*a2[2]*u[1]*u[3] + 2*a1[2]*a2[1]*u[2]*u[3] + 2*a1[1]*a2[2]*u[2]*u[3] -
   a1[0]*a2[0]*(u[3]*u[3]) - a1[1]*a2[1]*(u[3]*u[3]) + a1[2]*a2[2]*(u[3]*u[3]) );
}

/* Calculate the hopping term for a SU(2) triplet at site i,
* in the "forward" direction. Specifically, calculates
*		-2 Tr A(x) U_j(x) A(x+j) U_j(x)^+
* for j = dir.
*/
double hopping_triplet_forward(fields f, params p, long i, int dir) {
	double *a1 = f.su2triplet[i];
	double *a2 = NULL;
	double *U = NULL;
	double tot = 0.0;

	a2 = f.su2triplet[p.next[i][dir]];
	U = f.su2link[i][dir];

	tot -= 2.0 * hopping_trace_triplet(a1, U, a2);

	return tot;
}

/* Calculate the hopping term for a SU(2) triplet at site i,
* in the "backwards" direction. Specifically, calculates
*		-2 Tr A(x-j) U_j(x-j) A(x) U_j(x-j)^+
* for j = dir.
*/
double hopping_triplet_backward(fields f, params p, long i, int dir) {
	double *a1 = NULL;
	double *a2 = f.su2triplet[i];
	double *U = NULL;
	double tot = 0.0;

	long previous = p.prev[i][dir];

	a1 = f.su2triplet[previous];
	U = f.su2link[previous][dir];

	tot -= 2.0 * hopping_trace_triplet(a1, U, a2);

	return tot;
}

/* Calculate the covariant derivative for a SU(2) triplet at site i,
* in the "forward" directions. Specifically, calculates
*		2 \sum_j [ Tr A^2 - Tr A(x) U_j(x) A(x+j) U_j(x)^+ ]
* using hopping_triplet_forward().
*/
double covariant_triplet(fields f, params p, long i) {
	double tot = 0.0;
	double mod = tripletsq(f.su2triplet[i]);
	for (int dir=0; dir<p.dim; dir++) {
		tot += 2.0 * mod + hopping_triplet_forward(f, p, i, dir);
	}
	
	return tot;
}

/* Calculate the action due to single su2triplet field at site i.
* This includes the potential, as well as hopping terms
* in "forward" and "backwards" directions.
* Used in metropolis update.
*/
double localact_triplet(fields f, params p, long i) {

	double tot = 0.0;
	// Full covariant derivative with the local 2*Tr A^2 included,
	// and summed over directions:
	tot += covariant_triplet(f, p, i);

	for (int dir=0; dir<p.dim; dir++) {
		// contribution from backwards hopping terms:
		tot += hopping_triplet_backward(f, p, i, dir);
	}

	// add potential
	double mod = tripletsq(f.su2triplet[i]);
	
	tot += p.msq_triplet * mod + p.b4 * mod * mod;
	#ifdef HIGGS
		// add term 0.5 Tr \Phi^+ \Phi Tr A^2
		tot += p.a2 * doubletsq(f.su2doublet[i]) * mod;
	#endif

	return tot;
}
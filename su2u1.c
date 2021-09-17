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
*	Usually in continuum we have phi = (H1, H2)^T ; the matrix parametrization is then
*
* Phi = (H2* , H1)
*				(-H1*, H2)
*
*	Real SU(2) triplets are parametrized as:
*		A = \frac12 A[a] sigma[a],  a = 1,2,3.
* This matches ref. hep-ph/9704416 but is different from that of hep-lat/9504001,
* which uses A = A[a] sigma[a] and has the action written in a slightly different,
* but equal, form.
*
*	My potential is  V = m^2 Tr A^2 + b4 (Tr A^2)^2 + 0.5 a2 Tr(\Phi^+ \Phi) Tr A^2,
*	where the normalization matches what is usually used in continuum for the A[a].
*
*	U(1) links are written simply as U_j(x) = exp(i a_j(x)), where a_j(x) is real
* and stored in u1link[x][j]. I use a compact formulation (gauge group really is U(1)).
* Higgs hypercharge is normalized to Y=1 in the code.
*
* For 2 Higgs doublets (flag HIGGS2), I assume a doublet basis where the kinetic
* terms are diagonal and canonically normalized.
*/

#include "su2.h"


// calculate the determinant, or norm, of a matrix in our SU(2) parametrization
double su2sqr(double *u) {
	return u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
}


/* rotate U1 from the right by U2 and store again in U1: U1 <- U1.U2.
* If conjugate != 0, take hermitian conjugate of U2 first */
void su2rot(double *u1, double const* u2, int conjugate) {

	double new1, new2, new3, new4;
	if (conjugate == 0) {
		new1 = u1[0]*u2[0] - u1[1]*u2[1] - u1[2]*u2[2] - u1[3]*u2[3];
		new2 = u1[1]*u2[0] + u1[0]*u2[1] + u1[3]*u2[2] - u1[2]*u2[3];
		new3 = u1[2]*u2[0] - u1[3]*u2[1] + u1[0]*u2[2] + u1[1]*u2[3];
		new4 = u1[3]*u2[0] + u1[2]*u2[1] - u1[1]*u2[2] + u1[0]*u2[3];
	} else {
		// Conjugate corresponds to flipping signs of u2[k] for k != 0
		// mini-optimize by hardcoding the result (as we can't modify u2 here)
		new1 = u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2] + u1[3]*u2[3];
		new2 = u1[1]*u2[0] - u1[0]*u2[1] - u1[3]*u2[2] + u1[2]*u2[3];
		new3 = u1[2]*u2[0] + u1[3]*u2[1] - u1[0]*u2[2] - u1[1]*u2[3];
		new4 = u1[3]*u2[0] - u1[2]*u2[1] + u1[1]*u2[2] - u1[0]*u2[3];
	}
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

	su2rot(u1, u2, 0); // u1 <- u1.u2
	su2rot(u1, u3, 1); // u1 <- u1.u3^+
	su2rot(u1, u4, 1); // u1 <- u1.u4^+
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
	/*
	//same using staples, for debugging. Both work.
	double staple[4];
	double V[4];
	memcpy(V, f.su2link[i][dir], SU2LINK*sizeof(double));
	su2staple_wilson(f, i, dir, staple);
	su2rot(V,staple, 0);

	tot = p.betasu2 * (1.0 * p.dim * 0.5 - 0.5*2*V[0]);
	*/

	// hopping terms:
	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) tot += hopping_doublet_forward(l, f, i, dir, db);
	#endif

	#ifdef TRIPLET
		tot += hopping_triplet_forward(l, f, p, i, dir);
	#endif

	return tot;
}

/* Calculate a simple plaquette "clover" for SU(2) at a given site, see
*	fig 1 in hep-lat/0106023. This calculates their eq. 12 in the (d1, d2) plane.
* Heuristically, this reproduces the field strength tensor at the lattice site
* at O(a), while the standard plaquette gives F_ij in
* the middle of the plaquette. If O_munu is the clover, then for SU(N)
* g F_munu(x) = -i/8 [(O_munu(x) - O^+_munu(x)) - 1/N Tr(O_munu(x) - O^+_munu(x)) ],
* i.e. the antihermitian part of O_munu projected to the Lie algebra. */
void clover_su2(lattice const* l, fields const* f, long i, int d1, int d2, double* clover) {

	su2plaquette(l, f, i, d1, d2, clover); // standard plaquette

	// add other plaquettes in (d1,d2) plane that begin at site i
	double u1[SU2LINK], u2[SU2LINK], u3[SU2LINK], u4[SU2LINK];
	long site;

	// U_nu(x) U^+_mu(x+nu-mu) U^+_nu(x-mu) U_mu(x-mu)
	memcpy(u1, f->su2link[i][d2], SU2LINK * sizeof(*u1));
	site = l->next[i][d2];
	site = l->prev[site][d1]; // x + nu - mu
	memcpy(u2, f->su2link[site][d1], SU2LINK * sizeof(*u1));
	site = l->prev[i][d1]; // x - mu
	memcpy(u3, f->su2link[site][d2], SU2LINK * sizeof(*u1));
	memcpy(u4, f->su2link[site][d1], SU2LINK * sizeof(*u1));
	// take conjugates and multiply
	for (int k=1; k<SU2LINK; k++) {
		u2[k] = -1.0*u2[k];
		u3[k] = -1.0*u3[k];
	}
	su2rot(u3, u4, 0);
	su2rot(u2, u3, 0);
	su2rot(u1, u2, 0);
	for (int k=0; k<SU2LINK; k++) {
		clover[k] += u1[k];
	}

	// U^+_mu(x-mu) U^+_nu(x-nu-mu) U_mu(x-mu-nu) U_nu(x-nu
	site = l->prev[i][d1]; // x - mu
	memcpy(u1, f->su2link[site][d1], SU2LINK * sizeof(*u1));
	site = l->prev[site][d2]; // x - nu - mu
	memcpy(u2, f->su2link[site][d2], SU2LINK * sizeof(*u1));
	memcpy(u3, f->su2link[site][d1], SU2LINK * sizeof(*u1));
	site = l->prev[i][d2]; // x - nu
	memcpy(u4, f->su2link[site][d2], SU2LINK * sizeof(*u1));
	// take conjugates and multiply
	for (int k=1; k<SU2LINK; k++) {
		u1[k] = -1.0*u1[k];
		u2[k] = -1.0*u2[k];
	}
	su2rot(u3, u4, 0);
	su2rot(u2, u3, 0);
	su2rot(u1, u2, 0);
	for (int k=0; k<SU2LINK; k++) {
		clover[k] += u1[k];
	}

	// U^+_nu(x-nu) U_mu(x-nu) U_nu(x+mu-nu)U^+_mu(x)
	site = l->prev[i][d2]; // x - nu
	memcpy(u1, f->su2link[site][d2], SU2LINK * sizeof(*u1));
	memcpy(u2, f->su2link[site][d1], SU2LINK * sizeof(*u1));
	site = l->next[site][d1]; // x + mu - nu
	memcpy(u3, f->su2link[site][d2], SU2LINK * sizeof(*u1));
	memcpy(u4, f->su2link[i][d1], SU2LINK * sizeof(*u1));
	// take conjugates and multiply
	for (int k=1; k<SU2LINK; k++) {
		u1[k] = -1.0*u1[k];
		u4[k] = -1.0*u4[k];
	}
	su2rot(u3, u4, 0);
	su2rot(u2, u3, 0);
	su2rot(u1, u2, 0);
	for (int k=0; k<SU2LINK; k++) {
		clover[k] += u1[k];
	}

	// done
}


#ifdef U1
/**********************************
* 	Routines for U(1) fields		*
***********************************/

/* The link variable is u_mu(x) = e^(i alpha_\mu(x)), and my
* u1link[i][mu] = alpha_\mu at site i. In compact formulation
* the alpha's are angular variables, so restrict to  ]-pi, pi] */

/* My U(1) action is S = betau1 * sum_{x, i<j} (1 - Re p_{ij}^r),
* which is simply a plaquette in r-representation of U(1). r = integer.
* My r = 1/gamma of hep-lat/9705003 and hep-lat/9612006.
* Note that all matter fields need to be in irreps of U(1), so that the charge
* is always an integer. When converting from continuum theory this may require
* a rescaling of the gauge coupling. For example if Higgs covariant derivative is
* D_i \phi = partial_i \phi + i g'/2 B_i \phi and there are no other matter fields,
* then beta_G' = 4 / (a g'^2 r^2) gives the correct continuum action.
*/

/* Calculate plaquette trace in the (dir1, dir2) plane
*	of U(1) link at lattice site i, i.e. same as su2ptrace()
* but for U(1). */
double u1ptrace(lattice const* l, fields const* f, long i, int dir1, int dir2) {

	double u1 = f->u1link[i][dir1];
	double u2 = f->u1link[ l->next[i][dir1] ][dir2];
	double u3 = f->u1link[ l->next[i][dir2] ][dir1];
	double u4 = f->u1link[i][dir2];

	return u1 + u2 - u3 - u4;

}

/* Calculate Wilson action for a single U(1) link. Construction is so that
* a loop over sites gives the total gauge action. This is NOT the full contribution
* due to a single link, so do NOT use this in update algorithms.
*  Specifically, calculates:
* beta_U1 * \Sum_{i < j} {1 - cos[r * (a_i(x) + a_j(x+i) - a_i(x+j) - a_j(x))]} */
double local_u1wilson(lattice const* l, fields const* f, params const* p, long i) {

	double res = 0.0;

	for (int dir1 = 0; dir1 < l->dim; dir1++) {
		for (int dir2 = 0; dir2 < dir1; dir2++ ) {
			double plaq = u1ptrace(l, f, i, dir2, dir1);
			res += (1.0 - cos(p->r_u1 * plaq) );
		}
	}

	return p->betau1 * res;
}


/* Calculate the contribution to the total action due to a single U(1) link.
* This requires calculating the plaquette trace in two plaquettes in each plane
* that contain the link at site i. Used in metropolis updates.
*
* See localact_su2link() for SU(2) version. */
double localact_u1link(lattice const* l, fields const* f, params const* p, long i, int dir) {

	double tot = 0.0;

	for (int dir2 = 0; dir2<l->dim; dir2++) {
		if (dir2 != dir) {
			tot += (1.0 - cos(p->r_u1 * u1ptrace(l, f, i, dir, dir2)));
			tot += (1.0 - cos(p->r_u1 * u1ptrace(l, f, l->prev[i][dir2], dir, dir2)));
		}
	}
	tot *= p->betau1;

	// hopping terms:
	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) tot += hopping_doublet_forward(l, f, i, dir, db);
	#endif

	return tot;
}

#endif // U1

/**********************************
* 	Routines for SU(2) doublets		*
***********************************/

/* Calculate 0.5 Tr \Phi^+ \Phi
*	 which corresponds to \phi^\dagger \phi in continuum notation. */
double doubletsq(double* a) {
	return 0.5 * (a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);
}

/* Calculate product of two doublets f1, f2 in matrix parametrization
* and store again in f1: f1 <- f1.f2. If conj==1, takes conjugate of f1 first */
void phiproduct(double* f1, double const* f2, int conj) {

	double a[SU2DB];
	memcpy(a, f1, SU2DB * sizeof(a[0]));
	if (conj) {
		for (int k=1; k<SU2DB; k++) a[k] *= -1.0;
	}

	f1[0] = (a[0]*f2[0] - a[1]*f2[1] - a[2]*f2[2] - a[3]*f2[3]) / sqrt(2.0);
	f1[1] = (a[1]*f2[0] + a[0]*f2[1] + a[3]*f2[2] - a[2]*f2[3]) / sqrt(2.0);
	f1[2] = (a[2]*f2[0] - a[3]*f2[1] + a[0]*f2[2] + a[1]*f2[3]) / sqrt(2.0);
	f1[3] = (a[3]*f2[0] + a[2]*f2[1] - a[1]*f2[2] + a[0]*f2[3]) / sqrt(2.0);
}

/* Calculate trace of two SU(2) doublets and a SU(2) link. Used in hopping terms.
* Specifically, calculates:
*		Tr \Phi_1^+ U \Phi_2.
* Note that this is always real. */
double hopping_trace(double* phi1, double* u, double* phi2) {

	return phi1[0]*phi2[0]*u[0] + phi1[1]*phi2[1]*u[0] + phi1[2]*phi2[2]*u[0] +
   phi1[3]*phi2[3]*u[0] + phi1[1]*phi2[0]*u[1] - phi1[0]*phi2[1]*u[1] -
   phi1[3]*phi2[2]*u[1] + phi1[2]*phi2[3]*u[1] + phi1[2]*phi2[0]*u[2] +
   phi1[3]*phi2[1]*u[2] - phi1[0]*phi2[2]*u[2] - phi1[1]*phi2[3]*u[2] +
   phi1[3]*phi2[0]*u[3] - phi1[2]*phi2[1]*u[3] + phi1[1]*phi2[2]*u[3] -
   phi1[0]*phi2[3]*u[3];
}

#ifdef U1
/* Same as hopping_trace(), but includes hypercharge.
* Specifically, calculates:
*		Re Tr \Phi_1^+ U \Phi_2 exp[-i Y alpha sigma_3], where alpha is stored in u1link
* Here the Higgs hypercharge is normalized as Y = 1 ALWAYS (if changing this, change also staples.c etc)*/
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
#endif

/* Routines below access f->su2doublet directly, so protect with preprocessor if */
#if (NHIGGS > 0 )

/* Calculate the hopping term for a SU(2) doublet 'phi' at site i,
* in the "forward" direction. Specifically, calculates
*		-Tr \Phi(x)^+ U_j(x) \Phi(x+j) exp(-i Y \alpha_j(x) \sigma_3)
* for j = dir and \alpha_j(x) = 0 if hypercharge is neglected. */
double hopping_doublet_forward(lattice const* l, fields const* f, long i, int dir, int higgs_id) {

	double **higgs = f->su2doublet[higgs_id];

	double *phi1 = higgs[i];
	double *phi2 = higgs[l->next[i][dir]];
	double *U = f->su2link[i][dir];

	double tot = 0.0;

	#ifndef U1
		tot -= hopping_trace(phi1, U, phi2);
	#else
		tot -= hopping_trace_su2u1(phi1, U, phi2, f->u1link[i][dir]);
	#endif

	return tot;
}

/* Calculate the hopping term for a SU(2) doublet 'phi' at site i,
* in the "backwards" direction. Specifically, calculates
*		- Tr \Phi(x-j)^+ U_j(x-j) \Phi(x) exp(-i Y \alpha_j(x-j) \sigma_3)
* for j = dir and \alpha_j(x-j) = 0 if hypercharge is neglected. */
double hopping_doublet_backward(lattice const* l, fields const* f, long i, int dir, int higgs_id) {

	double **higgs = f->su2doublet[higgs_id];

	long prev = l->prev[i][dir];

	double *phi1 = higgs[prev];
	double *phi2 = higgs[i];
	double *U = f->su2link[prev][dir];
	double tot = 0.0;

	#ifndef U1
		// no U(1)
		tot -= hopping_trace(phi1, U, phi2);
	#else
		// include U(1)
		tot -= hopping_trace_su2u1(phi1, U, phi2, f->u1link[prev][dir]);
	#endif

	return tot;

}

/* Calculate the full covariant derivative for an SU(2) doublet at site i,
* in the "forward" directions. Specifically, calculates
*		\sum_j [ Tr\Phi(x)^+ \Phi(x) - Tr \Phi(x)^+ U_j(x) \Phi(x+j) exp(-i Y \alpha_j(x) \sigma_3) ] */
double covariant_doublet(lattice const* l, fields const* f, long i, int higgs_id) {

	double tot = 0.0;
	double **phi = f->su2doublet[higgs_id];
	double mod = doubletsq(phi[i]);
	for (int dir=0; dir<l->dim; dir++){
		// multiply by 2 here because doubletsq gives 0.5 Tr Phi^+ Phi
		tot += 2.0 * mod + hopping_doublet_forward(l, f, i, dir, higgs_id);
	}

	return tot;
}


/* Calculate the action due to single su2doublet field at site i.
* This includes the potential, as well as hopping terms
* in "forward" and "backwards" directions.
* Used in metropolis update. */
double localact_doublet(lattice const* l, fields const* f, params const* p, long i, int higgs_id) {

	double tot = 0.0;
	// Full covariant derivative with the local Tr Phi^+ Phi included,
	// and summed over directions:
	tot += covariant_doublet(l, f, i, higgs_id);

	for (int dir=0; dir<l->dim; dir++) {
		// contribution from backwards hopping terms:
		tot += hopping_doublet_backward(l, f, i, dir, higgs_id);
	}

	tot += higgspotential(f, p, i);

	return tot;
}

/* Calculate phi_{12} = phi1^\dagger phi2, where the fields are now in vector parametrization */
complex get_phi12(double const* h1, double const* h2) {
	complex res;
	/* R = Re phi12, I = Im phi12 = -0.5 i Tr \Phi_1 \sigma_3 \he\Phi_2 */
	res.re = 0.5*(h1[0]*h2[0] + h1[1]*h2[1] + h1[2]*h2[2] + h1[3]*h2[3]);
	res.im = 0.5*(h1[3]*h2[0] + h1[2]*h2[1] - h1[1]*h2[2] - h1[0]*h2[3]);

	return res;
}

#endif // if (NHIGGS > 0)


/* Calculate the full scalar potential at site i, including all scalar fields.
* Specifically: V = 0.5 m^2 Tr(\Phi^+ \Phi) + 0.25 * \lambda (Tr(\Phi^+ \Phi))^2 in the SM.
* Used in localact_doublet(), which is used in metropolis update. */
double higgspotential(fields const* f, params const* p, long i) {

	double pot = 0.0;

	#if (NHIGGS > 0)
		double *h1 = f->su2doublet[0][i];
		double mod = doubletsq(h1);
		pot += p->msq_phi * mod + p->lambda_phi * mod*mod;

		#if (NHIGGS == 2)
			/* Two-Higgs doublet potential is
			* 	V = m1^2 f11 + m2^2 f22 + 0.5(m12^2 f12 + h.c.)
			*			 + lam1 f11^2 + lam2 f22^2 + lam3 f11 f22 + lam4 f12 f21
			*  		 + 0.5(lam5 f12^2 + lam6 f11 f12 + lam7 f22 f21 + h.c. )
			* where f11 = phi1^+.phi1 etc. see documentation for the alternative form used below */

			double *h2 = f->su2doublet[1][i];
			double f11 = mod;
			double f22 = doubletsq(h2);

			complex f12 = get_phi12(h1, h2);
			double R = f12.re; double I = f12.im;

			pot += p->msq_phi2 * f22 + p->m12sq.re * R - p->m12sq.im * I + p->lam2 * f22*f22
					  + p->lam3 * f11*f22 + p->lam4 * (R*R + I*I) + p->lam5.re*(R*R - I*I) - 2.0*p->lam5.im*R*I
						+ f11 * (p->lam6.re*R - p->lam6.im*I) + f22*(p->lam7.re*R + p->lam7.im*I);
		#endif

	#endif

	#ifdef TRIPLET
		// add 0.5 m^2 Tr Sigma^2 + b4 (0.5 Tr Sigma^2)^2, plus couplings to other scalars
		double mod_trip = tripletsq(f->su2triplet[i]); // 0.5 Tr Sigma^2
		pot += p->msq_triplet * mod_trip + p->b4 * mod_trip * mod_trip;
		#if (NHIGGS > 0)
			pot += p->a2 * mod * mod_trip;
		#endif
	#endif

	#ifdef SINGLET
		pot += potential_singlet(f, p, i);
	#endif

	return pot;
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
* Note that this is always real. */
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
* for j = dir. */
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
* for j = dir. */
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
* using hopping_triplet_forward(). */
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
* Used in metropolis update. */
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
	#if (NHIGGS > 0)
		// add term 0.5 Tr Phi^+ Phi Tr A^2
		tot += p->a2 * doubletsq(f->su2doublet[0][i]) * mod;
	#endif

	return tot;
}


#ifdef SINGLET
/**********************************
* 	Routines for singlets		*
***********************************/

/* Singlet has just 1 component at each site, but is stored as a 2D array f.singlet[site][0] */
double get_singlet(double* s) {
	return s[0];
}

double singletsq(double* s) {
	return 0.5 * s[0] * s[0];
}

/* Calculate the action due to singlet field at site i. */
double localact_singlet(lattice const* l, fields const* f, params const* p, long i) {

	double res = 0.0;
	double S = f->singlet[i][0];
	/* kinetic term: \sum_{x,i} [S(x)^2 - S(x)S(x+i)] */
	res += l->dim * S*S;

	for (int dir=0; dir<l->dim; dir++) {
		long next = l->next[i][dir];
		long prev = l->prev[i][dir];

		res -= S * (f->singlet[next][0] + f->singlet[prev][0]);
	}

	res += potential_singlet(f, p, i);

	return res;
}

/* Contributions to the scalar potential from singlet field */
double potential_singlet(fields const* f, params const* p, long i) {

	double pot = 0.0;
	/* V(S) = b1 S + 1/2 msq_s S^2 + 1/3 b3 S^3 + 1/4 b4 S^4 + 1/2 a1 S \he\phi\phi + 1/2 a2 S^2 \he\phi\phi */
	double S = f->singlet[i][0];
	pot += p->b1_s * S + 0.5*p->msq_s * S*S + 1.0/3.0 * p->b3_s * S*S*S + 0.25*p->b4_s * S*S*S*S;

	#if (NHIGGS==1)
		double *h1 = f->su2doublet[0][i];
		double mod = doubletsq(h1);

		pot += 0.5*p->a1_s * S * mod + 0.5*p->a2_s * S*S * mod;
	#endif

	return pot;
}

#endif // end SINGLET routines

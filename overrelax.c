/** @file overrelax.c
*
* Routines for implementing the overrelaxation algorithm for fields.
*	Actual updating of the lattice is performed in update.c
*
* TODO figure out why Higgs overrelaxation blows up if lambda is very small
*
*/

#include "su2.h"


// Find real root to the equation a*x^3 + b*x^2 + c*x + d = 0
double polysolve3(long double a, long double b, long double c, long double d) {
	// discriminant. NB! the numbers here can get quite large/small...
	long double discr = 18.0*a*b*c*d - 4.0*b*b*b*d + b*b*c*c - 4.0*a*c*c*c - 27.0*a*a*d*d;
	long double del1 =  2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d;
	double x;

	if (discr < 0.0) {
		// there is only one real root and the solution is easy
		long double A = 0.5 * (del1 + a * sqrt(-27.0 * discr));
		A = cbrt(A); // ~ O(1) number

		x = -1.0 * (b + A + (b*b - 3.0*a*c) / A) / (3.0 * a);
	} else {
		// now there exist more than 1 real root, and the solution is more complicated.
		// this should not happen in realistic cases! so just return the trivial solution.
		// note that our b/a = old X
		x = b / a;
		printf("Warning: Did not find solution to X overrelax!\n");
	}
	return x;
}


#if (NHIGGS > 0)


/* Higgs overrelaxation as described in hep-lat/9510020 and hep-lat/9804019 (more recent).
* The doublet is Phi(x) = 1/sqrt(2) (phi_0 I + i sig_i phi_i), i = 1,2,3,
* and the local Higgs action is written as
* 	S = phi_a F_a + B phi_a phi_a + C (phi_a phi_a)^2. (note sign of F_a! Kari has diff.)
* Then define:
*		F = sqrt(F_a F_a), f_a = F_a/F, X = -phi_a f_a, Y = phi_a + X f_a (note signs again)
*	The action becomes:
*		S = X F + B (Y^2 + X^2) + C (Y^4 + X^4 + 2 X^2 Y^2).
*	We then overrelax X,Y individually: Y -> -Y
* and new X' is solved from S(X') - S(X) = 0 that is accepted with probability
* p(X') = min(p0, 1), p0 = (dS(X)/dX) / (dS(X')/dX'). If new X is accepted,
* the Y overrelaxation reads: phi'_a = -phi_a - f_a (X' + X).
* Note however that the Y overrelaxation is not necessary at all, X update is enough.
* So if the action is not invariant under Y -> -Y, can choose to keep Y constant instead. */
int overrelax_doublet(lattice const* l, fields* f, params const* p, long i) {

	// 1 Higgs doublet only!! For 2 Higgses see overrelax_higgs2()
	double s[SU2DB] = {0.0};
	int higgs_id = 0;

	double** higgs = f->su2doublet[higgs_id];
	double mod = doubletsq(higgs[i]);

	// calculate hopping staple s_a (denote s_a = F_a)
	staple_doublet(s, l, f, p, i, higgs_id);
	// staple normalization
	double F = 0.0;
	for (int k=0; k<SU2DB; k++) {
		F += s[k] * s[k];
	}
	F = sqrt(F);

	// Cartesian X and Y coordinates for the Higgs.
	double X = 0.0;
	for (int k=0; k<SU2DB; k++) {
		X += higgs[i][k] * s[k]; // this "contains" a minus sign
	}
	X /= F;
	double Y[4], Ysq = 0.0;
	for (int k=0; k<SU2DB; k++) {
		Y[k] = higgs[i][k] - X * s[k] / F;
		Ysq += Y[k] * Y[k];
	}

	// remaining terms in the local action
	double B = 0.5 * p->msq_phi + 1.0 * l->dim;
	#ifdef TRIPLET
		B += 0.5 * p->a2 * tripletsq(f->su2triplet[i]);
	#endif
	#ifdef SINGLET
		// V = 1/2 a1 S \he\phi\phi + 1/2 a2 S^2 \he\phi\phi + ...
		double S = f->singlet[i][0];
		B += 0.25 * p->a1_s * S + 0.25 * p->a2_s * S*S;
	#endif

	double C = 0.25 * p->lambda_phi;

	// we need to solve V(X') - V(X) = 0, where Y is kept constant. Write this as
	// (x - y) (alpha x^3 + beta * x^2 + gamma * x + delta) = 0
	// for x = X', y = X, and find a nontrivial real root.
	// Note: need to be careful about large numbers here. Best to not rescale by C
	long double alpha = C;
	long double beta = C*X;
	long double cc = B + 2.0 * Ysq * C; // quadratic term coefficient in V
	long double gamma = (cc + C * X*X);
	long double delta = (F + cc * X + C * X*X*X);

	double newX = polysolve3(alpha, beta, gamma, delta);
	// now accept/reject based on the derivatives
	// dV/dX = F + 2B*X + 4C(X^3 + Y^2 X), with again same Y in both cases
	double dV = F + 2.0*B*X + 4.0*C*(X*X*X + Ysq*X);
	double dV_new = F + 2.0*B*newX + 4.0*C*(newX*newX*newX + Ysq*newX);

	beta = fabs(dV/dV_new);

	if (beta >= drand48()) {
		// accept, so overrelax Y' = -Y using the new X
		// phi'_a = Y' - X' f_a = -phi_a - (X' + X) f_a,
		// but my X calculated above has diff sign already. */
		for (int k=0; k<SU2DB; k++) {
			f->su2doublet[0][i][k] = -1.0*higgs[i][k] + (newX + X) * s[k] / F;
		}
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}

#endif // NHIGGS > 0


#ifdef TRIPLET

/* Same Cartesian overrelax as overrelax_doublet(), but for adjoint scalar.
*/
int overrelax_triplet(lattice const* l, fields* f, params const* p, long i) {

	double s[3] = {0.0, 0.0, 0.0};
	double u[SU2LINK];
	double b[SU2TRIP];

	// calculate hopping staple s_a (denote. s_a = F_a)
	for (int dir=0; dir<l->dim; dir++) {
		long next = l->next[i][dir];
		// link variable
		for (int d=0; d<SU2LINK; d++) {
			u[d] = f->su2link[i][dir][d];
		}
		// Sigma at next site
		for (int d=0; d<SU2TRIP; d++) {
			b[d] = f->su2triplet[next][d];
		}
		s[0] += -(b[0]*(u[0]*u[0])) - b[0]*(u[1]*u[1]) + 2*b[2]*u[0]*u[2] -
		 				2*b[1]*u[1]*u[2] + b[0]*(u[2]*u[2]) - 2*b[1]*u[0]*u[3]
						- 2*b[2]*u[1]*u[3] + b[0]*(u[3]*u[3]);
		s[1] += -(b[1]*(u[0]*u[0])) - 2*b[2]*u[0]*u[1] + b[1]*(u[1]*u[1]) -
						2*b[0]*u[1]*u[2] - b[1]*(u[2]*u[2]) + 2*b[0]*u[0]*u[3] -
						2*b[2]*u[2]*u[3] + b[1]*(u[3]*u[3]);
		s[2] += -(b[2]*(u[0]*u[0])) + 2*b[1]*u[0]*u[1] + b[2]*(u[1]*u[1])
						- 2*b[0]*u[0]*u[2] + b[2]*(u[2]*u[2]) - 2*b[0]*u[1]*u[3]
						- 2*b[1]*u[2]*u[3] - b[2]*(u[3]*u[3]);
		// same for backwards directions
		long prev = l->prev[i][dir];
		for (int d=0; d<SU2LINK; d++) {
			u[d] = f->su2link[prev][dir][d];
		}
		for (int d=0; d<SU2TRIP; d++) {
			b[d] = f->su2triplet[prev][d];
		}
		s[0] += -(b[0]*(u[0]*u[0])) - b[0]*(u[1]*u[1]) - 2*b[2]*u[0]*u[2]
						- 2*b[1]*u[1]*u[2] + b[0]*(u[2]*u[2]) + 2*b[1]*u[0]*u[3]
						- 2*b[2]*u[1]*u[3] + b[0]*(u[3]*u[3]);
		s[1] += -(b[1]*(u[0]*u[0])) + 2*b[2]*u[0]*u[1] + b[1]*(u[1]*u[1])
						- 2*b[0]*u[1]*u[2] - b[1]*(u[2]*u[2]) - 2*b[0]*u[0]*u[3]
						- 2*b[2]*u[2]*u[3] + b[1]*(u[3]*u[3]);
		s[2] += -(b[2]*(u[0]*u[0])) - 2*b[1]*u[0]*u[1] + b[2]*(u[1]*u[1])
						+ 2*b[0]*u[0]*u[2] + b[2]*(u[2]*u[2]) - 2*b[0]*u[1]*u[3]
						- 2*b[1]*u[2]*u[3] - b[2]*(u[3]*u[3]);
	}
	// staple normalization
	double F = 0.0;
	for (int k=0; k<SU2TRIP; k++) {
		F += s[k] * s[k];
	}
	F = sqrt(F);

	// Cartesian X and Y coordinates
	double X = 0.0;
	for (int k=0; k<SU2TRIP; k++) {
		X += f->su2triplet[i][k] * s[k];
	}
	X /= F;
	double Y[SU2TRIP], Ysq = 0.0;
	for (int k=0; k<SU2TRIP; k++) {
		Y[k] = f->su2triplet[i][k] - X * s[k] / F;
		Ysq += Y[k] * Y[k];
	}

	// remaining terms in the local action
	double B = 0.5 * p->msq_triplet + 1.0 * l->dim;
	#ifdef HIGGS
		B += 0.5 * p->a2 * doubletsq(f->su2doublet[i]);
	#endif
	double C = 0.25 * p->b4;

	// we need to solve V(X') - V(X) = 0, where Y is kept constant. Write this as
	// (x - y) (alpha x^3 + beta * x^2 + gamma * x + delta) = 0
	// for x = X', y = X, and find a nontrivial real root.
	// Note: need to be careful about large numbers here. Best to not rescale by C
	long double alpha = C;
	long double beta = C*X;
	long double cc = B + 2.0 * Ysq * C; // quadratic term coefficient in V
	long double gamma = (cc + C * X*X);
	long double delta = (F + cc * X + C * X*X*X);

	double newX = polysolve3(alpha, beta, gamma, delta);
	// now accept/reject based on the derivatives
	// dV/dX = F + 2B*X + 4C(X^3 + Y^2 X), with again same Y in both cases
	double dV = F + 2.0*B*X + 4.0*C*(X*X*X + Ysq*X);
	double dV_new = F + 2.0*B*newX + 4.0*C*(newX*newX*newX + Ysq*newX);

	beta = fabs(dV/dV_new);

	if (beta >= drand48()) {
		// accept, so overrelax Y' = -Y using the new X
		for (int k=0; k<SU2TRIP; k++) {
			f->su2triplet[i][k] = -1.0*f->su2triplet[i][k] + (newX + X) * s[k] / F;
		}
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}

#endif // TRIPLET


#ifdef SINGLET

int overrelax_singlet(lattice const* l, fields* f, params const* p, long i) {

	double S = f->singlet[i][0];
	/* Local action due to S(x): act = c1 S + c2 S^2 + c3 S^3 + c4 S^4 */
	double c1 = p->b1_s;
	for (int dir=0; dir<l->dim; dir++) {
		long next = l->next[i][dir];
		long prev = l->prev[i][dir];
		c1 -= (f->singlet[next][0] + f->singlet[prev][0]);
	}

	double c2 = l->dim + 0.5*p->msq_s;
	double c3 = 1.0/3.0 * p->b3_s;
	double c4 = 0.25 * p->b4_s;

	#if (NHIGGS == 1)
		double phisq = doubletsq(f->su2doublet[0][i]);
		c1 += 0.5*p->a1_s*phisq;
		c2 += 0.5*p->a2_s*phisq;
	#endif

	/* Solve Y from act(S) - act(Y) = 0. Can factor out trivial solution Y=S, so
	* need to solve d0 + d1 Y + d2 Y^2 + d3 Y^3 = 0 */
	double d0 = c1 + c2*S + c3*S*S + c4*S*S*S;
	double d1 = c2 + c3*S + c4*S*S;
	double d2 = c3 + c4*S;
	double d3 = c4;

	double Y = polysolve3(d3, d2, d1, d0);

	/* Acc/rej to preserve detailed balance. Derivatives of the local action wrt. S and Y */
	double dV = c1 + 2.0*c2*S + 3.0*c3*S*S + 4.0*c4*S*S*S;
	double dV_new = c1 + 2.0*c2*Y + 3.0*c3*Y*Y + 4.0*c4*Y*Y*Y;

	if (fabs(dV/dV_new) >= drand48()) {
		f->singlet[i][0] = Y;
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}

#endif

#if (NHIGGS == 2)

/* XY-overrelax for multiple Higgs doublets. higgs_id specifies which doublet
* is updated. Here I am more consistent with signs and take the local action to be
* V[phi] ~ -F_a phi_a + others. Also, Y is not changed at all here, so the update works
* even for potentials that are not invariant under Y -> -Y. */
int overrelax_higgs2(lattice const* l, fields* f, params const* p, long i, int higgs_id) {

	double s[SU2DB] = {0.0};

	double** higgs = f->su2doublet[higgs_id];
	double mod = doubletsq(higgs[i]);

	// calculate hopping staple s_a (denote s_a = F_a)
	staple_doublet(s, l, f, p, i, higgs_id); // now S ~ f[a] s[a], does not include minus sign

	/* Now there is f12 = phi1^+ phi2 etc in the potential, add these to the "staple".
	* There are also f1[a] * f1^2 etc, but these do NOT contribute to F_a which should
	* remain constant in the local f1[a] -> f_new[a] update.
	* Here R = Re f12, I = Im f12, H = phi_other (4-vec), G = +/- i sig_2 H (4-dimensional Pauli matrix),
	* so in terms of f1[a] vectors: f11 = 0.5 * f1.f1, R = 0.5 * f1.f2, I = 0.5 * f1.G.
	* The sign in G is - if updating phi1 (so H = phi2) and + if updating phi2. */
	double* H;
	double G[SU2DB];
	int sign = -1; // sign of G

	// full staple is s[a] = s_hop[a] + A H + B G. m12^2 term and either lam6 or lam7 contribute
	complex lam67, lam67other;
	double lam12, msq;
	if (higgs_id == 0) {

		H = f->su2doublet[1][i];

		lam12 = p->lambda_phi;
		msq = p->msq_phi;
		lam67 = p->lam7;
		lam67other = p->lam6;
		lam67other.im *= -1.0;
	} else {

		H = f->su2doublet[0][i];
		sign = 1;

		lam12 = p->lam2;
		msq = p->msq_phi2;
		lam67 = p->lam6;
		lam67.im *= -1.0;
		lam67other = p->lam7;
	}

	G[0] = sign*H[3]; G[1] = sign*H[2];
	G[2] = -1.0*sign*H[1]; G[3] = -1.0*sign*H[0];

	double Hsq = doubletsq(H);

	for (int k=0; k<SU2DB; k++) {
		s[k] += 0.5*(p->m12sq.re + Hsq * lam67.re) * H[k];
		s[k] += 0.5*(-1.0*p->m12sq.im + Hsq * lam67.im) * G[k];
		s[k] = -1.0*s[k]; // change staple sign to match F_a
	}

	// staple normalization
	double F = 0.0;
	for (int k=0; k<SU2DB; k++) F += s[k] * s[k];
	F = sqrt(F);

	// Cartesian X and Y coordinates for the Higgs
	double X = 0.0;
	for (int k=0; k<SU2DB; k++) {
		s[k] /= F; // s <- s/F = f_a
		X += higgs[i][k] * s[k];
	}

	double Y[SU2DB], Ysq = 0.0;
	for (int k=0; k<SU2DB; k++) {
		Y[k] = higgs[i][k] - X * s[k];
		Ysq += Y[k] * Y[k];
	}

	/* Write local action as V(X) = b4 X^4 + b3 X^3 + b2 X^2 + b1 X + b0.
	* See Mathematica notebook overrelax.c for the expressions */

	// some dot products
	double Hf = 0.0, Gf = 0.0, HY = 0.0, GY = 0.0;
	for (int k=0; k<SU2DB; k++) {
		Hf += H[k] * s[k];
		Gf += G[k] * s[k];
		HY += H[k] * Y[k];
		GY += G[k] * Y[k];
	}

	// terms that differ for phi1 and phi2 (b2 includes term from covariant der.)
	double b4 = 0.25 * lam12;
	double b3 = 0.25 * (Hf * lam67other.re + Gf * lam67other.im);
	double b2 = 0.5 * msq + 1.0*l->dim + 0.5*Ysq*lam12 + 0.25*(HY*lam67other.re + GY*lam67other.im);
	double b1 = -1.0*F + 0.25*Ysq*(Hf*lam67other.re + Gf*lam67other.im); // F contains terms from the "staple"

	// then mutual terms for both phi1,2
	b2 += 0.5*Hsq * p->lam3 + 0.25*p->lam4 * (Hf*Hf + Gf*Gf)
			+ 0.25*p->lam5.re * (Hf*Hf - Gf*Gf) - 0.5*p->lam5.im * Gf*Hf;

	b1 += 0.5*Hf * (HY*p->lam4 + HY*p->lam5.re - GY*p->lam5.im)
			+ 0.5*Gf * (GY*p->lam4 - GY*p->lam5.re - HY*p->lam5.im);

	/* Solve V(X') - V(X) = 0: write this as (x - X) (ax^3 + bx^2 + cx + d) = 0
	*  and find the real root of the 3rd degree polynomial */
	long double aa = b4;
	long double bb = b3 + b4*X;
	long double cc = b2 + b3*X + b4*X*X;
	long double dd = b1 + b2*X + b3*X*X + b4*X*X*X;

	double Xn = polysolve3(aa, bb, cc, dd);
	/* now accept/reject based on the change in measure */
	double dV = 4.0*b4*X*X*X + 3.0*b3*X*X + 2.0*b2*X + b1;
	double dV_new = 4.0*b4*Xn*Xn*Xn + 3.0*b3*Xn*Xn + 2.0*b2*Xn + b1;

	b4 = fabs(dV/dV_new);
	if (b4 >= drand48()) {
		// accept, so change phi_a so that Y is unchanged.
		for (int k=0; k<SU2DB; k++) {
			higgs[i][k] = higgs[i][k] + s[k] * (Xn - X);
		}
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}
#endif // NHIGGS == 2

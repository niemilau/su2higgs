/** @file overrelax.c
*
* Routines for implementing the overrelaxation algorithm for fields.
*	Actual updating of the lattice is performed in update.c
*
* TODO figure out why Higgs overrelaxation blows up if lambda is very small
*
*/

#include "su2.h"


// Find real root to the equation x^3 + b*x^2 + c*x + d = 0
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


/* Higgs overrelaxation as described in hep-lat/9510020 and hep-lat/9804019 (more recent).
* The doublet is Phi(x) = 1/sqrt(2) (phi_0 I + i sig_i phi_i), i = 1,2,3,
* and the local Higgs action is written as
* 	S = phi_a F_a + B phi_a phi_a + C (phi_a phi_a)^2.
* Then define:
*		F = sqrt(F_a F_a), f_a = F_a/F, X = phi_a f_a, Y = phi_a - X f_a
*	The action becomes:
*		S = X F + B (Y^2 + X^2) + C (Y^4 + X^4 + 2 X^2 Y^2).
*	We then overrelax X,Y individually: Y -> -Y
* and new X' is solved from S(X') - S(X) = 0 that is accepted with probability
* p(X') = min(p0, 1), p0 = (dS(X)/dX) / (dS(X')/dX'). If new X is accepted,
* the Y overrelaxation reads: phi'_a = -phi_a + f_a (X' + X).
* (could probably just ignore the Y reflection a la Kari? So that
* 	phi'_a = Y + X' f_a )
*
*/
int overrelax_doublet(fields f, params p, long i) {

	double s[4] = {0.0, 0.0, 0.0, 0.0};
	double u[4];
	double b[4];

	// calculate hopping staple s_a (denote. s_a = F_a)
	for (int dir=0; dir<p.dim; dir++) {
		u[0] = f.su2link[i][dir][0];
		u[1] = f.su2link[i][dir][1];
		u[2] = f.su2link[i][dir][2];
		u[3] = f.su2link[i][dir][3];
		long next = p.next[i][dir];
		// Phi at next site
		b[0] = f.su2doublet[next][0];
		b[1] = f.su2doublet[next][1];
		b[2] = f.su2doublet[next][2];
		b[3] = f.su2doublet[next][3];
		s[0] += -(b[0]*u[0]) + b[1]*u[1] + b[2]*u[2] + b[3]*u[3];
		s[1] += -(b[1]*u[0]) - b[0]*u[1] + b[3]*u[2] - b[2]*u[3];
		s[2] += -(b[2]*u[0]) - b[3]*u[1] - b[0]*u[2] + b[1]*u[3];
		s[3] += -(b[3]*u[0]) + b[2]*u[1] - b[1]*u[2] - b[0]*u[3];
		// same for backwards directions
		long prev = p.prev[i][dir];
		u[0] = f.su2link[prev][dir][0];
		u[1] = f.su2link[prev][dir][1];
		u[2] = f.su2link[prev][dir][2];
		u[3] = f.su2link[prev][dir][3];
		// Phi at previous site
		b[0] = f.su2doublet[prev][0];
		b[1] = f.su2doublet[prev][1];
		b[2] = f.su2doublet[prev][2];
		b[3] = f.su2doublet[prev][3];
		s[0] += -(b[0]*u[0]) - b[1]*u[1] - b[2]*u[2] - b[3]*u[3];
		s[1] += -(b[1]*u[0]) + b[0]*u[1] - b[3]*u[2] + b[2]*u[3];
		s[2] += -(b[2]*u[0]) + b[3]*u[1] + b[0]*u[2] - b[1]*u[3];
		s[3] += -(b[3]*u[0]) - b[2]*u[1] + b[1]*u[2] + b[0]*u[3];
	}
	// staple normalization
	double F = 0.0;
	for (int k=0; k<SU2DB; k++) {
		F += s[k] * s[k];
	}
	F = sqrt(F);

	// Cartesian X and Y coordinates for the Higgs.
	double X = 0.0;
	for (int k=0; k<SU2DB; k++) {
		X += f.su2doublet[i][k] * s[k];
	}
	X /= F;
	double Y[4], Ysq = 0.0;
	for (int k=0; k<SU2DB; k++) {
		Y[k] = f.su2doublet[i][k] - X * s[k] / F;
		Ysq += Y[k] * Y[k];
	}

	// remaining terms in the local action
	double B = 0.5 * p.msq_phi + 1.0 * p.dim;
	#ifdef TRIPLET
		B += 0.5 * p.a2 * tripletsq(f.su2triplet[i]);
	#endif
	double C = 0.25 * p.lambda_phi;

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
		for (int k=0; k<SU2DB; k++) {
			f.su2doublet[i][k] = -1.0*f.su2doublet[i][k] + (newX + X) * s[k] / F; // this works well
			//f.su2doublet[i][k] = Y[k] + newX * s[k] / F; // Kari seems to have this instead, no Y -> -Y?
		}
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}


/* Same Cartesian overrelax as overrelax_doublet(), but for adjoint scalar.
*/
int overrelax_triplet(fields f, params p, long i) {

	double s[3] = {0.0, 0.0, 0.0};
	double u[SU2LINK];
	double b[SU2TRIP];

	// calculate hopping staple s_a (denote. s_a = F_a)
	for (int dir=0; dir<p.dim; dir++) {
		long next = p.next[i][dir];
		// link variable
		for (int d=0; d<SU2LINK; d++) {
			u[d] = f.su2link[i][dir][d];
		}
		// Sigma at next site
		for (int d=0; d<SU2TRIP; d++) {
			b[d] = f.su2triplet[next][d];
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
		long prev = p.prev[i][dir];
		for (int d=0; d<SU2LINK; d++) {
			u[d] = f.su2link[prev][dir][d];
		}
		for (int d=0; d<SU2TRIP; d++) {
			b[d] = f.su2triplet[prev][d];
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
		X += f.su2triplet[i][k] * s[k];
	}
	X /= F;
	double Y[SU2TRIP], Ysq = 0.0;
	for (int k=0; k<SU2TRIP; k++) {
		Y[k] = f.su2triplet[i][k] - X * s[k] / F;
		Ysq += Y[k] * Y[k];
	}

	// remaining terms in the local action
	double B = 0.5 * p.msq_triplet + 1.0 * p.dim;
	#ifdef HIGGS
		B += 0.5 * p.a2 * doubletsq(f.su2doublet[i]);
	#endif
	double C = 0.25 * p.b4;

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
			f.su2triplet[i][k] = -1.0*f.su2triplet[i][k] + (newX + X) * s[k] / F; // this works well
		}
		return 1;
	} else {
		// reject, no changes to the field
		return 0;
	}

}

/* OLD overrelaxation for the doublet.
*
*	The doublet is Phi(x) = 1/sqrt(2) sig_i phi_i , i = 0,1,2,3,
*	and we perform the update for phi(x)_a following Kari's simple recipe:
* Write the action in the form
* 	S = phi_a s_a + B phi_a phi_a + C (phi_a phi_a)^2,
* where s_a is a staple from the hopping terms. Completing the square gives:
* 	S = B (phi_a + s_a /(2B))^2 + C (phi_a phi_a)^2 - s_a^2 /(4B),
*	where the last term is constant in phi_a.
* phi_a is overrelaxed by choosing new p_a so that
*		(p_a + s_a / (2B) ) = -(phi_a + s_a/(2B)), i.e.
* 	p_a = -phi_a - s_a / B
*
*	Finally we apply a accept/reject step with the quartic term C as in metropolis.
* The acceptance rate should be high if C is small.
* Returns 1 if update was accepted and 0 if rejected.
*/
int overrelax_doublet_old(fields f, params p, long i) {

	double oldfield[SU2DB];
	for (int k=0; k<SU2DB; k++) {
		oldfield[k] = f.su2doublet[i][k];
	}
	double sa[4] = {0.0, 0.0, 0.0, 0.0};
	double u[4];
	double b[4];
	// calculate hopping staple s_a
	for (int dir=0; dir<p.dim; dir++) {
		u[0] = f.su2link[i][dir][0];
		u[1] = f.su2link[i][dir][1];
		u[2] = f.su2link[i][dir][2];
		u[3] = f.su2link[i][dir][3];
		long next = p.next[i][dir];
		// Phi at next site
		b[0] = f.su2doublet[next][0];
		b[1] = f.su2doublet[next][1];
		b[2] = f.su2doublet[next][2];
		b[3] = f.su2doublet[next][3];
		sa[0] += -(b[0]*u[0]) + b[1]*u[1] + b[2]*u[2] + b[3]*u[3];
		sa[1] += -(b[1]*u[0]) - b[0]*u[1] + b[3]*u[2] - b[2]*u[3];
		sa[2] += -(b[2]*u[0]) - b[3]*u[1] - b[0]*u[2] + b[1]*u[3];
		sa[3] += -(b[3]*u[0]) + b[2]*u[1] - b[1]*u[2] - b[0]*u[3];
		// same for backwards directions
		long prev = p.prev[i][dir];
		u[0] = f.su2link[prev][dir][0];
		u[1] = f.su2link[prev][dir][1];
		u[2] = f.su2link[prev][dir][2];
		u[3] = f.su2link[prev][dir][3];
		// Phi at previous site
		b[0] = f.su2doublet[prev][0];
		b[1] = f.su2doublet[prev][1];
		b[2] = f.su2doublet[prev][2];
		b[3] = f.su2doublet[prev][3];
		sa[0] += -(b[0]*u[0]) - b[1]*u[1] - b[2]*u[2] - b[3]*u[3];
		sa[1] += -(b[1]*u[0]) + b[0]*u[1] - b[3]*u[2] + b[2]*u[3];
		sa[2] += -(b[2]*u[0]) + b[3]*u[1] + b[0]*u[2] - b[1]*u[3];
		sa[3] += -(b[3]*u[0]) - b[2]*u[1] + b[1]*u[2] + b[0]*u[3];
	}

	// calculate coefficient of quadratic terms, B.
	// B = mass term + local term from covariant derivative + extra from other fields
	double B = 0.5 * p.msq_phi + 1.0 * p.dim;
	#ifdef TRIPLET
		B += 0.5 * p.a2 * tripletsq(f.su2triplet[i]);
	#endif
	// add check for B != 0?

	// calculate quartic term before update
	double mod = doubletsq(oldfield);
	double oldquartic = p.lambda_phi * mod * mod;

	// Overrelax all the components and then apply a accept/reject based on the quartic term.
	// This should be OK since overrelaxation for phi_a is independent of the other components
	for (int dof=0; dof<4; dof++) {
		f.su2doublet[i][dof] = -1.0*f.su2doublet[i][dof] - sa[dof] / B;
	}
	mod = doubletsq(f.su2doublet[i]);
	double newquartic = p.lambda_phi * mod * mod;

	double diff = newquartic - oldquartic;

	int accept;

	if (diff < 0) {
		accept = 1;
	}
	else if (diff >= 0 && ( exp(-(diff)) > drand48() )) {
		accept = 1;
	}
	else {
		accept = 0;
		for (int dof=0; dof<SU2DB; dof++) {
			f.su2doublet[i][dof] = oldfield[dof];
		}
	}

	return accept;
}


/* Same as overrelax_doublet_old() but for the real triplet.
* Logic is the same: reflect the Gaussian part, acc/rej on the quartic part.
*/
int overrelax_triplet_old(fields f, params p, long i) {

	double oldfield[3];
	oldfield[0] = f.su2triplet[i][0];
	oldfield[1] = f.su2triplet[i][1];
	oldfield[2] = f.su2triplet[i][2];

	double sa[3] = {0.0, 0.0, 0.0};
	double u[4];
	double b[3];
	// calculate hopping staple s_a
	for (int dir=0; dir<p.dim; dir++) {
		u[0] = f.su2link[i][dir][0];
		u[1] = f.su2link[i][dir][1];
		u[2] = f.su2link[i][dir][2];
		u[3] = f.su2link[i][dir][3];
		long next = p.next[i][dir];
		// triplet at next site
		b[0] = f.su2triplet[next][0];
		b[1] = f.su2triplet[next][1];
		b[2] = f.su2triplet[next][2];
		sa[0] += -(b[0]*(u[0]*u[0])) - b[0]*(u[1]*u[1]) + 2*b[2]*u[0]*u[2] - 2*b[1]*u[1]*u[2] +
							b[0]*(u[2]*u[2]) - 2*b[1]*u[0]*u[3] - 2*b[2]*u[1]*u[3] + b[0]*(u[3]*u[3]);
		sa[1] += -(b[1]*(u[0]*u[0])) - 2*b[2]*u[0]*u[1] + b[1]*(u[1]*u[1]) - 2*b[0]*u[1]*u[2] -
   						b[1]*(u[2]*u[2]) + 2*b[0]*u[0]*u[3] - 2*b[2]*u[2]*u[3] + b[1]*(u[3]*u[3]);
		sa[2] += -(b[2]*(u[0]*u[0])) + 2*b[1]*u[0]*u[1] + b[2]*(u[1]*u[1]) - 2*b[0]*u[0]*u[2] +
   						b[2]*(u[2]*u[2]) - 2*b[0]*u[1]*u[3] - 2*b[1]*u[2]*u[3] - b[2]*(u[3]*u[3]);
		// same for backwards directions
		long prev = p.prev[i][dir];
		u[0] = f.su2link[prev][dir][0];
		u[1] = f.su2link[prev][dir][1];
		u[2] = f.su2link[prev][dir][2];
		u[3] = f.su2link[prev][dir][3];
		// Phi at previous site
		b[0] = f.su2triplet[prev][0];
		b[1] = f.su2triplet[prev][1];
		b[2] = f.su2triplet[prev][2];
		sa[0] += -(b[0]*(u[0]*u[0])) - b[0]*(u[1]*u[1]) - 2*b[2]*u[0]*u[2] - 2*b[1]*u[1]*u[2] +
   						b[0]*(u[2]*u[2]) + 2*b[1]*u[0]*u[3] - 2*b[2]*u[1]*u[3] + b[0]*(u[3]*u[3]);
		sa[1] += -(b[1]*(u[0]*u[0])) + 2*b[2]*u[0]*u[1] + b[1]*(u[1]*u[1]) - 2*b[0]*u[1]*u[2] -
   						b[1]*(u[2]*u[2]) - 2*b[0]*u[0]*u[3] - 2*b[2]*u[2]*u[3] + b[1]*(u[3]*u[3]);
		sa[2] += -(b[2]*(u[0]*u[0])) - 2*b[1]*u[0]*u[1] + b[2]*(u[1]*u[1]) + 2*b[0]*u[0]*u[2] +
   						b[2]*(u[2]*u[2]) - 2*b[0]*u[1]*u[3] - 2*b[1]*u[2]*u[3] - b[2]*(u[3]*u[3]);
	}

	// calculate coefficient of quadratic terms, B.
	// B = mass term + local term from covariant derivative + extra from other fields
	double B = 0.5 * p.msq_triplet + 1.0 * p.dim;
	#ifdef HIGGS
		B += 0.5 * p.a2 * doubletsq(f.su2doublet[i]);
	#endif
	// add check for B != 0?

	// calculate quartic term before update
	double mod = tripletsq(oldfield);
	double oldquartic = p.b4 * mod * mod;

	// Overrelax all the components and then apply a accept/reject based on the quartic term.
	// This should be OK since overrelaxation for A_a is independent of the other components
	for (int dof=0; dof<3; dof++) {
		f.su2triplet[i][dof] = -1.0*f.su2triplet[i][dof] - sa[dof] / B;
	}
	mod = tripletsq(f.su2triplet[i]);
	double newquartic = p.b4 * mod * mod;

	double diff = newquartic - oldquartic;

	if (diff < 0) {
		return 1;
	}
	else if (diff >= 0 && ( exp(-(diff)) > drand48() )) {
		return 1;
	}
	else {
		f.su2triplet[i][0] = oldfield[0];
		f.su2triplet[i][1] = oldfield[1];
		f.su2triplet[i][2] = oldfield[2];
		return 0;
	}
}

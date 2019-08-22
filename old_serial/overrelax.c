/** @file overrelax.c
*
* Routines for implementing the overrelaxation algorithm for fields.
*	Actual updating of the lattice is performed in update.c
*
* TODO figure out why Higgs overrelaxation blows up if lambda is very small
*
*/

#include "su2.h"


/*
* Update a single SU(2) scalar doublet using overrelaxation.
*	The doublet is Phi(x) = 1/sqrt(2) sig_i phi_i , i = 0,1,2,3,
*	and we perform the update for phi(x)_a following Kari's recipe:
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
int overrelax_doublet(fields f, params p, ulong i) {

	double oldfield[4];
	oldfield[0] = f.su2doublet[i][0];
	oldfield[1] = f.su2doublet[i][1];
	oldfield[2] = f.su2doublet[i][2];
	oldfield[3] = f.su2doublet[i][3];

	double sa[4] = {0.0, 0.0, 0.0, 0.0};
	double u[4];
	double b[4];
	// calculate hopping staple s_a
	for (int dir=0; dir<p.dim; dir++) {
		u[0] = f.su2link[i][dir][0];
		u[1] = f.su2link[i][dir][1];
		u[2] = f.su2link[i][dir][2];
		u[3] = f.su2link[i][dir][3];
		ulong next = p.next[i][dir];
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
		ulong prev = p.prev[i][dir];
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

	//debug
	//oldquartic = localact_doublet(f,p,i);

	// Overrelax all the components and then apply a accept/reject based on the quartic term.
	// This should be OK since overrelaxation for phi_a is independent of the other components
	for (int dof=0; dof<4; dof++) {
		f.su2doublet[i][dof] = -1.0*f.su2doublet[i][dof] - sa[dof] / B;
	}
	mod = doubletsq(f.su2doublet[i]);
	double newquartic = p.lambda_phi * mod * mod;
	//debug
	//newquartic = localact_doublet(f, p, i);

	double diff = newquartic - oldquartic;


	if (diff < 0) {
		return 1;
	}
	else if (diff >= 0 && ( exp(-(diff)) > drand48() )) {
		return 1;
	}
	else {
		f.su2doublet[i][0] = oldfield[0];
		f.su2doublet[i][1] = oldfield[1];
		f.su2doublet[i][2] = oldfield[2];
		f.su2doublet[i][3] = oldfield[3];
		return 0;
	}
}

/* Same as overrelax_doublet() but for the real triplet.
* Logic is the same: reflect the Gaussian part, acc/rej on the quartic part.
*/
int overrelax_triplet(fields f, params p, ulong i) {

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
		ulong next = p.next[i][dir];
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
		ulong prev = p.prev[i][dir];
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

	//debug
	//oldquartic = localact_doublet(f,p,i);

	// Overrelax all the components and then apply a accept/reject based on the quartic term.
	// This should be OK since overrelaxation for A_a is independent of the other components
	for (int dof=0; dof<3; dof++) {
		f.su2triplet[i][dof] = -1.0*f.su2triplet[i][dof] - sa[dof] / B;
	}
	mod = tripletsq(f.su2triplet[i]);
	double newquartic = p.b4 * mod * mod;
	//debug
	//newquartic = localact_doublet(f, p, i);

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

/** @file init.c
*
* Routines for initializing a simulation.


//TODO
*
*/

#include "su2.h"


/*
****************** Field initializations ******************
*/

// set SU(2) links to unity
void setsu2(fields* f, lattice const* l) {

	for (long i=0; i<l->sites; i++) {
		for (int dir=0; dir<l->dim; dir++) {
			f->su2link[i][dir][0] = 1.0;
			f->su2link[i][dir][1] = 0.0;
			f->su2link[i][dir][2] = 0.0;
			f->su2link[i][dir][3] = 0.0;
		}
	}
}

// set U(1) links to unity
void setu1(fields* f, lattice const* l) {

	for (long i=0; i<l->sites; i++) {
		for (int dir=0; dir<l->dim; dir++) {
			f->u1link[i][dir] = 0.0;
		}
	}
}

/* generate a random SU(2) matrix and store it in the argument */
void random_su2link(double *su2) {

	double u[4];
	// this number can be adjusted for e.g. metro update
	u[0] = 10.0 * (dran() - 0.5);
	u[1] = 1.0 * (dran() - 0.5);
	u[2] = 1.0 * (dran() - 0.5);
	u[3] = 1.0 * (dran() - 0.5);

	double norm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);

	su2[0] = u[0]/norm;
	su2[1] = u[1]/norm;
	su2[2] = u[2]/norm;
	su2[3] = u[3]/norm;
}

// Initialize SU(2) doublets to unity * 1/sqrt(2).
void setdoublets(fields* f, lattice const* l, params const* p) {

	#if (NHIGGS > 0)

	for (int db=0; db<NHIGGS; db++) {
		for (long i=0; i<l->sites_total; i++) {
			f->su2doublet[db][i][0] = p->phi0;
			for (int dof=1; dof<SU2DB; dof++) f->su2doublet[db][i][dof] = 0.0;
		}
	}

	#endif
}


// Initialize SU(2) adjoint scalars
void settriplets(fields* f, lattice const* l, params const* p) {
	for (long i=0; i<l->sites; i++) {
		f->su2triplet[i][0] = p->sigma0;
		f->su2triplet[i][1] = 0.0;
		f->su2triplet[i][2] = 0.0;
	}
}

#ifdef SINGLET
	void set_singlets(fields* f, lattice const* l, params const* p) {
		for (long i=0; i<l->sites; i++) f->singlet[i][0] = p->singlet0;
	}
#endif

// Initialize all fields used in the simulation
void setfields(fields* f, lattice* l, params const* p) {

	setsu2(f, l);
	#ifdef U1
		setu1(f, l);
	#endif

	#if (NHIGGS > 0)
		setdoublets(f, l, p);
	#endif

	#ifdef TRIPLET
		settriplets(f, l, p);
	#endif

	#ifdef SINGLET
		set_singlets(f, l, p);
	#endif

	sync_halos(l, f);
}


/* Copy all fields from "fields" struct f_old to f_new.
* Used when e.g. performing Wilson flow for renormalization
* without wanting to lose the equilibrium configuration.
* Does not modify the old fields struct */
void copy_fields(lattice const* l, fields const* f_old, fields* f_new) {

	for (long i=0; i<l->sites_total; i++) {

		// gauge links
		for (int dir=0; dir<l->dim; dir++) {
			memcpy(f_new->su2link[i][dir], f_old->su2link[i][dir], SU2LINK*sizeof(f_old->su2link[i][dir][0]));

			#ifdef U1
				f_new->u1link[i][dir] = f_old->u1link[i][dir];
			#endif
		}

		// scalars
		#if (NHIGGS > 0 )
			for (int db=0; db<NHIGGS; db++) {
				 memcpy(f_new->su2doublet[db][i], f_old->su2doublet[db][i], SU2DB*sizeof(f_old->su2doublet[db][i][0]));
			 }
		#endif

		#ifdef TRIPLET
			memcpy(f_new->su2triplet[i], f_old->su2triplet[i], SU2TRIP*sizeof(f_old->su2triplet[i][0]));
		#endif

		#ifdef SINGLET
			f_new->singlet[i][0] = f_old->singlet[i][0];
		#endif

	}

}

/* Copy an individual field with 'dofs' components per site (does not copy halos) */
void cp_field(lattice const* l, double** field, double** new, int dofs, int parity) {

	long max, offset = 0;
	if (parity == EVENODD) max = l->sites;
	else if (parity == EVEN) max = l->evensites;
	else if (parity == ODD) {
		max = l->sites;
		offset = l->evensites;
	}

	for (long i=offset; i<max; i++) memcpy(new[i], field[i], dofs * sizeof(double));
}


// Initialize accept/reject counters and time
void init_counters(counters* c) {

	c->iter = 1;

	c->higgs_sweeps = 0;
	c->triplet_sweeps = 0;

	#ifdef SINGLET
		c->singlet_sweeps = 0;
	#endif

	c->accepted_su2link = 0;
	c->accepted_u1link = 0;

	c->accepted_triplet = 0;
	c->acc_overrelax_triplet = 0;

	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) {
			c->accepted_doublet[db] = 0;
			c->acc_overrelax_doublet[db] = 0;
			c->total_overrelax_doublet[db] = 0;
			c->total_doublet[db] = 0;
		}
	#endif

	#ifdef SINGLET
		c->accepted_singlet = 0;
		c->total_singlet = 0;
		c->acc_overrelax_singlet = 0;
		c->total_overrelax_singlet = 0;
	#endif

	c->total_su2link = 0;
	c->total_u1link = 0;
	c->total_triplet = 0;
	c->total_overrelax_triplet = 0;
	c->accepted_muca = 0;
	c->total_muca = 0;
}

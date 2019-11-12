
/** @file alloc.c
*
* Routines for memory allocation.

//TODO
*
*/

#include "su2.h"
#include "comms.h"
#ifdef WALL
  #include "wallprofile.h"
#endif

/* Allocates memory for a field with a single degree of freedom per site.
*
*/
double *make_singletfield(long sites) {

  double *field = malloc(sites * sizeof(*field));

	if (field == NULL) {
		printf("Failed to allocate memory for a singlet field!\n");
		die(12);
	}
  return field;
}


/* Allocates contiguous memory for a field with dofs degrees of freedom.
* It can be accessed as field[x][dof], where x is the lattice site.
*/
double **make_field(long sites, int dofs) {

	double **field = malloc(sites * sizeof(**field));

	if (field == NULL) {
		printf("Failed to allocate memory for a field!\n");
		die(11);
	}

	field[0] = malloc(sites * dofs * sizeof(field[0]));

	for (long i=0; i<sites; i++) {
		field[i] = field[0] + i * dofs;
	}

	return field;
}


/* Allocates contiguous memory for a gauge field with dofs degrees of freedom.
*	 It can be accessed as field[x][dir][dof], where x is the lattice site
*	 and dir is the direction, pointing towards positive coordinate axes.
*  Original implementation by David.
*/
double ***make_gaugefield(long sites, int dim, int dofs) {

	double *f = malloc(sites * dim * dofs * sizeof(*f));

  double ***field = malloc(sites * sizeof(***field));

	if (field == NULL || f == NULL) {
		printf("Failed to allocate memory for a gauge field!\n");
		die(10);
	}

	for (long i=0; i<sites; i++) {
		field[i] = malloc(dim * sizeof(field[i]));
		for (int dir=0; dir<dim; dir++) {
			field[i][dir] = &f[i * dim * dofs + dir * dofs];
		}
	}

	field[0][0] = f;

  return field;
}


/* Free the memory allocated by make_singletfield()
*
*/
void free_singletfield(double *field) {
  free(field);
}

/* Free the contiguous memory allocated by make_field()
*
*/
void free_field(double **field) {
	free(field[0]);
  free(field);
}

/* Free the contiguous memory allocated by make_gaugefield()
*
*/
void free_gaugefield(long sites, double ***field) {

	free(field[0][0]);

	for (long i=0; i<sites; i++) {
		free(field[i]);
	}

  free(field);
}


/** Allocate all the fields needed for a simulation.
 */
void alloc_fields(params p, fields *f) {

	long sites = p.sites_total;
	f->su2link = make_gaugefield(sites, p.dim, SU2LINK);

	#ifdef HIGGS
		f->su2doublet = make_field(sites, SU2DB);
	#endif
	#ifdef TRIPLET
		f->su2triplet = make_field(sites, SU2TRIP);
	#endif
  #ifdef U1
    /* allocate U(1) gauge field, accessed as u1link[i][dir]. Note that memory
    * wise this is essentially just a non-gauge field with p.dim components.
    */
    f->u1link = make_field(sites, p.dim);
  #endif
	if (!p.rank)
		printf("Allocated memory for fields.\n");
}


/** Free all the fields allocated by alloc_fields().
 */
void free_fields(params p, fields *f) {

	long sites = p.sites_total;
	free_gaugefield(sites, f->su2link);
	#ifdef HIGGS
		free_field(f->su2doublet);
	#endif
	#ifdef TRIPLET
		free_field(f->su2triplet);
	#endif
  #ifdef U1
    free_field(f->u1link);
  #endif
	if (!p.rank)
		printf("Freed memory allocated for fields.\n");

}


/* Allocate contiguous memory for long-valued 2D array.
* These are mainly used for navigating on the lattice
*/
long **alloc_latticetable(int dim, long sites) {

	long **array = malloc(sites * sizeof(*(array)));
	array[0] = malloc(sites * dim * sizeof(array[0]));
	for (long i=0; i<sites; i++) {
		array[i] = array[0] + i * dim;
	 }

  return array;
}


/* Allocate all needed lookup lists (both lattice sites and MPI nodes)
* also include halos
*/
void alloc_lattice_arrays(params *p) {

	p->next = alloc_latticetable(p->dim, p->sites_total);
	p->prev = alloc_latticetable(p->dim, p->sites_total);

	// allocate parity arrays
	p->parity = malloc(p->sites_total * sizeof(*(p->parity)));

	if (!p->rank)
		printf("Allocated memory for lookup tables.\n");
}


/* Free all memory reserved for layouting and lookup tables.
*
*/
void free_lattice_arrays(params *p) {

	// lookup tables and parity:
	free(p->parity);
	free_latticetable(p->next);
	free_latticetable(p->prev);
  #ifdef WALL
    free_latticetable(wallcoord);
  #endif
	// slicing:
	free(p->sliceL);
	free(p->nslices);
	// then finally the full lattice size, allocated in get_parameters()
	free(p->L);

	if (!p->rank)
		printf("Freed memory allocated for layouting.\n");
}

// Free memory allocated by allocate_latticetable()
void free_latticetable(long** list) {
	free(list[0]);
	free(list);
}

// Free memory allocated for comlist and its substructures
void free_comlist(comlist_struct* comlist) {

	for (int k=0; k<comlist->neighbors; k++) {
		free(comlist->recv_from[k].sitelist);
		free(comlist->send_to[k].sitelist);
	}
	free(comlist->recv_from);
	free(comlist->send_to);
}

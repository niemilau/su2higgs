
/** @file alloc.c
*
* Routines for memory allocation.

//TODO
*
*/

#include "su2.h"
#include "comms.h"

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
void alloc_fields(lattice const* l, fields *f) {

	long sites = l->sites_total;
	f->su2link = make_gaugefield(sites, l->dim, SU2LINK);

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
    f->u1link = make_field(sites, l->dim);
  #endif
}


/** Free all the fields allocated by alloc_fields().
 */
void free_fields(lattice const* l, fields *f) {

	long sites = l->sites_total;
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

}


/* Allocate contiguous memory for long-valued 2D array.
* These are mainly used for navigating on the lattice
*/
long **alloc_latticetable(int dim, long sites) {

	long **array = malloc(sites * sizeof(*(array)));
  if (array == NULL) {
    printf("Failed to allocate memory for a lattice table!\n");
    die(1001);
  }
	array[0] = malloc(sites * dim * sizeof(array[0]));
  if (array[0] == NULL) {
    printf("Failed to allocate memory for a lattice table!\n");
    die(1002);
  }
	for (long i=0; i<sites; i++) {
		array[i] = array[0] + i * dim;
	 }

  return array;
}

// "Reallocate" a lattice table
long **realloc_latticetable(long** arr, int dim, long oldsites, long newsites) {

  long** new_arr = alloc_latticetable(dim, newsites);

  // now just copy to the new array and free the old one
  for (long i=0; i<newsites; i++) {
    for (int dir=0; dir<dim; dir++) {
      if (i >= oldsites) {
        new_arr[i][dir] = 0;
      } else {
        new_arr[i][dir] = arr[i][dir];
      }
    }
  }

  free_latticetable(arr);
  return new_arr;
}


/* Allocate all needed lookup tables needed for layouting
*/
void alloc_lattice_arrays(lattice *l, long sites) {

  l->coords = alloc_latticetable(l->dim, sites);
	l->next = alloc_latticetable(l->dim, sites);
	l->prev = alloc_latticetable(l->dim, sites);

	// allocate parity arrays
	l->parity = malloc(sites * sizeof(*(l->parity)));

	if (!l->rank)
		printf("Allocated memory for lookup tables.\n");
}

// Realloc everything originally allocated in alloc_lattice_arrays
void realloc_lattice_arrays(lattice *l, long oldsites, long newsites) {

  l->coords = realloc_latticetable(l->coords, l->dim, oldsites, newsites);
  l->next = realloc_latticetable(l->next, l->dim, oldsites, newsites);
  l->prev = realloc_latticetable(l->prev, l->dim, oldsites, newsites);

	l->parity = realloc(l->parity, newsites * sizeof(*(l->parity)));
}


/* Free all memory reserved for layouting, lookup tables
* and comlists.
*/
void free_lattice(lattice *l) {

  // comlist
  free_comlist(&l->comlist);

	// lookup tables and parity:
	free(l->parity);
  free_latticetable(l->coords);
	free_latticetable(l->next);
	free_latticetable(l->prev);

	// slicing:
	free(l->sliceL);
	free(l->nslices);
	// then finally the full lattice size, allocated in get_parameters()
	free(l->L);
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

/** @file alloc.c
*
* Routines for memory allocation.

//TODO allocate contagious memory, modify for parallelization
*
*/

#include "su2.h"



/* Allocates memory for a field with a single degree of freedom per site.
*
*/
double *make_singletfield(params p) {

  double *field = malloc(p.vol*sizeof(field));

  return field;
}


/* Allocates memory for a field with dofs degrees of freedom.
* It can be accessed as field[x][dof], where x is the lattice site.
*/
double **make_field(params p, int dofs) {

  double **field = (double **)malloc(p.vol*sizeof(*field));

	for (int i=0; i<p.vol; i++) {
		field[i] = malloc(dofs * sizeof(*(field[i]) ) );
	}

  return field;
}


/* Allocates memory for a gauge field with dofs degrees of freedom.
*	 It can be accessed as field[x][dir][dof], where x is the lattice site
*	 and dir is the direction, pointing towards positive coordinate axes.
*	 the number of directions is taken from p.dim.
*/
double ***make_gaugefield(params p, int dofs) {

  double ***field = (double ***)malloc(p.vol*sizeof(**field));

	for (int x=0; x<p.vol; x++) {
		field[x] = (double **)malloc(p.dim * sizeof(*(field[x]) ) );
		for (int dir=0; dir<p.dim; dir++) {
			field[x][dir] = (double *)malloc(SU2LINK * sizeof(*(field[x][dir])) );
		}
	}

  return field;
}


/* Free the memory allocated by make_singletfield()
*
*/
void free_singletfield(params p, double *field) {
  free(field);
}

/* Free the memory allocated by make_field(p, dofs)
*
*/
void free_field(params p, double **field) {

	for (int i=0; i<p.vol; i++) {
		free(field[i]);
	}

  free(field);
}

/* Free the memory allocated by make_gaugefield(p, dofs)
*
*/
void free_gaugefield(params p, double ***field) {

	for (int i=0; i<p.vol; i++) {
		for (int dir=0; dir<p.dim; dir++) {
			free(field[i][dir]);
		}
		free(field[i]);
	}

  free(field);
}


/** Allocate all the fields needed for a simulation.
 */
void alloc_fields(params p, fields *f) {

	f->su2link = make_gaugefield(p, SU2LINK);

	#ifdef HIGGS
		f->su2doublet = make_field(p, SU2DB);
	#endif
	#ifdef TRIPLET
		f->su2triplet = make_field(p, SU2TRIP);
	#endif

  printf("Successfully allocated memory for fields.\n");
}


/** Free all the fields allocated by alloc_fields().
 */
void free_fields(params p, fields *f) {

	free_gaugefield(p, f->su2link);
	#ifdef HIGGS
		free_field(p, f->su2doublet);
	#endif
	#ifdef TRIPLET
		free_field(p, f->su2triplet);
	#endif

	printf("Successfully freed memory allocated for fields.\n");

}


// Allocate memory for global neighbor pointers
ulong **alloc_neighborList(params p) {

  ulong **ptr = (ulong **)malloc(p.vol*sizeof(*ptr));

	for (ulong i=0; i<p.vol; i++) {
		ptr[i] = malloc(p.dim * sizeof(*(ptr[i]) ) );
	}

  return ptr;
}


// Allocate all needed neighbor lists
void alloc_neighbors(params *p) {

	p->next = alloc_neighborList(*p);
	p->prev = alloc_neighborList(*p);

	// allocate the parity pointers
	p->parity = malloc(p->vol * sizeof(p->parity));

	printf("Successfully allocated memory for neighboring sites.\n");
}

// Free the memory allocated by allocate_neighbors()
void free_neighbors(params *p) {

  for (int i=0; i<p->dim; i++) {
		free(p->next[i]);
		free(p->prev[i]);
	}

  free(p->next);
	free(p->prev);
	printf("Successfully freed memory allocated for neighboring sites.\n");
}

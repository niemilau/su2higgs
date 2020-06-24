
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
}

// Realloc everything originally allocated in alloc_lattice_arrays
void realloc_lattice_arrays(lattice *l, long oldsites, long newsites) {

  l->coords = realloc_latticetable(l->coords, l->dim, oldsites, newsites);
  l->next = realloc_latticetable(l->next, l->dim, oldsites, newsites);
  l->prev = realloc_latticetable(l->prev, l->dim, oldsites, newsites);

	l->parity = realloc(l->parity, newsites * sizeof(*(l->parity)));
}

/* Allocate substructures in a comlist.
* Needs to know how many sendrecv structs to alloc;
* use realloc_comlist() to reallocate these.
* Does not alloc sitelists! */
void alloc_comlist(comlist_struct* comlist, int nodes) {

  comlist->recv_from = malloc(nodes * sizeof(*(comlist->recv_from)));
  comlist->send_to = malloc(nodes * sizeof(*(comlist->send_to)));
  comlist->sends = 0;
  comlist->recvs = 0;
  /*
  for (int k=0; k<nodes; k++) {
    comlist->recv_from[k].sitelist = malloc(sites * sizeof(*(comlist->recv_from[k].sitelist)));
    comlist->send_to[k].sitelist = malloc(sites * sizeof(*(comlist->send_to[k].sitelist)));
  }
  */
}

/* Realloc all substructures in a comlist to only the needed size
* Assumes that site lists in sendrecv_structs have not been allocated
* if there's nothing to send/receive.
* Does not free sendrecv_structs even if they are not needed,
* this is left for free_comlist().
* sendrecv = 0: touch only send_to structs,
* sendrecv = 1: touch only recv_from structs */
void realloc_comlist(comlist_struct* comlist, int sendrecv) {

  int new_size;
  if (sendrecv == SEND) {

    if (comlist->sends > 0) {
      new_size = comlist->sends;
    } else {
      // nothing to send to anyone. Realloc as small as possible, but don't free
      new_size = 1;
    }

    comlist->send_to = realloc(comlist->send_to, new_size * sizeof(*comlist->send_to));
    if (comlist->send_to == NULL) {
      printf("Failed to realloc blocklists (send)!\n");
      die(-9120);
    }
  }

  if (sendrecv == RECV) {
    if (comlist->recvs > 0) {
      new_size = comlist->recvs;
    } else {
      // nothing to receive from anyone. Realloc as small as possible, but don't free
      new_size = 1;
    }

    comlist->recv_from = realloc(comlist->recv_from, new_size * sizeof(*comlist->recv_from));
    if (comlist->recv_from == NULL) {
      printf("Failed to realloc blocklists (recv)!\n");
      die(-9121);
    }
  }

  // then sitelists
  sendrecv_struct* sr;
  int n;

  if (sendrecv == SEND) {
    n = comlist->sends;
  } else if (sendrecv == RECV) {
    n = comlist->recvs;
  }

  for (int k=0; k<n; k++) {

    if (sendrecv == SEND) {
      sr = &comlist->send_to[k];
    } else if (sendrecv == RECV) {
      sr = &comlist->recv_from[k];
    }

    sr->sitelist = realloc(sr->sitelist, sr->sites * sizeof(*sr->sitelist));

    if (sr->sitelist == NULL) {
      printf("Failed to realloc sitelist in blocklist!\n");
      die(-9122);
    }
  }

}

/* Free all memory reserved for layouting, lookup tables
* and comlists.
*/
void free_lattice(lattice *l) {

  // comlist
  free_comlist(&l->comlist);

  #ifdef BLOCKING
    free_comlist(&l->blocklist);
  #endif

	// lookup tables and parity:
	free(l->parity);
  free_latticetable(l->coords);
	free_latticetable(l->next);
	free_latticetable(l->prev);

  #ifdef MEASURE_Z
    free_latticetable(l->site_at_z);
  #endif

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

	for (int k=0; k<comlist->sends; k++) {
		free(comlist->send_to[k].sitelist);
	}
  for (int k=0; k<comlist->recvs; k++) {
		free(comlist->recv_from[k].sitelist);
	}

	free(comlist->recv_from);
	free(comlist->send_to);
}

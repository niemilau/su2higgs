/** @file layout.c
*
* Routines for MPI parallelization on a hypercubic (=rectangle) lattice
*	with periodic boundaries.
*
*/

#include "su2.h"
#include "comms.h"

#ifdef MPI

/* layout()
* Perform all required steps to get sites and halos in place,
* initialize lookup tables and comlists.
*/
void layout(lattice *l, int do_prints, int run_checks) {

	clock_t start, end;
	double time;

	// these are needed for make_slices():
	l->sliceL = malloc(l->dim * sizeof(*(l->sliceL)));
	l->nslices = malloc(l->dim * sizeof(*(l->nslices)));

	l->offset = malloc(l->dim * sizeof(*(l->offset))); // needed for sitemap()

	make_slices(l, do_prints);

	// slicing done, allocate lattice tables
	alloc_lattice_arrays(l, l->sites_total);
	// here both l->halos and l->sites_total still contain self halos
	if (!l->rank && do_prints) {
		printf("Allocated memory for lookup tables.\n");
	}

	long tot_old = l->sites_total;

	/* calculate l->coords and neighbor site tables.
	* This also moves real sites to indices before halo sites and self halos to the end,
	* and updates l->halos, l->sites_total.
	* Realloc the tables here so that all traces of self halos are removed */
	sitemap(l);

	long tot_new = l->sites_total;
	realloc_lattice_arrays(l, tot_old, tot_new);
	// ordering OK and self halos gone.

	if (do_prints)
		printf0(*l, "Each node needs %lu additional halo sites.\n", l->halos);

	// site parity:
	set_parity(l);

	/* Do one last remapping of site indices. We want sites
	* with EVEN parity to come first,
	* then ODD, for both real and halo sites. Ordering is:
	* 1. EVEN real 2. ODD real 3. EVEN halo 4. ODD halo. */
	l->reorder_parity = 1;
	if (l->reorder_parity) {
		long newindex[l->sites_total];
		paritymap(l, newindex);
		remap_lattice_arrays(l, newindex, l->sites_total);
	}

	// --- Site ordering not changed beyond this point ---

	/* find index of the lattice site residing at node coordinates
	* (x,y,z) = (0,0,0). First need the physical coords of this site */
	long xnode[l->dim];
	long x[l->dim];
	indexToCoords(l->dim, l->nslices, l->rank, xnode);
	for (int dir=0; dir<l->dim; dir++) {
		x[dir] = xnode[dir] * l->sliceL[dir];
	}
	l->firstsite = findsite(l, x, 0);


	barrier(l->comm);
	// construct communication tables
	make_comlists(l, &(l->comlist));

	// Run sanity checks on lattice layout and comms?
	barrier(l->comm);
	if (run_checks) {

		start = clock();

		test_coords(l);
		test_neighbors(l);
		test_comms(l);
		test_comms_individual(l);

		end = clock();
		time = ((double) (end - start)) / CLOCKS_PER_SEC;

		if (do_prints) {
			printf0(*l, "All tests OK! Time taken: %lf seconds.\n", time);
		}
	}

	#ifdef MEASURE_Z
		init_z_coord(l);
	#endif

	barrier(l->comm);

}


/* Lay out the lattice on available nodes
* The shapes of the sublattices are not fixed manually:
* the tricks for factorization are heavily inspired by respective routines in Kari's code.
* In fact we should always obtain the same factorization as what Kari has in his susy code.
* The routine checks that the sites can always be factorized evenly among the nodes, i.e, all nodes
* get the same amount of sites. This will be assumed in all other routines in the file.
*/
void make_slices(lattice *l, int do_prints) {

	// first check that the total number of lattice sites is possible
	// to factorize in the given number of nodes (MPI processes)
	int nodes = l->size;
	if (l->vol % nodes != 0) {
		printf0(*l, "Can\'t lay out %lu sites using %d nodes!\n", l->vol, nodes);
		die(100);
	}

	static int n_primes = 10;
  static int prime[] = {2,3,5,7,11,13,17,19,23,29};

	// factorize the number of nodes into primes, i.e. calculate integers a, b, c, d,... >= 0
	// so that number of nodes = (2^a) * (3^b) * (5^c) * (7^d) * ...
	// Since number of sites = Lx * Ly * Lz * ... = integer * nodes,
	// the same primes can be used to factorize the lattice.
  int nfactors[n_primes];
	int i = nodes;

  for (int n=0; n<n_primes; n++) {
    nfactors[n] = 0;
    while (i % prime[n] == 0) {
      i = i/prime[n];
      nfactors[n]++;
    }
  }
	if(i != 1) {
		// this can happen if we don't have listed enough primes
    printf0(*l, "Can\'t factorize %d nodes with primes up to %d! List more primes in layout.c.\n",
	      l->vol, prime[n_primes-1]);

    die(101);
  }

	// Now factorize the lattice side lengths in terms of these primes.
	// We start with full lengths in l->L[dir] and slice the largest length to get a smaller slice,
	// then just repeat. So the steps are: 1) find currently largest slice 2) slice it 3) repeat.
	long slice[l->dim]; // how many lattice sites in one slice in direction dir
	int nslices[l->dim]; // how many hypercubes can we fit in one direction
	int dir;
	// initialize these
	for (dir=0; dir<l->dim; dir++) {
		slice[dir] = l->L[dir];
		nslices[dir] = 1;
	}

	int firstdim = 0;
	/*
	#ifdef WALL
		// for wall profiling, ONLY slice the last direction! (OUTDATED, works with general slices now)
		if (l->L[l->dim-1] % nodes != 0) {
			printf0(*l, "Wall routines failed: cannot split z-direction into %d pieces! in layout.c\n", nodes);
			die(419);
		}
		firstdim = l->dim-1;
	#endif
	*/

	// start with the largest prime!
	for (int n=n_primes-1; n>=0; n--) {
		// whatever we do for one prime, repeat as many times as it appears in node factorization
		for (int f=0; f<nfactors[n]; f++) {

			// find largest slice direction that we can still chop with this prime.
			int largest=1;
			for (dir=firstdim; dir<l->dim; dir++) {
				if ( (slice[dir] > largest) && (slice[dir] % prime[n] == 0) ) {
					largest = slice[dir];
				}
			}

			// now chop the longest direction
			for (dir=firstdim; dir<l->dim; dir++) {
				if (slice[dir] == largest) {
					slice[dir] = slice[dir] / prime[n];
					nslices[dir] *= prime[n];
					break;
				}
			}
			// this then repeats as long as it is possible to still chop some direction into smaller slices.
		}
	}

	/* We have now sliced the full lattice into smaller hypercubes with side lengths given by sliceL.
	*	p.rank is to be understood as the index of the sublattice where the node operates (comparable to lattice site index).
	* This is implemented in sitemap() below.
	*/


	// print the obtained processor layout and node sizes
	if (!l->rank && do_prints) {
		printf("Processor layout: ");
		for (dir=0; dir<l->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%d", nslices[dir]);
		}
		printf("\nSites on each node: ");
		for (dir=0; dir<l->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%ld", slice[dir]);
		}
	}


	// store these.
	long sites = 1;
	long halosites = 1; // 2 extra sites in each direction
	for (dir=0; dir<l->dim; dir++) {
		l->sliceL[dir] = slice[dir];
		l->nslices[dir] = nslices[dir];
		sites *= slice[dir];
		halosites *= (slice[dir]+2);
	}
	halosites = halosites - sites;
	if (do_prints)
		printf0(*l, " = %lu.\n",sites);

	l->sites = sites;
	// these will be modified later when we remove unneeded halos:
	l->halos = halosites;
	l->sites_total = sites + halosites;
}


/* Construct lookup tables for adjacent sites, including halos.
* Coordinates on the full lattice are stored in l->coords.
* This is first done using indices in the full haloed node,
* which are then remapped so that "self halos" are mapped away
* and real sites come before halos in indexing.
*/
void sitemap(lattice *l) {

	long maxindex = l->sites + l->halos;
	long i;
	int dir;
	// work arrays
	long xnode[l->dim]; // coordinates of the MPI nodes
	long x[l->dim]; // (x,y,z,...) coords on the slice
	char ishalo[maxindex]; // 1 if site is in the halo, 0 otherwise
	int whichnode[maxindex]; // rank of the node where site i resides in

	// where is the node located?
	indexToCoords(l->dim, l->nslices, l->rank, xnode);

	// store coordinate offset wrt. to the full lattice in each direction:
	for (dir=0; dir<l->dim; dir++) {
		l->offset[dir] = xnode[dir] * l->sliceL[dir];
	}

	/* now build a halo of thickness 1 around the node. Sites in the halo actually live
	* in other nodes and the first step is to separate these from "real" sites */
	int haloL[l->dim];
	for (dir=0; dir<l->dim; dir++) {
		haloL[dir] = l->sliceL[dir] + 2;
	}

	/* Find neighbors and l->coords in terms of indices in the haloed node.
	* For sites in the positive (negative) halo, we only need to know their neighbors
	* in the negative (positive) directions. Other neighbor indices are set to -1.
	* If we ever try to access these neighbors the program should crash, but this is
	* better than getting wrong results.
	*/

	// count how many "fake" halo sites we have. these are sites
	// in the halo that actually live in the same node because of periodicity
	long selfhalos = 0;

	for (i=0; i<maxindex; i++) {

		indexToCoords(l->dim, haloL, i, x);

		ishalo[i] = 0;

		for (dir=0; dir<l->dim; dir++) {
			int prev_done = 0, next_done = 0;

			// Assuming we did not cross the lattice boundary:
			l->coords[i][dir] = x[dir] + xnode[dir] * l->sliceL[dir] - 1; // halos OK too

			// figure out if we are in the halo, and check if we actually crossed the boundary
			if (x[dir] == l->sliceL[dir] + 1) {
				ishalo[i] = 1;
				// we are in the next node in positive direction.
				l->next[i][dir] = -1;
				next_done = 1;

				//Did we go over the lattice boundary?
				if (xnode[dir]+1 >= l->nslices[dir]) {
					// went over, force periodicity
					l->coords[i][dir] = 0;
				}
			} else if (x[dir] == 0) {
				ishalo[i] = 1;
				// we are in the previous node in this direction
				l->prev[i][dir] = -1;
				prev_done = 1;

				if (xnode[dir]-1 < 0) {
					// went over lattice boundary
					l->coords[i][dir] = l->L[dir]-1;
				}
			}

			// find neighbor sites
			if (!next_done) {
				x[dir]++;
				l->next[i][dir] = coordsToIndex(l->dim, haloL, x);
				x[dir]--;
			}
			if (!prev_done) {
				x[dir]--;
				l->prev[i][dir] = coordsToIndex(l->dim, haloL, x);
				x[dir]++;
			}

		} // end dir loop

		// which node does the site live in?
		whichnode[i] = coordsToRank(l, l->coords[i]);
		if (ishalo[i] && whichnode[i] == l->rank) {
			selfhalos++;
		} else if (!ishalo[i] && whichnode[i] != l->rank) {
			printf("Error in layout.c! mapping failed in node %d \n", l->rank);
			die(552);
		}

	} // end i loop, physical coords and neighbors done

	// get rid of the self halos
	l->halos -= selfhalos;
	// actual number of lattice sites that we need to care about:
	l->sites_total = l->sites + l->halos;

	/* Rearrange indices: non-halo sites should come before halos.
	* In case of self halos, move those to the end so that we can forget about them. */
	long newindex[maxindex];
	long n=0, k=0, j=0;

	for (i=0; i<maxindex; i++) {
		if (ishalo[i] && whichnode[i] == l->rank) {
			// self halo site, so map it to something we don't care about anymore
			newindex[i] = l->sites_total + n;
			n++;
		}
		else if (!ishalo[i]) {
			// non-halo
			newindex[i] = j;
			j++;
		} else {
			// normal halo site
			newindex[i] = l->sites + k;
			k++;
		}
	}

	// remap all
	remap_lattice_arrays(l, newindex, maxindex);

	/* Almost done, but in p.next and p.prev some sites can still point to a
	* neighboring self halo site index instead of the real site. To change these,
	* search for matching coordinates in l->coords
	*/
	for (i=0; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			// positive dirs
			long next = l->next[i][dir];
			if (next >= l->sites_total) {
				// points to a self halo site, find the matching real site
				for (int d=0; d<l->dim; d++) x[d] = l->coords[next][d];
				// no need to search in halos here
				l->next[i][dir] = findsite(l, x, 0);
			}
			// negative dirs
			long prev = l->prev[i][dir];
			if (prev >= l->sites_total) {
				for (int d=0; d<l->dim; d++) x[d] = l->coords[prev][d];

				l->prev[i][dir] = findsite(l, x, 0);
			}
		}
	}

	// tables done and self halos moved to indices > l->sites_total. realloc in layout()
}


/* Set checkerboard parity for all sites (excl. self halos,
* which are assumed to be mapped out already).
* Parity of a site is defined to be EVEN if the sum of
* its physical coordinates x + y + z + ... is an even number,
* ODD otherwise. These are chars defined in a header file (where?)
*/
void set_parity(lattice *l) {

	long tot;
	l->evensites = 0; l->oddsites = 0;
	l->evenhalos = 0; l->oddhalos = 0;

	// real sites
	for (long i=0; i<l->sites; i++) {
		tot = 0;
		for (int dir=0; dir<l->dim; dir++) {
			tot += l->coords[i][dir];
		}
		if (tot % 2 == 0) {
			l->parity[i] = EVEN;
			l->evensites++;
		} else {
			l->parity[i] = ODD;
			l->oddsites++;
		}
	}
	// halo sites
	for (long i=l->sites; i<l->sites_total; i++) {
		tot = 0;
		for (int dir=0; dir<l->dim; dir++) {
			tot += l->coords[i][dir];
		}
		if (tot % 2 == 0) {
			l->parity[i] = EVEN;
			l->evenhalos++;
		} else {
			l->parity[i] = ODD;
			l->oddhalos++;
		}
	}

}


/* Quick routine for finding site index from given
* physical coordinates. Assumes ordered l->coords.
*/
long findsite(lattice const* l, long* x, int include_halos) {

	long max;
	if (include_halos) {
		max = l->sites_total;
	} else {
		max = l->sites;
	}

	for (long i=0; i<max; i++) {
		int match = 1;
		for (int dir=0; dir<l->dim; dir++) {
			if (x[dir] != l->coords[i][dir]) {
				match = 0;
				break;
			}
		}

		if (match) {
			return i;
		}
	}

	// no match after searching through all sites?
	printf("Node %d: Failed to find matching site in l->coords!\n", l->rank);
	return -1;
}


/* Test that our mapping from site index to physical coordinates makes sense.
* This checks that basic indexing is OK, but does not check that neighboring sites
* have adjacent indices (which they should, apart from sites at hypercube sides).
*/
void test_coords(lattice const* l) {

	long* xnode = malloc(l->dim * sizeof(*xnode));

	indexToCoords(l->dim, l->nslices, l->rank, xnode);

	int dir;
	long i;


	// first site on the node?
	// we remapped sites by parity, so first site is either i=0 or i=p.evensites
	for (dir=0; dir<l->dim; dir++) {
		if (l->coords[0][dir] != xnode[dir] * l->sliceL[dir] && l->coords[l->evensites][dir] != xnode[dir] * l->sliceL[dir]) {
			printf("Node %d: Error in test_coords! First site not where it should be \n", l->rank);
			die(-111);
		}
	}
	// last real site on the node?
	for (dir=0; dir<l->dim; dir++) {
		if (l->coords[l->sites-1][dir] != (xnode[dir] + 1)* l->sliceL[dir] - 1 && l->coords[l->evensites-1][dir] != (xnode[dir] + 1)* l->sliceL[dir] - 1) {
			printf("Node %d: Error in test_coords! Last site not where it should be \n", l->rank);
			die(-112);
		}
	}


	// non-halo sites: coordinate x[j] should be in range
	// xnode[j] * l->sliceL[j] <= x[j] <= xnode[j] * l->sliceL[j] + l->sliceL[j] - 1
	for (i=0; i<l->sites; i++) {
		for (dir=0; dir<l->dim; dir++) {
			if (l->coords[i][dir] > (xnode[dir] + 1) * l->sliceL[dir] - 1 || xnode[dir] * l->sliceL[dir] > l->coords[i][dir]) {
				printf("Node %d: Error in test_coords! Real site %ld not indexed properly \n", l->rank, i);
				die(-113);
			}
		}
	}

	// halo sites: should not find a halo site whose all coordinates x[j]
	// match a real site on our node
	for (i=l->sites; i<l->sites_total; i++) {
		int ok = 0;
		for (dir=0; dir<l->dim; dir++) {
			if (l->coords[i][dir] > (xnode[dir] + 1) * l->sliceL[dir] - 1 || xnode[dir] * l->sliceL[dir] > l->coords[i][dir]) {
				ok = 1;
			}
		}
		if (!ok) {
			printf("Node %d: Error in test_coords! Halo site %ld not indexed properly \n", l->rank, i);
			die(-114);
		}
	}

	free(xnode);
}

/* Perform strong checks on l.next and l.prev, and l.parity.
* Uses l->coords to check that the physical coordinates of the neighbor sites
* are what they should for each site and direction, including halos.
*/
void test_neighbors(lattice const* l) {

	int skip_parity = 0;
	for (int dir=0; dir<l->dim; dir++) {
		if (l->L[dir] % 2 != 0) {
			// side lengths not even numbers, skip parity checks...
			skip_parity = 1;
		}
	}
	// first site in master node should have even parity
	if (!skip_parity && l->rank == 0) {
		if (l->parity[0] != EVEN) {
			printf("Node 0: Error in test_neighbors! First site parity not is not EVEN \n");
			// don't die
		}
	}

	long x[l->dim];

	for (long i=0; i<l->sites_total; i++) {

		// first check neighbors in the positive directions
		for (int dir=0; dir<l->dim; dir++) {
			long next = l->next[i][dir];
			if (next == -1) {
				if (i < l->sites) {
					// no neighbor assigned to a real lattice site, error!
					printf("Node %d: Error in test_neighbors, site %ld! Next site in direction %d not assigned \n", l->rank, i, dir);
					die(-116);
				}
			} else {
				// check that coordinates of the next site make sense.
				// What we think they should be:
				for (int d=0; d<l->dim; d++) {
					x[d] = l->coords[i][d];
				}
				x[dir]++;
				// Periodicity?
				if (x[dir] >= l->L[dir]) {
					x[dir] = 0;
				}

				// what they are according to l->coords
				for (int d=0; d<l->dim; d++) {
					if (l->coords[next][d] != x[d]) {
						printf("Node %d: Error in test_neighbors, site %ld! Coordinates of next site in direction %d do not match! \n", l->rank, i, dir);
						die(-117);
					}
				}

				// parity check
				if (!skip_parity) {
					if (l->parity[i] == l->parity[next]) {
						printf("Node %d: Error in test_neighbors, site %ld! Next site in direction %d has same parity! \n", l->rank, i, dir);
					}
				}
			}
		}

		// then same checks for negative directions
		for (int dir=0; dir<l->dim; dir++) {
			long prev = l->prev[i][dir];
			if (prev == -1) {
				if (i < l->sites) {
					// no neighbor assigned to a real lattice site, error!
					printf("Node %d: Error in test_neighbors, site %ld! Previous site in direction %d not assigned \n", l->rank, i, dir);
					die(-118);
				}
			} else {
				// check that coordinates of the next site make sense.
				// What we think they should be:
				for (int d=0; d<l->dim; d++) {
					x[d] = l->coords[i][d];
				}
				x[dir]--;
				// Periodicity?
				if (x[dir] < 0) {
					x[dir] = l->L[dir] - 1;
				}

				// what they are according to l->coords
				for (int d=0; d<l->dim; d++) {
					if (l->coords[prev][d] != x[d]) {
						printf("Node %d: Error in test_neighbors, site %ld! Coordinates of previous site in direction %d do not match! \n", l->rank, i, dir);
						die(-119);
					}
				}

				// parity check
				if (!skip_parity) {
					if (l->parity[i] == l->parity[prev]) {
						printf("Node %d: Error in test_neighbors, site %ld! Previous site in direction %d has same parity! \n", l->rank, i, dir);
					}
				}

			}
		}
	}

}

/** Formatted printing function, root node only.
 *
 * Drop-in replacement for `fprintf(stderr,...)` that only
 * prints stuff if called by the master node, rank 0.
 *
 * Taken from David.
 */
void printf0(lattice l, char *msg, ...) {
  va_list fmtargs;
  va_start(fmtargs, msg);

  if(!l.rank)
    vprintf(msg, fmtargs);
}


/** Call `MPI_Finalize` and quit.
 *
 * The return value is 'howbad'.  Note that all nodes have to call
 * this, or the code may hang.
 *
 * Taken from David
 */
void die(int howbad) {
  MPI_Finalize();
  exit(howbad);
}

#else // No MPI, dummy routines for layouting

void layout(lattice *l, int do_prints, int run_checks) {

	l->rank = 0;
	l->size = 1;
	l->sites = l->vol;
	l->sites_total = l->vol;
	l->halos = 0;

	l->sliceL = malloc(l->dim * sizeof(*l->sliceL));
	l->nslices = malloc(l->dim * sizeof(*l->nslices));

	l->offset = malloc(l->dim * sizeof(*(l->offset)));

	for (int dir=0; dir<l->dim; dir++) {
		l->sliceL[dir] = l->L[dir];
		l->nslices[dir] = 1;
		l->offset[dir] = 0;
	}

	alloc_lattice_arrays(l, l->sites_total);

	// construct lookup tables for site neighbors and parity
	sitemap(l);
	set_parity(l);

	// reorder by parity
	l->reorder_parity = 1;
	if (l->reorder_parity) {
		long newindex[l->sites_total];
		paritymap(l, newindex);
		remap_lattice_arrays(l, newindex, l->sites_total);
	}

	alloc_comlist(&l->comlist, 1); // empty comlist

	if (do_prints) {
		printf("Site lookup tables constructed succesfully.\n");
	}

	#ifdef MEASURE_Z
		init_z_coord(l);
	#endif
}

void sitemap(lattice* l) {

	l->oddsites = 0; l->oddhalos = 0;
	l->evensites = 0; l->evenhalos = 0;

	for (long i=0; i<l->sites_total; i++) {

		// get the physical coordinates of site i and store in p.coords[i]
		indexToCoords(l->dim, l->L, i, l->coords[i]);

		for (int dir=0; dir < l->dim; dir++) {
			l->coords[i][dir]++;
			l->next[i][dir] = coordsToIndex(l->dim, l->L, l->coords[i]);
			l->coords[i][dir] -= 2;
			l->prev[i][dir] = coordsToIndex(l->dim, l->L, l->coords[i]);
			// return to the original value
			l->coords[i][dir]++;
		}

	}

}


void set_parity(lattice *l) {

	l->oddsites = 0; l->oddhalos = 0;
	l->evensites = 0; l->evenhalos = 0;
	long* x = malloc(l->dim * sizeof(*x));

	for (long i=0; i<l->sites_total; i++) {

		// get the physical coordinates of site i and store in x
		indexToCoords(l->dim, l->L, i, x);
		// calculate and store the parity of site i
		long coord = 0;
		for (int j=0; j<l->dim; j++) {
			coord += x[j];
		}

		if (coord % 2 == 0) {
			l->parity[i] = EVEN;
			l->evensites++;
		} else {
			l->parity[i] = ODD;
			l->oddsites++;
		}
	}

	free(x);
}


void printf0(lattice l, char *msg, ...) {
  va_list fmtargs;
  va_start(fmtargs, msg);

  vprintf(msg, fmtargs);
}

void die(int howbad) {
  exit(howbad);
}

#endif


/******************************************
*	Mutual routines for both serial and MPI *
******************************************/

/* Order sites according to their parity.
* Assumes that self halos have been removed and sites ordered so that
* real sites come before halos. Does not reorder EVEN (ODD) sites among themselves.
* Note that the routine overrides newindex table.
*/
void paritymap(lattice* l, long* newindex) {

	long even=0, odd=0, evenhalo=0, oddhalo=0;

	// order real sites
	for (long i=0; i<l->sites; i++) {
		if (l->parity[i] == EVEN) {
			newindex[i] = even;
			even++;
		} else {
			newindex[i] = l->evensites + odd;
			odd++;
		}
	}
	// order halo sites
	for (long i=l->sites; i<l->sites_total; i++) {
		if (l->parity[i] == EVEN) {
			newindex[i] = l->sites + evenhalo;
			evenhalo++;
		} else {
			newindex[i] = l->sites + l->evenhalos + oddhalo;
			oddhalo++;
		}
	}

}

/* Reorder a given lattice table using the mapping given in newindex.
* The table should be l.dim * maxindex sized.
*/
void remap_latticetable(lattice* l, long** arr, long* newindex, long maxindex) {
	long** temp = alloc_latticetable(l->dim, maxindex);

	int dir;
	for (long i=0; i<maxindex; i++) {
		for (dir=0; dir<l->dim; dir++) {
			temp[i][dir] = arr[i][dir];
		}
	}

	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		for (dir=0; dir<l->dim; dir++) {
			arr[new][dir] = temp[i][dir];
		}
	}

	free_latticetable(temp);
}

// Same as remap_latticetable() but specifically for neighbor lookup tables
void remap_neighbor_table(lattice* l, long** arr, long* newindex, long maxindex) {
	long** temp = alloc_latticetable(l->dim, maxindex);

	int dir;
	for (long i=0; i<maxindex; i++) {
		for (dir=0; dir<l->dim; dir++) {
			temp[i][dir] = arr[i][dir];
		}
	}
	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		for (dir=0; dir<l->dim; dir++) {
			// next[i][dir] can point to index -1 if the neighbor is outside the haloed node
			int next = temp[i][dir];
			if (next == -1) {
				arr[new][dir] = -1;
			} else {
				arr[new][dir] = newindex[next];
			}
		}
	}

	free_latticetable(temp);
}

/* Remap all lattice tables using mapping given in newindex.
*/
void remap_lattice_arrays(lattice* l, long* newindex, long maxindex) {

	// backup parity
	char par[maxindex];
	for (long i=0; i<maxindex; i++) {
		par[i] = l->parity[i];
	}

	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		l->parity[new] = par[i];
	}

	remap_neighbor_table(l, l->next, newindex, maxindex);
	remap_neighbor_table(l, l->prev, newindex, maxindex);

	remap_latticetable(l, l->coords, newindex, maxindex);
}


/***********************************************************************
* 	Routines for converting from site index to coordinates vice versa  *
************************************************************************/

// Calculate product L[0]*L[1]*...*L[max-1]. In practice L is the lattice or slice length.
inline long Lprod(int* L, int max) {
	long res = 1;
	for (int i=0; i<max; i++) {
		res *= L[i];
	}

	return res;
}

/* From the site index, get the physical coordinates (x, y, z,...) on a dim-dimensional
* lattice with side lengths defined as in L, and store in x.
* Can be used for the full lattice, or a single slice!
* see Mathematica notebook coordinates.nb for analytical relations.
* This only works in a "natural" ordering of site indices
*/
void indexToCoords(short dim, int* L, long i, long* x) {

	x[0] = i % L[0];
	// we want the floor() of integer division here. This is automatic in C, but I'm being explicit here.
	// possible problem here: floor(0.9999) returns 0, when the analytical value should be 1?
	x[dim-1] = (int)floor(i / (Lprod(L, dim-1)) );

	for (int dir=dim-2; dir>0; dir--) {
		int a=0;
			for (int k=dir+1; k<dim; k++) {
				a += x[k] * Lprod(L, k);
  		}
			// a is at most i here so can't get a negative x.
		x[dir] = (int)floor((i - a) / Lprod(L, dir));
	}

}


// Same as indexToCoords(), but convert from the (x, y, z, ...) to the "natural" site index.
long coordsToIndex(short dim, int* L, long* x) {

		long res = 0;
		for (int dir=0; dir<dim; dir++) {
			res += ( (x[dir] + L[dir]) % L[dir] ) * Lprod(L, dir);
		}
		return res;
}

/* Convert physical (x, y, z, ...) coordinates on the full lattice
* to MPI node index (which is the same as MPI rank in this layout).
*/
int coordsToRank(lattice const* l, long* coords) {
	// first find the (x, y, z, ...) coordinates of the node
	long* x = malloc( (l->dim) * sizeof(x));
	for (int dir=0; dir<l->dim; dir++) {
		x[dir] = floor(coords[dir] / l->sliceL[dir]);
	}
	// then convert the coordinates to node index
	int i = coordsToIndex(l->dim, l->nslices, x);
	free(x);

	if (i >= l->size) {
		printf("WARNING (Node %d): error in layout.c, coordsToRank()\n", l->rank);
	}
	return i;
}


/* Traverse the (sub-)lattice in a "natural" order,
* i.e. take coords (x,y,z) in my node. Then this moves as
* (0,0,0) -> (1,0,0) -> (2,0,0) -> ... (0,1,0) -> (1,1,0) -> ... (0,0,1) -> ...
* Not used for anything atm.
*/
void traverse_natural(lattice* l) {

	long steps[l->dim];
	for (int dir=0; dir<l->dim; dir++) {
		steps[dir] = 0;
	}

	// start from (x,y,z) = (0,0,0) in my node
	long start = l->firstsite;
	long site = start;

	for (long j=0; j<l->sites; j++) {

		// Debug
		/*
		printf("( ");
		long x[l->dim];
		for (int dir=0; dir<l->dim; dir++) {
			printf("%ld, ", l->coords[site][dir]);
		}
		printf(" )\n");
		*/

		// first attempt moving one step in the x0 direction
		if (steps[0] < l->sliceL[0] - 1) {
			site = l->next[site][0];
			steps[0]++;

		} else {
			int movenext = 0;
			// could not move, so attempt x1, x2 etc (in that order)
			for (int dir=0; dir<l->dim; dir++) {

				if (steps[dir] == l->sliceL[dir] - 1) {
					movenext = 1; // in next iteration, move in the next dir

					// backtrack the start point to x_dir = 0
					// note that start point for x_0 is never changed
					if (dir != 0) {
						for (long k=0; k<steps[dir]; k++) {
							start = l->prev[start][dir];
						}
					}

					steps[dir] = 0;
					continue;
				} else if (movenext) {
					start = l->next[start][dir]; // move once in this new direction
					site = start;
					steps[dir]++;
					movenext = 0;
					break;
				}
			} // end dir

		} // end else

	} // end j

}

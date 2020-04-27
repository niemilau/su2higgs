/** @file layout.c
*
* Routines for MPI parallelization on a hypercubic (=rectangle) lattice
*	with periodic boundaries.
*
*/

#include "su2.h"
#include "comms.h"
#ifdef WALL
		#include "wallprofile.h"
#endif

#ifdef MPI

/* layout()
* Perform all required steps to get sites and halos in place,
* initialize lookup tables and comlists.
*/
void layout(params *p, comlist_struct *comlist) {

	clock_t start, end;
	double time;

	// these are needed for make_slices():
	p->sliceL = malloc(p->dim * sizeof(p->sliceL));
	p->nslices = malloc(p->dim * sizeof(p->nslices));

	make_slices(p);

	// slicing done, allocate lattice tables
	alloc_lattice_arrays(p, p->sites_total);
	// here both p->halos and p->sites_total still contain self halos

	long tot_old = p->sites_total;

	/* calculate p->coords and neighbor site tables.
	* This also moves real sites to indices before halo sites and self halos to the end,
	* and updates p->halos, p->sites_total.
	* Realloc the tables here so that all traces of self halos are removed */
	sitemap(p);

	long tot_new = p->sites_total;
	realloc_lattice_arrays(p, tot_old, tot_new);
	// ordering OK and self halos gone.

	printf0(*p, "Each node needs %lu additional halo sites.\n", p->halos);

	// site parity:
	set_parity(p);

	/* Do one last remapping of site indices. We want sites
	* with EVEN parity to come first,
	* then ODD, for both real and halo sites. Ordering is:
	* 1. EVEN real 2. ODD real 3. EVEN halo 4. ODD halo. */
	p->reorder_parity = 1;
	if (p->reorder_parity) {
		long newindex[p->sites_total];
		paritymap(p, newindex);
		remap_lattice_arrays(p, newindex, p->sites_total);
	}

	// --- Site ordering not changed beyond this point ---

	#ifdef MEASURE_Z
		//


	#endif

	MPI_Barrier(MPI_COMM_WORLD);

	// construct communication tables
	make_comlists(p, comlist);

	// Run sanity checks on lattice layout and comms?
	MPI_Barrier(MPI_COMM_WORLD);
	if (p->run_checks) {

		start = clock();

		test_coords(*p);
		test_neighbors(*p);
		test_comms_individual(*p, *comlist);
		test_comms(*p, *comlist);

		end = clock();
		time = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf0(*p, "All tests OK! Time taken: %lf seconds.\n", time);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	#ifdef WALL
		printf0(*p, "Initializing wall profile structures\n");
		// initialize stuff defined in wallprofile.h

		sites_per_z = 1;
		if (p->dim == 1) {
			sites_per_z = p->sliceL[0];
		}
		else {
			for (int d=0; d<p->dim-1; d++) {
				sites_per_z *= p->sliceL[d];
			}
		}

		// allocate table so that wallcoord[z][n] gives the site index i
		// of site with z-coordinate z, and n runs over all sites with that z index.
		wallcoord = alloc_latticetable(sites_per_z, p->sliceL[p->dim-1]);

		long* xnode = malloc(p->dim * sizeof(xnode)); // coordinates of the MPI nodes
		indexToCoords(p->dim, p->nslices, p->rank, xnode);

		offset_z = xnode[p->dim-1] * p->sliceL[p->dim-1];

		for (long nz=0; nz<p->sliceL[p->dim-1]; nz++) {
			long tot = 0;
			for (long i=0; i<p->sites; i++) {
				if (nz + offset_z == p->coords[i][p->dim-1]) {
					wallcoord[nz][tot] = i;
					tot++;
				}
			}
			if (tot != sites_per_z) {
				printf("Error in wall profile routines in layout.c!\n");
				die(420);
			}
		}

		free(xnode);

	#endif // end WALL

}


/* Lay out the lattice on available nodes
* The shapes of the sublattices are not fixed manually:
* the tricks for factorization are heavily inspired by respective routines in Kari's code.
* In fact we should always obtain the same factorization as what Kari has in his susy code.
* The routine checks that the sites can always be factorized evenly among the nodes, i.e, all nodes
* get the same amount of sites. This will be assumed in all other routines in the file.
*/
void make_slices(params *p) {

	// first check that the total number of lattice sites is possible
	// to factorize in the given number of nodes (MPI processes)
	int nodes = p->size;
	if (p->vol % nodes != 0) {
		printf0(*p, "Can\'t lay out %lu sites using %d nodes!\n", p->vol, nodes);
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
    printf0(*p, "Can\'t factorize %d nodes with primes up to %d! List more primes in layout.c.\n",
	      p->vol, prime[n_primes-1]);

    die(101);
  }

	// Now factorize the lattice side lengths in terms of these primes.
	// We start with full lengths in p->L[dir] and slice the largest length to get a smaller slice,
	// then just repeat. So the steps are: 1) find currently largest slice 2) slice it 3) repeat.
	long slice[p->dim]; // how many lattice sites in one slice in direction dir
	int nslices[p->dim]; // how many hypercubes can we fit in one direction
	int dir;
	// initialize these
	for (dir=0; dir<p->dim; dir++) {
		slice[dir] = p->L[dir];
		nslices[dir] = 1;
	}

	int firstdim;
	#ifdef WALL
		// for wall profiling, ONLY slice the last direction!
		if (p->L[p->dim-1] % nodes != 0) {
			printf0(*p, "Wall routines failed: cannot split z-direction into %d pieces! in layout.c\n", nodes);
			die(419);
		}
		firstdim = p->dim-1;
	#else
		firstdim = 0;
	#endif

	// start with the largest prime!
	for (int n=n_primes-1; n>=0; n--) {
		// whatever we do for one prime, repeat as many times as it appears in node factorization
		for (int f=0; f<nfactors[n]; f++) {

			// find largest slice direction that we can still chop with this prime.
			int largest=1;
			for (dir=firstdim; dir<p->dim; dir++) {
				if ( (slice[dir] > largest) && (slice[dir] % prime[n] == 0) ) {
					largest = slice[dir];
				}
			}

			// now chop the longest direction
			for (dir=firstdim; dir<p->dim; dir++) {
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
	if (!p->rank) {
		printf("Processor layout: ");
		for (dir=0; dir<p->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%d", nslices[dir]);
		}
		printf("\nSites on each node: ");
		for (dir=0; dir<p->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%ld", slice[dir]);
		}
	}


	// store these.
	long sites = 1;
	long halosites = 1; // 2 extra sites in each direction
	for (dir=0; dir<p->dim; dir++) {
		p->sliceL[dir] = slice[dir];
		p->nslices[dir] = nslices[dir];
		sites *= slice[dir];
		halosites *= (slice[dir]+2);
	}
	halosites = halosites - sites;
	printf0(*p, " = %lu.\n",sites);

	p->sites = sites;
	// these will be modified later when we remove unneeded halos:
	p->halos = halosites;
	p->sites_total = sites + halosites;
}


/* Construct lookup tables for adjacent sites, including halos.
* Coordinates on the full lattice are stored in p->coords.
* This is first done using indices in the full haloed node,
* which are then remapped so that "self halos" are mapped away
* and real sites come before halos in indexing.
*/
void sitemap(params *p) {

	long maxindex = p->sites + p->halos;
	long i;
	int dir;
	// work arrays
	long xnode[p->dim]; // coordinates of the MPI nodes
	long x[p->dim]; // (x,y,z,...) coords on the slice
	char ishalo[maxindex]; // 1 if site is in the halo, 0 otherwise
	int whichnode[maxindex]; // rank of the node where site i resides in

	// where is the node located?
	indexToCoords(p->dim, p->nslices, p->rank, xnode);

	/* now build a halo of thickness 1 around the node. Sites in the halo actually live
	* in other nodes and the first step is to separate these from "real" sites */
	int haloL[p->dim];
	for (dir=0; dir<p->dim; dir++) {
		haloL[dir] = p->sliceL[dir] + 2;
	}

	/* Find neighbors and p->coords in terms of indices in the haloed node.
	* For sites in the positive (negative) halo, we only need to know their neighbors
	* in the negative (positive) directions. Other neighbor indices are set to -1.
	* If we ever try to access these neighbors the program should crash, but this is
	* better than getting wrong results.
	*/

	// count how many "fake" halo sites we have. these are sites
	// in the halo that actually live in the same node because of periodicity
	long selfhalos = 0;

	for (i=0; i<maxindex; i++) {

		indexToCoords(p->dim, haloL, i, x);

		ishalo[i] = 0;

		for (dir=0; dir<p->dim; dir++) {
			int prev_done = 0, next_done = 0;

			// Assuming we did not cross the lattice boundary:
			p->coords[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1; // halos OK too

			// figure out if we are in the halo, and check if we actually crossed the boundary
			if (x[dir] == p->sliceL[dir] + 1) {
				ishalo[i] = 1;
				// we are in the next node in positive direction.
				p->next[i][dir] = -1;
				next_done = 1;

				//Did we go over the lattice boundary?
				if (xnode[dir]+1 >= p->nslices[dir]) {
					// went over, force periodicity
					p->coords[i][dir] = 0;
				}
			} else if (x[dir] == 0) {
				ishalo[i] = 1;
				// we are in the previous node in this direction
				p->prev[i][dir] = -1;
				prev_done = 1;

				if (xnode[dir]-1 < 0) {
					// went over lattice boundary
					p->coords[i][dir] = p->L[dir]-1;
				}
			}

			// find neighbor sites
			if (!next_done) {
				x[dir]++;
				p->next[i][dir] = coordsToIndex(p->dim, haloL, x);
				x[dir]--;
			}
			if (!prev_done) {
				x[dir]--;
				p->prev[i][dir] = coordsToIndex(p->dim, haloL, x);
				x[dir]++;
			}

		} // end dir loop

		// which node does the site live in?
		whichnode[i] = coordsToRank(*p, p->coords[i]);
		if (ishalo[i] && whichnode[i] == p->rank) {
			selfhalos++;
		} else if (!ishalo[i] && whichnode[i] != p->rank) {
			printf("Error in layout.c! mapping failed in node %d \n", p->rank);
		}

	} // end i loop, physical coords and neighbors done

	// get rid of the self halos
	p->halos -= selfhalos;
	// actual number of lattice sites that we need to care about:
	p->sites_total = p->sites + p->halos;

	/* Rearrange indices: non-halo sites should come before halos.
	* In case of self halos, move those to the end so that we can forget about them. */
	long newindex[maxindex];
	long n=0, k=0, j=0;

	for (i=0; i<maxindex; i++) {
		if (ishalo[i] && whichnode[i] == p->rank) {
			// self halo site, so map it to something we don't care about anymore
			newindex[i] = p->sites_total + n;
			n++;
		}
		else if (!ishalo[i]) {
			// non-halo
			newindex[i] = j;
			j++;
		} else {
			// normal halo site
			newindex[i] = p->sites + k;
			k++;
		}
	}

	// remap all
	remap_lattice_arrays(p, newindex, maxindex);

	/* Almost done, but in p.next and p.prev some sites can still point to a
	* neighboring self halo site index instead of the real site. To change these,
	* search for matching coordinates in p->coords
	*/
	for (i=0; i<p->sites_total; i++) {
		for (dir=0; dir<p->dim; dir++) {
			// positive dirs
			long next = p->next[i][dir];
			if (next >= p->sites_total) {
				// points to a self halo site, find the matching real site
				for (int d=0; d<p->dim; d++) x[d] = p->coords[next][d];
				// no need to search in halos here
				p->next[i][dir] = findsite(p, x, 0);
			}
			// negative dirs
			long prev = p->prev[i][dir];
			if (prev >= p->sites_total) {
				for (int d=0; d<p->dim; d++) x[d] = p->coords[prev][d];

				p->prev[i][dir] = findsite(p, x, 0);
			}
		}
	}

	// tables done and self halos moved to indices > p->sites_total. realloc in layout()
}


/* Set checkerboard parity for all sites (excl. self halos,
* which are assumed to be mapped out already).
* Parity of a site is defined to be EVEN if the sum of
* its physical coordinates x + y + z + ... is an even number,
* ODD otherwise. These are chars defined in a header file (where?)
*/
void set_parity(params *p) {

	long tot;
	p->evensites = 0; p->oddsites = 0;
	p->evenhalos = 0; p->oddhalos = 0;

	// real sites
	for (long i=0; i<p->sites; i++) {
		tot = 0;
		for (int dir=0; dir<p->dim; dir++) {
			tot += p->coords[i][dir];
		}
		if (tot % 2 == 0) {
			p->parity[i] = EVEN;
			p->evensites++;
		} else {
			p->parity[i] = ODD;
			p->oddsites++;
		}
	}
	// halo sites
	for (long i=p->sites; i<p->sites_total; i++) {
		tot = 0;
		for (int dir=0; dir<p->dim; dir++) {
			tot += p->coords[i][dir];
		}
		if (tot % 2 == 0) {
			p->parity[i] = EVEN;
			p->evenhalos++;
		} else {
			p->parity[i] = ODD;
			p->oddhalos++;
		}
	}

}


/* Quick routine for finding site index from given
* physical coordinates. Assumes ordered p->coords.
*/
long findsite(params* p, long* x, int include_halos) {

	long max;
	if (include_halos) {
		max = p->sites_total;
	} else {
		max = p->sites;
	}

	for (long i=0; i<max; i++) {
		int match = 1;
		for (int dir=0; dir<p->dim; dir++) {
			if (x[dir] != p->coords[i][dir]) {
				match = 0;
				break;
			}
		}

		if (match) {
			return i;
		}
	}

	// no match after searching through all sites?
	printf("Node %d: Failed to find matching site in p->coords!\n", p->rank);
	return -1;
}


/* Test that our mapping from site index to physical coordinates makes sense.
* This checks that basic indexing is OK, but does not check that neighboring sites
* have adjacent indices (which they should, apart from sites at hypercube sides).
*/
void test_coords(params p) {

	long* xnode = malloc(p.dim * sizeof(*xnode));

	indexToCoords(p.dim, p.nslices, p.rank, xnode);

	int dir;
	long i;


	// first site on the node?
	// we remapped sites by parity, so first site is either i=0 or i=p.evensites
	for (dir=0; dir<p.dim; dir++) {
		if (p.coords[0][dir] != xnode[dir] * p.sliceL[dir] && p.coords[p.evensites][dir] != xnode[dir] * p.sliceL[dir]) {
			printf("Node %d: Error in test_coords! First site not where it should be \n", p.rank);
			die(-111);
		}
	}
	// last real site on the node?
	for (dir=0; dir<p.dim; dir++) {
		if (p.coords[p.sites-1][dir] != (xnode[dir] + 1)* p.sliceL[dir] - 1 && p.coords[p.evensites-1][dir] != (xnode[dir] + 1)* p.sliceL[dir] - 1) {
			printf("Node %d: Error in test_coords! Last site not where it should be \n", p.rank);
			die(-112);
		}
	}


	// non-halo sites: coordinate x[j] should be in range
	// xnode[j] * p.sliceL[j] <= x[j] <= xnode[j] * p.sliceL[j] + p.sliceL[j] - 1
	for (i=0; i<p.sites; i++) {
		for (dir=0; dir<p.dim; dir++) {
			if (p.coords[i][dir] > (xnode[dir] + 1) * p.sliceL[dir] - 1 || xnode[dir] * p.sliceL[dir] > p.coords[i][dir]) {
				printf("Node %d: Error in test_coords! Real site %ld not indexed properly \n", p.rank, i);
				die(-113);
			}
		}
	}

	// halo sites: should not find a halo site whose all coordinates x[j]
	// match a real site on our node
	for (i=p.sites; i<p.sites_total; i++) {
		int ok = 0;
		for (dir=0; dir<p.dim; dir++) {
			if (p.coords[i][dir] > (xnode[dir] + 1) * p.sliceL[dir] - 1 || xnode[dir] * p.sliceL[dir] > p.coords[i][dir]) {
				ok = 1;
			}
		}
		if (!ok) {
			printf("Node %d: Error in test_coords! Halo site %ld not indexed properly \n", p.rank, i);
			die(-114);
		}
	}

	free(xnode);
}

/* Perform strong checks on p.next and p.prev, and p.parity.
* Uses p->coords to check that the physical coordinates of the neighbor sites
* are what they should for each site and direction, including halos.
*/
void test_neighbors(params p) {

	int skip_parity = 0;
	for (int dir=0; dir<p.dim; dir++) {
		if (p.L[dir] % 2 != 0) {
			// side lengths not even numbers, skip parity checks...
			skip_parity = 1;
		}
	}
	// first site in master node should have even parity
	if (!skip_parity && p.rank == 0) {
		if (p.parity[0] != EVEN) {
			printf("Node 0: Error in test_neighbors! First site parity not is not EVEN \n");
			// don't die
		}
	}

	long x[p.dim];

	for (long i=0; i<p.sites_total; i++) {

		// first check neighbors in the positive directions
		for (int dir=0; dir<p.dim; dir++) {
			long next = p.next[i][dir];
			if (next == -1) {
				if (i < p.sites) {
					// no neighbor assigned to a real lattice site, error!
					printf("Node %d: Error in test_neighbors, site %ld! Next site in direction %d not assigned \n", p.rank, i, dir);
					die(-116);
				}
			} else {
				// check that coordinates of the next site make sense.
				// What we think they should be:
				for (int d=0; d<p.dim; d++) {
					x[d] = p.coords[i][d];
				}
				x[dir]++;
				// Periodicity?
				if (x[dir] >= p.L[dir]) {
					x[dir] = 0;
				}

				// what they are according to p->coords
				for (int d=0; d<p.dim; d++) {
					if (p.coords[next][d] != x[d]) {
						printf("Node %d: Error in test_neighbors, site %ld! Coordinates of next site in direction %d do not match! \n", p.rank, i, dir);
						die(-117);
					}
				}

				// parity check
				if (!skip_parity) {
					if (p.parity[i] == p.parity[next]) {
						printf("Node %d: Error in test_neighbors, site %ld! Next site in direction %d has same parity! \n", p.rank, i, dir);
					}
				}
			}
		}

		// then same checks for negative directions
		for (int dir=0; dir<p.dim; dir++) {
			long prev = p.prev[i][dir];
			if (prev == -1) {
				if (i < p.sites) {
					// no neighbor assigned to a real lattice site, error!
					printf("Node %d: Error in test_neighbors, site %ld! Previous site in direction %d not assigned \n", p.rank, i, dir);
					die(-118);
				}
			} else {
				// check that coordinates of the next site make sense.
				// What we think they should be:
				for (int d=0; d<p.dim; d++) {
					x[d] = p.coords[i][d];
				}
				x[dir]--;
				// Periodicity?
				if (x[dir] < 0) {
					x[dir] = p.L[dir] - 1;
				}

				// what they are according to p->coords
				for (int d=0; d<p.dim; d++) {
					if (p.coords[prev][d] != x[d]) {
						printf("Node %d: Error in test_neighbors, site %ld! Coordinates of previous site in direction %d do not match! \n", p.rank, i, dir);
						die(-119);
					}
				}

				// parity check
				if (!skip_parity) {
					if (p.parity[i] == p.parity[prev]) {
						printf("Node %d: Error in test_neighbors, site %ld! Previous site in direction %d has same parity! \n", p.rank, i, dir);
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
void printf0(params p, char *msg, ...) {
  va_list fmtargs;
  va_start(fmtargs, msg);

  if(!p.rank)
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

void layout(params *p, comlist_struct *comlist) {

	p->rank = 0;
	p->size = 1;
	p->sites = p->vol;
	p->sites_total = p->vol;
	p->halos = 0;

	p->sliceL = malloc(p->dim * sizeof(*p->sliceL));
	p->nslices = malloc(p->dim * sizeof(*p->nslices));

	for (int dir=0; dir<p->dim; dir++) {
		p->sliceL[dir] = p->L[dir];
		p->nslices[dir] = 1;
	}

	alloc_lattice_arrays(p, p->sites_total);

	// construct lookup tables for site neighbors and parity
	sitemap(p);
	set_parity(p);

	// reorder by parity
	p->reorder_parity = 1;
	if (p->reorder_parity) {
		long newindex[p->sites_total];
		paritymap(p, newindex);
		remap_lattice_arrays(p, newindex, p->sites_total);
	}


	comlist->neighbors = 0;
	comlist->send_to = malloc(sizeof(*comlist->send_to));
	comlist->recv_from = malloc(sizeof(*comlist->recv_from));

}

void sitemap(params* p) {

	p->oddsites = 0; p->oddhalos = 0;
	p->evensites = 0; p->evenhalos = 0;

	for (long i=0; i<p->sites_total; i++) {

		// get the physical coordinates of site i and store in p.coords[i]
		indexToCoords(p->dim, p->L, i, p->coords[i]);

		for (int dir=0; dir < p->dim; dir++) {
			p->coords[i][dir]++;
			p->next[i][dir] = coordsToIndex(p->dim, p->L, p->coords[i]);
			p->coords[i][dir] -= 2;
			p->prev[i][dir] = coordsToIndex(p->dim, p->L, p->coords[i]);
			// return to the original value
			p->coords[i][dir] ++;
		}

	}

	printf("Site lookup tables constructed succesfully.\n");
}


void set_parity(params *p) {

	p->oddsites = 0; p->oddhalos = 0;
	p->evensites = 0; p->evenhalos = 0;
	long* x = malloc(p->dim * sizeof(*x));

	for (long i=0; i<p->sites_total; i++) {

		// get the physical coordinates of site i and store in x
		indexToCoords(p->dim, p->L, i, x);
		// calculate and store the parity of site i
		long coord = 0;
		for (int j=0; j<p->dim; j++) {
			coord += x[j];
		}

		if (coord % 2 == 0) {
			p->parity[i] = EVEN;
			p->evensites++;
		} else {
			p->parity[i] = ODD;
			p->oddsites++;
		}
	}

	free(x);
}


void printf0(params p, char *msg, ...) {
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
void paritymap(params* p, long* newindex) {

	long even=0, odd=0, evenhalo=0, oddhalo=0;

	// order real sites
	for (long i=0; i<p->sites; i++) {
		if (p->parity[i] == EVEN) {
			newindex[i] = even;
			even++;
		} else {
			newindex[i] = p->evensites + odd;
			odd++;
		}
	}
	// order halo sites
	for (long i=p->sites; i<p->sites_total; i++) {
		if (p->parity[i] == EVEN) {
			newindex[i] = p->sites + evenhalo;
			evenhalo++;
		} else {
			newindex[i] = p->sites + p->evenhalos + oddhalo;
			oddhalo++;
		}
	}

}

/* Reorder a given lattice table using the mapping given in newindex.
* The table should be p.dim * maxindex sized.
*/
void remap_latticetable(params* p, long** arr, long* newindex, long maxindex) {
	long** temp = alloc_latticetable(p->dim, maxindex);

	int dir;
	for (long i=0; i<maxindex; i++) {
		for (dir=0; dir<p->dim; dir++) {
			temp[i][dir] = arr[i][dir];
		}
	}

	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		for (dir=0; dir<p->dim; dir++) {
			arr[new][dir] = temp[i][dir];
		}
	}

	free_latticetable(temp);
}

// Same as remap_latticetable() but specifically for neighbor lookup tables
void remap_neighbor_table(params* p, long** arr, long* newindex, long maxindex) {
	long** temp = alloc_latticetable(p->dim, maxindex);

	int dir;
	for (long i=0; i<maxindex; i++) {
		for (dir=0; dir<p->dim; dir++) {
			temp[i][dir] = arr[i][dir];
		}
	}
	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		for (dir=0; dir<p->dim; dir++) {
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
void remap_lattice_arrays(params* p, long* newindex, long maxindex) {

	// backup parity
	char par[maxindex];
	for (long i=0; i<maxindex; i++) {
		par[i] = p->parity[i];
	}

	// remap
	for (long i=0; i<maxindex; i++) {
		long new = newindex[i];
		p->parity[new] = par[i];
	}

	remap_neighbor_table(p, p->next, newindex, maxindex);
	remap_neighbor_table(p, p->prev, newindex, maxindex);

	remap_latticetable(p, p->coords, newindex, maxindex);
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
* This obviously only works in a "natural" ordering of site indices
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
int coordsToRank(params p, long* coords) {
	// first find the (x, y, z, ...) coordinates of the node
	long* x = malloc( (p.dim) * sizeof(x));
	for (int dir=0; dir<p.dim; dir++) {
		x[dir] = floor(coords[dir] / p.sliceL[dir]);
	}
	// then convert the coordinates to node index
	int i = coordsToIndex(p.dim, p.nslices, x);
	free(x);

	return i;
}

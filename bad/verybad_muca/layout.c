/** @file layout.c
*
* Routines for MPI parallelization on a hypercubic (=rectangle) lattice
*	with periodic boundaries.
*
*	Future improvement: it might be simpler to perform layouting and comms
* using MPI topologies?
*/

#include "su2.h"
#include "comms.h"

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

	long maxindex = p->sites + p->halos;

	// some work arrays:
	long* newindex = malloc(maxindex * sizeof(newindex));
	// (x,y,z,...) coords on the full lattice
	long** xphys = alloc_latticetable(p->dim, p->sites + p->halos);

	// calculate xphys and rearrange sites. This also sets p->sites_total
	// and moves self halos to indices >= p->sites_total.
	sitemap(p, xphys, newindex);

	// slicing and haloing done, so can alloc lookup tables on the slices
	alloc_lattice_arrays(p);

	// construct lookup tables for site neighbors:
	calculate_neighbors(p, xphys, newindex);

	// site parity:
	set_parity(p, xphys);


	MPI_Barrier(MPI_COMM_WORLD);

	// construct communication tables
	make_comlists(p, comlist, xphys);

	// Run sanity checks on lattice layout and comms?
	MPI_Barrier(MPI_COMM_WORLD);
	if (p->run_checks) {

		start = clock();

		test_xphys(*p, xphys);
		test_neighbors(*p, xphys);
		test_comms_individual(*p, *comlist, xphys);
		test_comms(*p, *comlist);

		end = clock();
		time = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf0(*p, "All tests OK! Time taken: %lf seconds.\n", time);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	free_latticetable(xphys);
	free(newindex);

}


/* Lay out the lattice on available nodes.
*	Tricks for factorizing the lattice are mostly inspired by respective routines
* in codes of Kari and David.
* We follow Kari in the sense that the shapes of the sublattices are not fixed manually.
* In fact we should always obtain the same factorization as what Kari has in his susy code for example.
* The routine checks that the sites can always be factorized evenly among the nodes, i.e, all nodes
* get the same amount of sites. This will be assumed in all other routines in the file.
*/
void make_slices(params *p) {

	// first check that the total number of lattice sites is possible
	// to factorize in number of nodes (MPI processes)
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
	uint slice[p->dim]; // how many lattice sites in one slice in direction dir
	int nslices[p->dim]; // how many hypercubes can we fit in one direction
	int dir;
	// initialize these
	for (dir=0; dir<p->dim; dir++) {
		slice[dir] = p->L[dir];
		nslices[dir] = 1;
	}

	// start with the largest prime!
	for (int n=n_primes-1; n>=0; n--) {
		// whatever we do for one prime, repeat as many times as it appears in node factorization
		for (int f=0; f<nfactors[n]; f++) {

			// find largest slice direction that we can still chop with this prime.
			int largest=1;
			for (dir=0; dir<p->dim; dir++) {
				if ( (slice[dir] > largest) && (slice[dir] % prime[n] == 0) ) {
					largest = slice[dir];
				}
			}

			// now chop the longest direction
			for (dir=0; dir<p->dim; dir++) {
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

	// output the obtained processor layout and node sizes
	if (!p->rank) {
		printf("Processor layout: ");
		for (dir=0; dir<p->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%d", nslices[dir]);
		}
		printf("\nSites on each node: ");
		for (dir=0; dir<p->dim; dir++) {
			if (dir>0) printf(" x ");
			printf("%d", slice[dir]);
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
	printf0(*p, "Each node needs %lu additional halo sites (self halos to be removed).\n", halosites);

	p->sites = sites;
	p->halos = halosites;
}


/* Find mapping between site indices in my node and the physical (x, y, z...)
* coordinates on the full lattice, and arrange sites accordingly.
* Note that on a periodic lattice some halo sites may actually live in the same node as real sites.
*
* Specifically, the routine arranges sites so that indices of real sites
* come before halo site indices, and removes possible duplicates (self halos).
* This ordering is stored in newindex, which maps unordered indices on the full node
* to ordered indices.
*
* The coordinate mapping is stored in xphys, accessed as xphys[i][dir] = physical coordinate x_dir.
*/
void sitemap(params *p, long** xphys, long* newindex) {

	int dir;
	long i;
	long maxindex = p->sites + p->halos;

	// temporary working arrays
	long* xnode = malloc(p->dim * sizeof(xnode)); // coordinates of the MPI nodes
	long* x = malloc(p->dim * sizeof(x)); // (x,y,z,...) coords on the slice
	int* siterank = malloc(maxindex * sizeof(siterank)); // rank of the node where site i resides in
	int* ishalo = malloc(maxindex * sizeof(ishalo)); // is the site a halo or not?


	// where is the node located?
	indexToCoords(p->dim, p->nslices, p->rank, xnode);

	// loop over sites in my node, plus all halos.
	// so the loop is over a hypercube of side lengths sliceL[dir]+2
	uint* haloL = malloc( p->dim * sizeof(haloL));
	for (dir=0; dir<p->dim; dir++) {
		haloL[dir] = p->sliceL[dir] + 2;
	}

	for (i=0; i<maxindex; i++) {

		ishalo[i] = 0; // innocent until proven otherwise

		// Site coordinates (x, y, z...) in the haloed node.
		indexToCoords(p->dim, haloL, i, x);

		// Next, calculate the corresponding physical (x, y, z...) on the full lattice
		// For this we need to know if the site is halo or not
		for (dir=0; dir<p->dim; dir++) {

			// is this a halo site?
			if (x[dir] == p->sliceL[dir] + 1) {
				ishalo[i] = 1;
				// we are in the next node in positive direction. Did we go over the lattice boundary?
				if (xnode[dir]+1 >= p->nslices[dir]) {
					// went over, force periodicity
					xphys[i][dir] = 0;
				}	else {
					// no periodicity.
					xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
				}
			} else if (x[dir] == 0) {
				ishalo[i] = 1;
				// we are in the previous node in this direction. Do we force periodicity?
				if (xnode[dir]-1 < 0) {
					// went over lattice boundary
					xphys[i][dir] = p->L[dir]-1;
				} else {
					// no periodicity.
					xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
				}
			} else {
				// not a halo site
				xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
			}

		} // end dir loop

		// which node does the site live in?
		siterank[i] = coordsToRank(*p, xphys[i]);
	} // end i loop, physical coordinates done.


	// how many self halo sites do we have?
	long selfhalos = 0;
	for (i=0; i<maxindex; i++) {
		if (ishalo[i] && siterank[i] == p->rank) {
			selfhalos++;
		}
	}
	// actual number of lattice sites that we need to care about:
	p->sites_total = p->sites + p->halos - selfhalos;


	// Rearrange indices: non-halo sites should come before halos.
	// In case of self halos, move those to the end.
	long j = 0;
	long k = 0;
	long n = 0;
	for (i=0; i<maxindex; i++) {
		if (ishalo[i] && siterank[i] == p->rank) {
			// self halo site, so map it to something we don't care about anymore (could also just realloc away)
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

	// before rearranging, copy xphys into temp array
	long** xphys_temp = alloc_latticetable(p->dim, maxindex);
	for (i=0; i<maxindex; i++) {
		for (dir=0; dir<p->dim; dir++) {
			xphys_temp[i][dir] = xphys[i][dir];
		}
	}

	for (i=0; i<maxindex; i++) {
		for (dir=0; dir<p->dim; dir++) {
			xphys[newindex[i]][dir] = xphys_temp[i][dir];
		}
	}


	// free temp arrays
	free_latticetable(xphys_temp);
	free(x);
	free(xnode);
	free(ishalo);
	free(haloL);
	free(siterank);
}


/* Construct lookup tables for site neighbors, including halos.
* This is first done in the full haloed node, and the indices are
* then mapped to new, ordered indices using the mapping given in newindex.
* References to self halos are remapped to real sites instead,
* so that the final p.next and p.prev only contain real sites + actual halos.
*/
void calculate_neighbors(params *p, long** xphys, long* newindex) {

	long maxindex = p->sites + p->halos;
	long i;
	int dir;
	// work arrays
	long** nextsite = alloc_latticetable(p->dim, maxindex);
	long** prevsite = alloc_latticetable(p->dim, maxindex);
	long* x = malloc(p->dim * sizeof(*x)); // (x,y,z,...) coords on the slice

	// haloed node
	uint* haloL = malloc(p->dim * sizeof(haloL));
	for (dir=0; dir<p->dim; dir++) {
		haloL[dir] = p->sliceL[dir] + 2;
	}

	// find neighbors in the "natural" index ordering
	for (i=0; i<maxindex; i++) {

		indexToCoords(p->dim, haloL, i, x);

		for (dir=0; dir<p->dim; dir++) {

			// positive directions
			if (x[dir] + 1 >= haloL[dir]) {
				// went over the boundary, so just set unrealistic value.
				// trying to access a field at this index will give an error, so we know something went wrong!
				nextsite[i][dir] = -1;
			} else {
				x[dir]++;
				nextsite[i][dir] = coordsToIndex(p->dim, haloL, x);
				x[dir]--;
			}
			// negative directions
			if (x[dir] - 1 < 0) {
				// went over the halo boundary
				prevsite[i][dir] = -1;
			} else {
				x[dir]--;
				prevsite[i][dir] = coordsToIndex(p->dim, haloL, x);
				x[dir]++;
			}

		}
	}

	// now use newindex[] to order these similarly to xphys in sitemap(), and store in p.next/p.prev
	long next, prev;
	for (i=0; i<maxindex; i++) {
		if (newindex[i] >= p->sites_total) {
			// skip self halos
			continue;
		}
		for (dir=0; dir<p->dim; dir++) {

			// positive directions
			next = nextsite[i][dir];
			if (next == -1) {
				p->next[newindex[i]][dir] = -1;
			}
			else if (newindex[next] >= p->sites_total) {
				// the neighbor is a self halo site,
				// so use xphys to find the real site corresponding to it.
				// coords of the site whose neighbor we are looking for?
				for (int d=0; d<p->dim; d++) {
					x[d] = xphys[newindex[i]][d];
				}
				// coords of the neighbor?
				if (x[dir] + 1 >= p->L[dir]) {
					x[dir] = 0;
				} else {
					x[dir]++;
				}

				// no need to look in halos here
				p->next[newindex[i]][dir] = findsite(*p, x, xphys, 0);

			} else {
				p->next[newindex[i]][dir] = newindex[next];
			}

			// negative directions
			prev = prevsite[i][dir];
			if (prev == -1) {
				p->prev[newindex[i]][dir] = -1;
			} else if (newindex[prev] >= p->sites_total) {
				// neighbor is a self halo site
				for (int d=0; d<p->dim; d++) {
					x[d] = xphys[newindex[i]][d];
				}
				// coords of the neighbor?
				if (x[dir] - 1 < 0) {
					x[dir] = p->L[dir] - 1;
				} else {
					x[dir]--;
				}

				// no need to look in halos here
				p->prev[newindex[i]][dir] = findsite(*p, x, xphys, 0);

			}	else {
				p->prev[newindex[i]][dir] = newindex[prev];
			}
		}
	}

	free_latticetable(prevsite);
	free_latticetable(nextsite);
	free(x);
	free(haloL);
}


/* Set checkerboard parity for all sites (excl. self halos).
* Parity of a site is defined to be EVEN if the sum of
* its physical coordinates x + y + z + ... is an even number,
* ODD otherwise. These are chars defined in a header file (where?)
*/
void set_parity(params *p, long** xphys) {

	long tot;
	p->evensites = 0;
	p->oddsites = 0;

	for (long i=0; i<p->sites_total; i++) {
		tot = 0;
		for (int dir=0; dir<p->dim; dir++) {
			tot += xphys[i][dir];
		}
		if (tot % 2 == 0) {
			p->parity[i] = EVEN;
			p->evensites++;
		} else {
			p->parity[i] = ODD;
			p->oddsites++;
		}
	}

}


/* Quick routine for finding site index from given
* physical coordinates. Assumes ordered xphys.
*/
long findsite(params p, long* x, long** xphys, int include_halos) {

	long max;

	if (include_halos) {
		max = p.sites_total;
	} else {
		max = p.sites;
	}

	for (long i=0; i<max; i++) {
		int match = 1;
		for (int dir=0; dir<p.dim; dir++) {
			if (x[dir] != xphys[i][dir]) {
				match = 0;
				break;
			}
		}

		if (match) {
			return i;
		}
	}

	// no match after searching through all sites?
	printf("Node %d: Failed to find matching site in xphys!\n", p.rank);
	return -1;
}


/* Test that our mapping from site index to physical coordinates makes sense.
* This checks that basic indexing is OK, but does not check that neighboring sites
* have adjacent indices (which they should, apart from sites at hypercube sides).
*/
void test_xphys(params p, long** xphys) {

	long* xnode = malloc(p.dim * sizeof(*xnode));

	indexToCoords(p.dim, p.nslices, p.rank, xnode);

	int dir;
	long i;

	// first site on the node?
	for (dir=0; dir<p.dim; dir++) {
		if (xphys[0][dir] != xnode[dir] * p.sliceL[dir]) {
			printf("Node %d: Error in test_xphys! First site not where it should be \n", p.rank);
			die(-111);
		}
	}
	// last real site on the node?
	for (dir=0; dir<p.dim; dir++) {
		if (xphys[p.sites - 1][dir] != (xnode[dir] + 1)* p.sliceL[dir] - 1) {
			printf("Node %d: Error in test_xphys! Last site not where it should be \n", p.rank);
			die(-112);
		}
	}

	// non-halo sites: coordinate x[j] should be in range
	// xnode[j] * p.sliceL[j] <= x[j] <= xnode[j] * p.sliceL[j] + p.sliceL[j] - 1
	for (i=0; i<p.sites; i++) {
		for (dir=0; dir<p.dim; dir++) {
			if (xphys[i][dir] > (xnode[dir] + 1) * p.sliceL[dir] - 1 || xnode[dir] * p.sliceL[dir] > xphys[i][dir]) {
				printf("Node %d: Error in test_xphys! Real site %ld not indexed properly \n", p.rank, i);
				die(-113);
			}
		}
	}

	// halo sites: should not find a halo site whose all coordinates x[j]
	// match a real site on our node
	for (i=p.sites; i<p.sites_total; i++) {
		int ok = 0;
		for (dir=0; dir<p.dim; dir++) {
			if (xphys[i][dir] > (xnode[dir] + 1) * p.sliceL[dir] - 1 || xnode[dir] * p.sliceL[dir] > xphys[i][dir]) {
				ok = 1;
			}
		}
		if (!ok) {
			printf("Node %d: Error in test_xphys! Halo site %ld not indexed properly \n", p.rank, i);
			die(-114);
		}
	}

	free(xnode);
}

/* Test for p.next and p.prev, and p.parity.
* Uses xphys to check that the physical coordinates of the neighbor sites
* are what they should for each site and direction, including halos.
*/
void test_neighbors(params p, long** xphys) {

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
			// don't die; this may not be fatal
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
					x[d] = xphys[i][d];
				}
				x[dir]++;
				// Periodicity?
				if (x[dir] >= p.L[dir]) {
					x[dir] = 0;
				}

				// what they are according to xphys
				for (int d=0; d<p.dim; d++) {
					if (xphys[next][d] != x[d]) {
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
					x[d] = xphys[i][d];
				}
				x[dir]--;
				// Periodicity?
				if (x[dir] < 0) {
					x[dir] = p.L[dir] - 1;
				}

				// what they are according to xphys
				for (int d=0; d<p.dim; d++) {
					if (xphys[prev][d] != x[d]) {
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

	p->sliceL = malloc(p->dim * sizeof(p->sliceL));
	p->nslices = malloc(p->dim * sizeof(p->nslices));

	for (int dir=0; dir<p->dim; dir++) {
		p->sliceL[dir] = p->L[dir];
		p->nslices[dir] = 1;
	}

	alloc_lattice_arrays(p);

	long** xphys;
	long* newindex;
	// construct lookup tables for site neighbors and parity:
	calculate_neighbors(p, xphys, newindex);

	comlist->neighbors = 0;
	comlist->send_to = malloc(sizeof(*comlist->send_to));
	comlist->recv_from = malloc(sizeof(*comlist->recv_from));
}

// arguments xphys and newindex are not actually used in serial
void calculate_neighbors(params* p, long** xphys, long* newindex) {

	p->oddsites = 0;
	p->evensites = 0;
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

		for (int dir=0; dir < p->dim; dir++) {
			x[dir]++;
			p->next[i][dir] = coordsToIndex(p->dim, p->L, x);
			x[dir] -= 2;
			p->prev[i][dir] = coordsToIndex(p->dim, p->L, x);
			// return to the original value
			x[dir]++;
		}

	}

	free(x);

	printf("Site lookup table constructed succesfully.\n");
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


/***********************************************************************
* 	Routines for converting from site index to coordinates vice versa  *
************************************************************************/

// Calculate product L[0]*L[1]*...*L[max-1]. In practice L is the lattice or slice length.
inline long Lprod(uint* L, int max) {
	long res = 1;
	for (int i=0; i<max; i++) {
		res *= L[i];
	}

	return res;
}

/* From the site index, get the physical coordinates (x, y, z,...) on a dim-dimensional
* lattice with side lengths defined as in L, and store in x.
* Can be used for the full lattice, or a single slice!
* see Mathematica notebook coordinates.nb for analytical relations
*/
void indexToCoords(ushort dim, uint* L, long i, long* x) {

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
		x[dir] = (uint)floor((i - a) / Lprod(L, dir));
	}

}


// Same as indexToCoords(), but convert from the (x, y, z, ...) to site index.
long coordsToIndex(ushort dim, uint* L, long* x) {

		long res = 0;
		for (int dir=0; dir<dim; dir++) {
			res += ( (x[dir] + L[dir]) % L[dir] ) * Lprod(L, dir);
		}
		return res;
}

/* Convert physical (x, y, z, ...) coordinates on the full lattice
* to MPI node index (which is the same as MPI rank in this layout).
*/
int coordsToRank(params p, long* xphys) {
	// first find the (x, y, z, ...) coordinates of the node
	long* x = malloc( (p.dim) * sizeof(x));
	for (int dir=0; dir<p.dim; dir++) {
		x[dir] = floor(xphys[dir] / p.sliceL[dir]);
	}
	// then convert the coordinates to node index
	int i = coordsToIndex(p.dim, p.nslices, x);
	free(x);

	return i;
}

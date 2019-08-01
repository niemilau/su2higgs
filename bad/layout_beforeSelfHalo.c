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
	// these are needed for make_slices():
	p->sliceL = malloc(p->dim * sizeof(p->sliceL));
	p->nslices = malloc(p->dim * sizeof(p->nslices));

	make_slices(p);
	// slicing done, so can alloc lookup tables on the slices
	alloc_lattice_arrays(p);

	// (x,y,z,...) coords on the full lattice
	long** xphys = alloc_latticetable(p->dim, p->sites + p->halos);

	/*
	// test memory
	long j = 0;
	for (long i=0; i<10; i++) {
		for (int dir=0; dir<p->dim; dir++) {
			xphys[i][dir] = j;
			j++;
		}
	}
	long *ptr = &xphys[0][0];
	for (long i=0; i<10; i++) {
			printf0(*p, "Memory address %d, value %ld \n", xphys+i, ptr[i]);
	}
	*/

	calculate_neighbors(p, xphys);
	// p.prev, p.next, p.parity and xphys done,
	// and sites indexed so that real sites come before halos.

	// Run sanity checks on lattice layout?
	MPI_Barrier(MPI_COMM_WORLD);
	if (p->run_checks) {
		if (p->dim == 2) {
			MPI_Barrier(MPI_COMM_WORLD);
			printf0(*p, "Printing visual description of lattice coordinates, enjoy!\n");
			printf0(*p, "Outer sites are halos.\n");
			//print_lattice_2D(p, xphys); // does not work with self halos!!!
		}
		test_layout(*p);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// construct communication tables
	make_comlists(p, comlist, xphys);

	MPI_Barrier(MPI_COMM_WORLD);
	free_latticetable(p->sites + p->halos, xphys);

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

	/*
	// Debug: show factorization
  if(!p->rank) {
    for(int n=0;n<n_primes;n++)
      printf("nfactors %d = %d\n", prime[n], nfactors[n]);
  }
	*/

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

	// begin with the largest prime!
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
	* This is implemented in calculate_neighbors() below where we specify neighbors and halos for lattice sites.
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
	printf0(*p, "Each node needs %lu additional halo sites.\n", halosites);

	p->sites = sites;
	p->halos = halosites;
	p->sites_total = sites + halosites;
}



/* Long routine for initializing lookup tables for neighboring sites and parity.
* Halo sites are given an "artificial" index i >= p.sites that is stored
* normally in p.next (or p.prev).
* Special care is required for sites at node boundary,
* and especially at the "ends" of the physical lattice,
* where we need to ensure periodicity.
* Physical (x, y, z, ...) coordinates on the full lattice are also calculated and stored in xphys.
*
* I can think of 2 complicated ways to implement this geometry; this function uses
* the more cryptic method that includes halos at the hypercube corners easily.
* The other method is used in test_layout(), which checks that what we do here makes sense.
*/
void calculate_neighbors(params *p, long** xphys) {

	int dir;
	long i;
	long maxindex = p->sites_total;

	// temporary working arrays
	long* xnode = malloc(p->dim * sizeof(xnode)); // coordinates of the MPI nodes
	long* x = malloc(p->dim * sizeof(x)); // (x,y,z,...) coords on the slice
	long** nextsite = alloc_latticetable(p->dim, maxindex);
	long** prevsite = alloc_latticetable(p->dim, maxindex);
	int* siterank = malloc(maxindex * sizeof(siterank)); // rank of the node where site i resides in
	int* ishalo = malloc(maxindex * sizeof(ishalo)); // is the site a halo or not?

	// location of my node?
	node_neighbors(p, xnode);

	// loop over sites in my node, plus all halos.
	// so the loop is over a hypercube of side lengths sliceL[dir]+2
	uint* haloL = malloc( p->dim * sizeof(haloL));
	for (dir=0; dir<p->dim; dir++) {
		haloL[dir] = p->sliceL[dir] + 2;
	}

	/* Coordinates are stored in:
	* xphys[i] for coords on the full lattice
	*	xnode for coords of the node
	* x for site coordinates on the node
	*/

	// Goal: find the physical coordinates of all sites in the node, including halo sites.
	// These are then used to figure out which MPI node the halo lies in, and to calculate site parity.
	for (i=0; i<maxindex; i++) {

		ishalo[i] = 0; // innocent until proven otherwise

		// Site coordinates (x, y, z...) in the haloed node.
		indexToCoords(p->dim, haloL, i, x);

		// Next, calculate the corresponding physical (x, y, z...) on the full lattice
		// For this we need to know if the site is halo or not
		long coord = 0;

		for (dir=0; dir<p->dim; dir++) {

			// is this a halo site?
			if (x[dir] == p->sliceL[dir] + 1) {
				ishalo[i] = 1;
				// we are in the next node in positive direction. Did we go over the lattice boundary?
				if (xnode[dir]+1 >= p->nslices[dir]) {
					// went over, force periodicity
					xphys[i][dir] = 0;
				}	else {
					// no periodicity. Triple check that this works in corners too!!
					xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
				}
			} else if (x[dir] == 0) {
				ishalo[i] = 1;
				// we are in the previous node in this direction. Do we force periodicity?
				if (xnode[dir]-1 < 0) {
					// went over lattice boundary
					xphys[i][dir] = p->L[dir]-1;
				} else {
					// no periodicity. Triple check that this works in corners too!!
					xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
				}
			} else {
				// not a halo site
				xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
			}
			coord += xphys[i][dir];

		} // end dir loop
		// physical coordinates done.

		// Parity:
		p->parity[i] = coord % 2;
		// which node does the site live in?
		siterank[i] = coordsToRank(*p, xphys[i]);


		/* Find neighboring sites with halos included.
		* Here neighbors over the halo boundary don't really matter as we shouldn't need them.
		* Note however that on a periodic lattice some halo sites may actually live in the same node as real sites.
		* We deal with these by simply referring p.next and p.prev to the real site instead of the "self halo",
		* and not touching them from there on.
		*/
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
		// neighbor lookup tables done

	} // end i loop

	// Rearrange indices: non-halo sites should come before halos
	long* newindex = malloc(maxindex * sizeof(newindex));
	long j = 0;
	long k = 0;
	for (i=0; i<maxindex; i++) {
		if (!ishalo[i]) {
			newindex[i] = j;
			j++;
		} else {
			newindex[i] = p->sites + k;
			k++;
		}
	}
	// before rearranging, copy parity and xphys into temp array
	char* parity_temp = malloc(maxindex * sizeof(parity_temp));
	memcpy(parity_temp, p->parity, maxindex * sizeof(parity_temp));
	long** xphys_temp = alloc_latticetable(p->dim, maxindex);
	for (i=0; i<maxindex; i++) {
		for (dir=0; dir<p->dim; dir++) {
			xphys_temp[i][dir] = xphys[i][dir];
		}
	}

	for (i=0; i<maxindex; i++) {
		p->parity[newindex[i]] = parity_temp[i];

		for (dir=0; dir<p->dim; dir++) {
			xphys[newindex[i]][dir] = xphys_temp[i][dir];
			// next/prev sites
			if (nextsite[i][dir] == -1) {
				p->next[newindex[i]][dir] = -1;
			} else {
				p->next[newindex[i]][dir] = newindex[ nextsite[i][dir] ];
			}
			if (prevsite[i][dir] == -1) {
				p->prev[newindex[i]][dir] = -1;
			} else {
				p->prev[newindex[i]][dir] = newindex[ prevsite[i][dir] ];
			}
		}
	}

	// Self halos. We need rearranged siterank[]; easire to just recalculate this.
	for (i=0; i<maxindex; i++) {
		siterank[i] = coordsToRank(*p, xphys[i]);
	}

	long next, prev;
	for (i=0; i<maxindex; i++) {

		for (dir=0; dir<p->dim; dir++) {
			next = p->next[i][dir];
			prev = p->prev[i][dir];
			// after rearranging, halos indices are >= p->sites
			if (next >= p->sites && siterank[next] == p->rank) {
				// find the real site from xphys
				for (long j=0; j<p->sites; j++) {
					int match = 1;
					for (int d=0; d<p->dim; d++) {
						if (xphys[next][d] != xphys[j][d]) {
							match = 0;
						}
					}
					if (match) {
						p->next[i][dir] = j;
						break;
					}
				}
			}
			// same for previous
			if (prev >= p->sites && siterank[prev] == p->rank) {
				// find the real site from xphys
				for (long j=0; j<p->sites; j++) {
					int match = 1;
					for (int d=0; d<p->dim; d++) {
						if (xphys[prev][d] != xphys[j][d]) {
							match = 0;
						}
					}
					if (match) {
						p->prev[i][dir] = j;
						break;
					}
				}
			}

		}
	}

	for (int rank=0; rank<p->size; rank++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (p->rank == rank) {
			printf("This is node %d\n", rank);
			for (i=0; i<maxindex; i++) {
				printf("i = %ld, prev[i][0] = %ld, prev[i][1] = %ld, next[i][0] = %ld, next[i][1] = %ld\n", i, p->prev[i][0], p->prev[i][1], p->next[i][0], p->next[i][1]);
				}
		}
	}
	// p.parity, p.next, p.prev now ready and arranged correctly, shouldn't touch these anymore!

	free(parity_temp);
	free_latticetable(p->sites + p->halos, xphys_temp);


	// free memory before finishing routine
	free_latticetable(p->halos + p->sites, prevsite);
	free_latticetable(p->halos + p->sites, nextsite);
	free(x);
	free(xnode);
	free(ishalo);
	free(haloL);
	free(siterank);
	free(newindex);

}


/* Routine for finding the ranks of all neighboring nodes
* and storing them in p.nextrank, p.prevrank.
* TODO: CORNERS ARE MISSING, but not strictly needed...
* Nodes are laid out in the same fashion as lattice sites,
*	so same routines can be used to figure out their physical location.
* Coordinates of the node are stored in xnode.
*/
void node_neighbors(params *p, long *xnode) {

	int dir;
	long i, oldx;

	// where is the node located?
	indexToCoords(p->dim, p->nslices, p->rank, xnode);
	// who are the neighbors?
	for (dir=0; dir<p->dim; dir++) {

		// positive directions
		oldx = xnode[dir];
		xnode[dir]++;
		// did we go over the lattice boundary?
		if (xnode[dir] >= p->nslices[dir]) {
			// went over, so account for peridicity
			xnode[dir] = 0;
		}
		p->nextrank[dir] = coordsToIndex(p->dim, p->nslices, xnode);
		// move back to where we were
		xnode[dir] = oldx;

		// same for negative directions
		if (xnode[dir] - 1 < 0) {
			// going over, so account for peridicity
			xnode[dir] = p->nslices[dir]-1;
		} else {
			xnode[dir]--;
		}
		p->prevrank[dir] = coordsToIndex(p->dim, p->nslices, xnode);
		// move back to where we were
		xnode[dir] = oldx;
	}
	// lookup table for neighboring MPI nodes ready

}


/* Routine for running some sanity checks on lattice layout.
* Assumes that sites are indiced "canonically", e.g. all real sites first, halos after, in the usual ordering.
* This routine essentially constructs the physical (x, y, z, ...) coordinates starting from the site indices,
* and uses those to check p.prev, p.next and p.parity, as well as halo ordering.
* Could replace the method in calculate_neighbors() if desired, but right now is just used for crosschecking.
*/
void test_layout(params p) {

	clock_t start, end;

	start = clock();
	printf0(p, "Running sanity checks on lattice layout in %d dimensions.\n", p.dim);

	long maxindex = p.sites + p.halos;
	long* xnode = malloc(p.dim * sizeof(xnode)); // coordinates of the MPI nodes
	long* x = malloc(p.dim * sizeof(x)); // (x,y,z,...) coords on the slice, without halos this time

	indexToCoords(p.dim, p.nslices, p.rank, xnode);

	long nextsite, prevsite;
	long* nextcoords = malloc(p.dim * sizeof(nextcoords));
	long** xphys = alloc_latticetable(p.dim, maxindex); // (x,y,z,...) coords on the full lattice

	long i;
	char par;
	char pass = 1;

	int skip_parity = 0;
	for (int dir=0; dir<p.dim; dir++) {
		if (p.L[dir] % 2 != 0) {
			printf0(p, "Lattice side lengths not even numbers! Skipping parity checks.\n");
			skip_parity = 1;
		}
	}

	// first find the physical coordinates on the full lattice.
	// this is easy because we can ignore the halos
	for (i=0; i<p.sites; i++) {
		indexToCoords(p.dim, p.sliceL, i, x);
		for (int dir=0; dir<p.dim; dir++) {
			xphys[i][dir] = x[dir] + p.sliceL[dir] * xnode[dir];
		}
	}


	// Run checks on the real lattice sites, ignoring halos
	for (i=0; i<p.sites; i++) {

		for (int dir=0; dir<p.dim; dir++) {
			nextsite = p.next[i][dir];
			prevsite = p.prev[i][dir];
			if (nextsite < 0 || prevsite < 0) {
				// went over the node boundary. But this can't happen in a loop over the real sites!
				printf("Node %d: Error at site %ld! Site not haloed properly, or ordering is wrong?\n", p.rank, i);
				pass = 0;
			}
			else {
				if (nextsite < p.sites) {
					// non-halo site, but may be self halo
					long nextphys = xphys[i][dir] + 1;
					if (nextphys >= p.L[dir]) {
						nextphys = 0;
					}
					indexToCoords(p.dim, p.sliceL, nextsite, x);
					if (nextphys != x[dir] + p.sliceL[dir] * xnode[dir]) {
						// shouldn't need to check periodicity here since nextsite is not in the halo
						printf("Node %d: Error at site %ld! Next site in positive dir %d wrong?\n", p.rank, i, dir);
						pass = 0;
					}
				} else {
					// halo site, update xphys. NB! this updates only some of the coordinates.
					// rest will be added later when we loop over the halos.
					xphys[nextsite][dir] = (xphys[i][dir] + 1) % p.L[dir];
				}
				if (prevsite < p.sites) {
					// non-halo site, but may be self halo
					long prevphys = xphys[i][dir] - 1;
					if (prevphys < 0) {
						prevphys = p.L[dir] - 1;
					}
					indexToCoords(p.dim, p.sliceL, prevsite, x);
					if (prevphys != x[dir] + p.sliceL[dir] * xnode[dir]) {
						printf("Node %d: Error at site %ld! Next site in negative direction %d wrong?\n", p.rank, i, dir);
						pass = 0;
					}
				} else {
					// update halo xphys
					if (xphys[i][dir] == 0) {
						xphys[prevsite][dir] = p.L[dir] - 1;
					}  else {
						xphys[prevsite][dir] = xphys[i][dir] - 1;
					}
				}
			}
		}
		// check if we're in a hypercube corner. If yes, find corner halo and update xphys for it
		indexToCoords(p.dim, p.sliceL, i, x);
		char corner = 1;
		long j = i;
		for (int dir=0; dir<p.dim; dir++) {
			if (x[dir] == 0) {
				j = p.prev[j][dir];
			} else if (x[dir] == p.sliceL[dir] - 1) {
				j = p.next[j][dir];
			} else {
				// not a corner
				corner = 0;
				break;
			}
		}
		if (corner) {
			for (int dir=0; dir<p.dim; dir++) {
				// need to account for periodicity here
				if (x[dir] == 0) {
					if (xphys[i][dir] == 0) {
						xphys[j][dir] = p.L[dir] - 1;
					} else {
						xphys[j][dir] = xphys[i][dir] - 1;
					}
				} else if (x[dir] == p.sliceL[dir] - 1) {
					xphys[j][dir] = (xphys[i][dir] + 1) % p.L[dir];
				}
			}
		}
		// coordinates for corner halos should now be in xphys
	} // end i loop

	/*
	*	Checks on haloing.
	*/

	i=0;
	// from the first real site, move to where we think the first halo site should be:
	for (int dir=0; dir<p.dim; dir++) {
		i = p.prev[i][dir];
	}
	if (i != p.sites) {
		printf("Node %d: Failed to reach the first halo site with p.prev! Check halo indexing.\n", p.rank);
		pass = 0;
	}
	// initialize xphys for the first halo
	for (int dir=0; dir<p.dim; dir++) {
		if (xphys[0][dir] == 0) {
			xphys[i][dir] = p.L[dir] - 1;
		} else {
			xphys[i][dir] = xphys[0][dir] - 1;
		}
		x[dir] = 0;
	}


	char skipnext, skipprev;
	// then just loop over the halos and check parity/neighbor sites, updating xphys on the move
	for (i=p.sites; i<maxindex; i++) {

		for (int dir=0; dir<p.dim; dir++) {
			nextsite = p.next[i][dir];
			prevsite = p.prev[i][dir];
			// check halo boundary and update xphys for next sites
			// but don't modify xphys for real sites.
			if (nextsite < 0) {
				skipnext = 1;
			} else {
				skipnext = 0;
				if (nextsite >= p.sites) {
					xphys[nextsite][dir] = (xphys[i][dir] + 1) % p.L[dir];
				}
			}
			if (prevsite < 0) {
				skipprev = 1;
			} else {
				skipprev = 0;
				if (prevsite >= p.sites) {
					if(xphys[i][dir] == 0) { // valgrind reported error here, xphys[i][dir] uninitialized?
						// periodicity
						xphys[prevsite][dir] = p.L[dir] - 1;
					} else {
						xphys[prevsite][dir] = xphys[i][dir] - 1;
					}
				}
			}

			// parity check:
			par = p.parity[i];
			if (skip_parity == 0 && !skipnext) {
				if (par == p.parity[nextsite]) {
					printf("Node %d: Inconsistent parity at site %ld, neighbors have same parity!\n", p.rank, i);
					pass = 0;
				}
			}
			if (skip_parity == 0 && !skipprev) {
				if (par == p.parity[prevsite]) {
					printf("Node %d: Inconsistent parity at site %ld, neighbors have same parity!\n", p.rank, i);
					pass = 0;
				}
			}
			// check that neighbors physical locations and indices agree.
			// Since we updated xphys for halos on the go, this is a strong check on halo layout.
			if (nextsite < p.sites && !skipnext) {
				// real site
				if (xphys[nextsite][dir] != (xphys[i][dir] + 1) % p.L[dir]) {
					printf0(p, "Node %d: Wrong neighbor in positive direction %d at halo site %ld!\n", p.rank, dir, i);
					pass = 0;
				}
			}
			if (prevsite < p.sites && !skipprev) {
				// real site
				long j;
				if (xphys[i][dir] == 0) {
					j = p.L[dir] - 1;
				} else {
					j = xphys[i][dir] - 1;
				}
				if (xphys[prevsite][dir] != j) {
					printf("Node %d: Wrong neighbor in negative direction %d at halo site %ld!\n", p.rank, dir, i);
					pass = 0;
				}
			}
		} // end dir loop

	} // end halo loop


	free_latticetable(p.halos + p.sites, xphys);
	free(nextcoords);
	free(xnode);
	free(x);

	end = clock();
	double t = (double) (end-start)/CLOCKS_PER_SEC;

	char allpass = pass;

	MPI_Barrier(MPI_COMM_WORLD);
	char buf;
	// send and receive to/from other nodes to see if they passed
	int tag = 42;
	for (int rank=0; rank<p.size; rank++) {
		if (p.rank != rank) {
			//printf("This is MPI process %d. I'm sending pass = %d to and receiving from process %d. \n", p.rank, pass, rank);
			MPI_Sendrecv(&pass, 1, MPI_BYTE, rank, tag, &buf, 1, MPI_BYTE, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (buf != 1) {
				allpass = 0;
				printf("Tests failed in node %d!\n", rank);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (allpass) {
		printf0(p, "All tests passed! Time taken: %lf seconds.\n", t);
	} else {
		printf0(p, "Some sanity checks failed! Check layout.c.\n");
		printf0(p, "Time taken for tests: %lf seconds\n", t);
		die(20);
	}

}


/* print 2D graphic of the sites + halos
	* Makes only sense for dim=2 really
	*/
void print_lattice_2D(params *p, long** xphys) {
	long nextindex, previndex;
	for (int rank = 0; rank < p->size; rank++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (p->rank == rank) {
			printf("This is MPI rank %d\n",p->rank);

			int dir = 0;
			// start from upper halo
			long i = p->prev[0][dir+1];
			if (i < 0) {
				printf("Unable to print 2D lattice, order the sites first! In layout.c\n");
				die(12);
			}
			long rowfirst = p->prev[i][dir];
			while (1==1) {

				// always print the previous site, apart from at last site of the row.
				// need some care here to account for self halos (TBA!)
				previndex = p->prev[i][dir];
				if (previndex != -1) {
					if (p->parity[previndex] == 0) {
						printf("-(");
						for (int d=0; d<p->dim; d++) {
							printf("%ld,", xphys[previndex][d]);
						}
						printf(")-");
					} else {
						printf("-[");
						for (int d=0; d<p->dim; d++) {
							printf("%ld,", xphys[previndex][d]);
						}
						printf("]-");
					}
				} else {
					rowfirst = i;
				}

				nextindex = p->next[i][dir];
				if (nextindex == -1 || nextindex == rowfirst) {
					nextindex = p->next[rowfirst][dir+1];

					if (p->parity[i] == 0) {
						printf("-(");
						for (int d=0; d<p->dim; d++) {
							printf("%ld,", xphys[i][d]);
						}
						printf(")-");
					} else {
						printf("-[");
						for (int d=0; d<p->dim; d++) {
							printf("%ld,", xphys[i][d]);
						}
						printf("]-");
					}

					printf("\n");
					if (nextindex == -1) {
						break;
					}
				}
				i = nextindex;
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
#warning Not using MPI!!!

void layout(params *p) {

	p->rank = 0;
	p->size = 1;
	p->sites = p->vol;
	p->halos = 0;
	// calculate neighbors here

}

void printf0(params p, char *msg, ...) {
  va_list fmtargs;
  va_start(fmtargs, msg);

  vfprintf(msg, fmtargs);
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
	x[dim-1] =  (int)floor(i / (Lprod(L, dim-1)) );


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

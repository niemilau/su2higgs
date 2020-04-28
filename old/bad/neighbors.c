/* Create lookup table for neighboring sites.
* Idea: each MPI node gets a hypercube with size defined by p.nslices[] and
* lattice sites on a node are indexed by i in range [0, p.sites-1].
* Here we do the geometrical computations to make this work.
* Special care is required for sites at node boundary,
* and especially at the "ends" of the physical lattice,
* where we need to ensure periodicity.
* For a site i at slice boundary, we define its halo to be a site with
* "artificial" index i >= p.sites and store it normally in p.next (or p.prev).
* This is convenient because then a loop over the
* physical sites is still from i=0 to i=p.sites.
* The actual halo index is figured out using the cartesian (x,y,z,...) coordinates.
* These halo indices need to be accounted for when allocating fields (in alloc.c).
* We also calculate the checkerboard parity of each site.
*/

void calculate_neighbors_old(params *p) {
	// cartesian coordinates 
	int* xnode = malloc(p->dim * sizeof(xnode)); // coordinates of the MPI nodes
	int* x = malloc(p->dim * sizeof(x)); // on the slice
	long** xphys = alloc_latticetable(p->dim, p->sites + p->halos); // on the full lattice
	
	// NB need to free these somewhere (now in main.c).
	// should put all of these into one function in malloc.c or smt
	p->halorank = malloc(p->halos * sizeof(p->halorank));
	p->halosite = malloc(p->halos * sizeof(p->halosite));


	int dir, slice;
	long i;
	long oldx;
	/* Lay out the nodes in the same fashion as lattice sites,
	* so same routines can be used to figure out their physical location.*/

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


	// PLAN: first treat the halos just as real lattice sites to find their neighbors 
	// in a natural manner. Then map site indices so that real sites always come before 
	// halos in index notation.
	
	// helper for checking if the current site is actually a halo site
	int* ishalo = malloc( (p->sites + p->halos) * sizeof(ishalo)); 
	
	// temporary lookup tables for the node + halos. 
	// we will later contiguously rearrange the contents to p->next, p->prev
	long** nextsite = alloc_latticetable(p->dim, p->sites + p->halos);
	long** prevsite = alloc_latticetable(p->dim, p->sites + p->halos);
	
	// now we loop over sites in my node, plus all halos. 
	// so the loop is over a hypercube of side lengths sliceL[dir]+2
	uint* haloL = malloc( p->dim * sizeof(haloL)); 
	for (dir=0; dir<p->dim; dir++) {
		haloL[dir] = p->sliceL[dir] + 2;
	}
	
	long maxindex = p->sites + p->halos;
	for (i=0; i<maxindex; i++) {
		ishalo[i] = 0; // innocent unless proven otherwise
		// site coordinates (x, y, z...) in the haloed node
		indexToCoords(p->dim, haloL, i, x);
		
		// Next, calculate the corresponding physical (x, y, z...) on the full lattice
		// For this we need to know if the site is halo or not 
		long coord = 0;
		
		for (dir=0; dir<p->dim; dir++) {
			
			// is this a halo site?
			if (x[dir] == p->sliceL[dir] + 1) {		
				ishalo[i] = 1;
				// we are in the next node in positive direction. Did we go over the lattice boundary?
				if (xnode[dir]+1 > p->nslices[dir]) {
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
				// not a halo site (YET!!!!), so the physical coordinates are obtained as 
				xphys[i][dir] = x[dir] + xnode[dir] * p->sliceL[dir] - 1;
			}
		}
		
		for (dir=0; dir<p->dim; dir++) {
			
			// sum of physical coordinates (for checkerboard parity)
			coord += xphys[i][dir];
			
			// find neighboring sites with halos included and store in temporary nextsite/prevsite.
			// here neighbors over the halo boundary don't really matter as we shouldn't need them.
			
			// positive directions
			if (x[dir] + 1 >= haloL[dir]) {
				// went over the boundary, so just set unrealistic value
				nextsite[i][dir] = -1;
			} else {
				x[dir]++;
				prevsite[i][dir] = coordsToIndex(p->dim, haloL, x);
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
		
		// parity
		p->parity[i] = coord % 2;
	} // end i loop
	


	/* print 2D graphic of the sites 
	* Makes only sense for dim=2 really 
	*/
		/*
	for (int rank = 0; rank < p->size; rank++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (p->rank == rank) {
			printf("This is MPI rank %d\n",p->rank);
		
			for (int y=0; y < haloL[1]; y++) {
				printf("|");
				
				for (int x1=0; x1 < haloL[0]; x1++) {
					printf("-(");
					int xx[] = {x1, y};
					for (dir=0; dir<p->dim; dir++) {
							printf("%ld,", xphys[coordsToIndex(2, haloL, xx)][dir]);
					}	
					printf(")-");
				}
				printf("\n");
			}
			
		}
	}
	*/
	
		
	// Now we rearrange the sites so that real sites on the node have index 0 < i < p.sites,
	// while halo sites have i >= p.sites. This is convenient for sweeping over the physical lattice.
	// Downside: memory addresses not contiguous (at least not for the halo parts) 
	// TODO try looping over halos too in update sweeps and just add a check: if (ishalo[i]) ?
	
	
	long* newindex = malloc( (maxindex) * sizeof(newindex));
	long haloindex = 0;
	long j = 0;
	for (i=0; i<maxindex; i++) {
		if (!ishalo[i]) {
			// physical site
			newindex[i] = j;
			j++;
		} else {
			// halo site
			newindex[i] = p->sites + haloindex;
			haloindex++;
		}
	}
	
	// update parity and site neighbors 
	for (i=0; i<maxindex; i++) {
		p->parity[newindex[i]] = p->parity[i];
		//printf("newindex = %lu\n",newindex[i]);
		for (dir=0; dir<p->dim; dir++) {		
			p->next[newindex[i]][dir] = nextsite[i][dir];
			p->prev[newindex[i]][dir] = prevsite[i][dir];
		}
	}
	
	
	/*
	// loop over sites in my node
	for (i=0; i<p->sites; i++) {
		// site coordinates (x, y, z...) in my node
		indexToCoords(p->dim, p->sliceL, i, x);
		// corresponding physical (x, y, z...) on the full lattice (for parity)
		uint coord = 0;
		for (dir=0; dir<p->dim; dir++) {
			xphys[dir] = x[dir] + xnode[dir] * p->sliceL[dir];
			coord += xphys[dir];

			/* neighbors of site i?
			* Idea: Just assign an index i >= p.sites to halo fields
			* in the order we encounter them. For later referencing,
			* figure out the real index of this site in the NEXT node
			* and store it in ??.
			* Possible problem: hard to allocate fast contiguous memory?
			

			// TODO corners!!!

			// positive directions
			oldx = x[dir];
			// do we go over the slice boundary?
			if (x[dir] + 1 >= p->sliceL[dir]) {
				// going over, so figure out site index in the next node.
				// this is easy because all of our nodes are of the same shape!
				if (haloindex >= p->halos) {
					printf("Halo index out of bounds! in layout.c.");
					die(102);
				}
				x[dir] = 0;
				newrank = p->nextrank[dir];
				newindex = coordsToIndex(p->dim, p->sliceL, x);
				p->halosite[haloindex] = newindex;
				
				p->next[i][dir] = p->sites + haloindex;
				printf0(*p,"haloindex = %lu\n",haloindex);
				haloindex++;
			} else {
				// still in the same slice
				x[dir]++;
				p->next[i][dir] = coordsToIndex(p->dim, p->sliceL, x);
			}
			// move back to where we were
			x[dir] = oldx;
			
			// repeat for negative directions
			if (x[dir] - 1 < 0) {
				// going over, so need to halo
				if (haloindex >= p->halos) {
					printf("Halo index out of bounds! in layout.c.");
					die(103);
				} 
				printf0(*p,"haloindex = %lu\n",haloindex);
				
				p->prev[i][dir] = p->sites + haloindex;	
				haloindex++;
			} else {
				// still in the same slice 
				x[dir]--;
				p->prev[i][dir] = coordsToIndex(p->dim, p->sliceL, x);
			}
			// move back to where we were
			x[dir] = oldx;
		}
		// parity:
		p->parity[i] = coord % 2;
	}
	*/

	// Do some tedious printing for debugging
	
	/*
	for (int ranks=0; ranks<p->size; ranks++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (p->rank == ranks) {
			printf("This is MPI process %d. My coordinates are:\n ", p->rank);
			for (dir=0; dir<p->dim; dir++) {
				printf("	");
				printf("xnode[%d] = %u ", dir, xnode[dir]);
			}
			printf("\nAnd ranks of my neighboring nodes are:\n ");
			for (dir=0; dir<p->dim; dir++) {
				printf("	");
				printf("nextrank[%d] = %u	 prevrank[%d] = %u \n", dir, p->nextrank[dir], dir, p->prevrank[dir]);
			}
			
			
			long phys;
			printf("\nMy sites and their physical coordinates are:\n");
			for (i=0; i<p->sites; i++) {
				indexToCoords(p->dim, p->sliceL, i, x);
				printf("i=%lu, claimed i=%lu, parity=%d, x = (", i, coordsToIndex(p->dim, p->sliceL, x), p->parity[i]);
				for (dir=0; dir<p->dim; dir++) {
					phys = x[dir] + xnode[dir] * p->sliceL[dir];
					printf("%ld,", phys);
				}
				printf("\b). Claimed x = (");
				for (dir=0; dir<p->dim; dir++) {
					printf("%ld,", xphys[newindex[i]][dir]);
				}
				printf("\b).\n");
			}
			
		}
	}
*/

	free(x);
	free_latticetable(p->dim, xphys);
	free(xnode);
	free(ishalo);
	free(haloL);
	free_latticetable(p->dim, nextsite);
	free_latticetable(p->dim, prevsite);
	free(newindex);
	
}
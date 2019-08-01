/*
*	For each site in my node, find its neighboring sites
* by brute forcing through xphys.
*
* Slow, so we apply some simple optimizations. 
* Namely, we specify an offset for the site loop 
* so that we don't search for the neighbors in hopeless locations.  
*
*	I had a faster version in old version of layout.c (see backups), 
* but it's harder to implement with self halos.
*/
void calculate_neighbors00(params *p, long** xphys) {
	
	long i, oldx, offset;
	int dir, match, matchfound;
	
	long* x = malloc(p->dim * sizeof(*x));
	// for optimization, make use of node coordinates
	long* xnode = malloc(p->dim * sizeof(*xnode));
	indexToCoords(p->dim, p->nslices, p->rank, xnode);
	
	// first neighboring sites in positive directions:
	for (i=0; i<p->sites_total; i++) {
		// my coords?
		for (dir=0; dir<p->dim; dir++) {
			x[dir] = xphys[i][dir];
		}
		for (dir=0; dir<p->dim; dir++) {
			offset = 0;
			oldx = x[dir];
			// next coordinate in direction dir?
			if (x[dir] + 1 >= p->L[dir]) {
				// periodicity, so either I am a halo, or neighbor is either halo or self halo 
				x[dir] = 0;
				if (p->nslices[dir] > 1 && i < p->sites) {
					// cannot be self halo
					offset = p->sites;
				} 
			} else {
				x[dir]++;
				// optimization:
				if (x[dir] > p->sliceL[dir] * (xnode[dir] + 1)) {
					offset = p->sites;
				}
			}
			
			matchfound = 0;
			// find matching site from xphys: 
			for (long j=offset; j<p->sites_total; j++) {
				match = 1;
				for (int d=0; d<p->dim; d++) {
					if (x[d] != xphys[j][d]) {
						match = 0;
						break;
					}
				}
				if (match) {
					// found next site from xphys
					matchfound = 1;
					p->next[i][dir] = j;
					break;
				}
			}
			if (!matchfound) {
				// no match found! next site not in my node (including my halo)
				p->next[i][dir] = -1;
			}
			
			// set x back to my own coordinate and repeat for next direction
			x[dir] = oldx;
		}
	}	// p.next done
	
	// then repeat for neighboring sites in negative directions:
	for (i=0; i<p->sites_total; i++) {
		// my coords?
		for (dir=0; dir<p->dim; dir++) {
			x[dir] = xphys[i][dir];
		}
		for (dir=0; dir<p->dim; dir++) {
			offset = 0;
			oldx = x[dir];
			// previous coordinate in direction dir?
			if (x[dir] - 1 < 0) {
				// periodicity, so either I'm a halo, or neighbor is either halo or self halo 
				x[dir] = p->L[dir] - 1;
				if (p->nslices[dir] > 1) {
					// cannot be self halo
					offset = p->sites;
				}
			} else {
				x[dir]--;
				// optimization:
				if (x[dir] < (p->sliceL[dir]) * xnode[dir]) {
					offset = p->sites;
				}
			}
			
			matchfound = 0;
			// find matching site from xphys: 
			for (long j=offset; j<p->sites_total; j++) {
				match = 1;
				for (int d=0; d<p->dim; d++) {
					if (x[d] != xphys[j][d]) {
						match = 0;
						break;
					}
				}
				if (match) {
					// found previous site from xphys
					matchfound = 1;
					p->prev[i][dir] = j;
					break;
				}
			}
			if (!matchfound) {
				// no match found! previous site not in my node (including my halo)
				p->prev[i][dir] = -1;
			}
			
			// set x back to my own coordinate and repeat for next direction
			x[dir] = oldx;
		}
	}	
	
	// test
	
	for (i=0; i<p->sites_total; i++) {
		for (dir=0; dir<p->dim; dir++) {
			printf0(*p, "xphys[%ld][%d] = %ld\n", i, dir, xphys[i][dir]);
		}
	}
	
	for (i=0; i<p->sites_total; i++) {
		for (dir=0; dir<p->dim; dir++) {
			printf0(*p, "next[%ld][%d] = %ld, prev[%ld][%d] = %ld \n", i, dir, p->next[i][dir], i, dir, p->prev[i][dir]);
		}
	}
	
	free(x);
	free(xnode);
	
}
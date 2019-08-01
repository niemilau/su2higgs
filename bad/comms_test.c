#include "su2.h"
#include "comms.h"


/* Routine for constructing MPI communication lookup tables.
* Performance goal: communications take max ~10% of all time (Kari has ~6% - 17% in the simulations I did)
*
* TODO SELF HALOS !!!!!!!!!!!!! and free comlist and its substructs
*/
void make_comlists(params *p, comlist_struct *comlist, long** xphys) {

	double comms_time = 0; // temp
	clock_t start, end;

	long i;
	int dir, k;

	long maxindex = p->sites + p->halos;

	int* siterank = malloc(maxindex * sizeof(siterank)); // rank of the node where site i resides in
	for (i=0; i<maxindex; i++) {
		siterank[i] = coordsToRank(*p, xphys[i]);
	}

	/* N-dimensional hypercube has N^3 - 1 neighboring hypercubes on a lattice (imagine a 3x3x3x... grid)
	* However, in some directions the neighbor may be the node itself, due to periodicity (small lattice/few nodes).
	* Here we assume that these have been dealt with already in layout.c, (references to "self halos" changed to point to real sites instead),
	* so we don't have to worry about including the these in comlists.
	*/
	
	int MAXNODES = p->dim * p->dim * p->dim - 1;
	int* neighbornodes = malloc(MAXNODES * sizeof(neighbornodes));
	for (k=0; k<MAXNODES; k++) {
		neighbornodes[k] = -1;
	}

	// allocate enough memory to host all needed receive structs, realloc later
	// note that comlist->from_node itself is a pointer to sendrecv_struct
	comlist->from_node = malloc(MAXNODES * sizeof(*(comlist->from_node)));

	// message size and tag
	long size = maxindex * p->dim; // check this
	int tag = 0;

	// Where to receive from?
	int nn = 0; // how many distinct neighbors have we found (nn = neighboring node)
	int newnode;
	for (i=p->sites; i<maxindex; i++) {
		// skip self halos
		if (siterank[i] == p->rank) {
			continue;
		}
		
		newnode = 1;
		// new neighbor node?
		for (k=0; k<MAXNODES; k++) {
			if (siterank[i] == neighbornodes[k]) {
				// node already in neighbors, so add the site in its comlists
				newnode = 0;
				comlist->from_node[k].even += p->parity[i];
				comlist->from_node[k].odd += p->parity[i];
				comlist->from_node[k].sitelist[ comlist->from_node[k].sites ] = i;
				comlist->from_node[k].sites++;
			}
		}
		if (newnode == 1) {
			// found new neighbor node, so initialize sendrecv_struct!
			comlist->from_node[nn].node = siterank[i];
			neighbornodes[nn] = siterank[i];
			comlist->from_node[nn].sites = 1;
			comlist->from_node[nn].even = p->parity[i];
			comlist->from_node[nn].odd = p->parity[i];
			
			// list of sites to receive from this node. Realloc later.
			comlist->from_node[nn].sitelist = malloc(p->halos * sizeof(comlist->from_node[nn].sitelist));
			comlist->from_node[nn].sitelist[0] = i;
			nn++;
		}
	}

	if (nn > MAXNODES) {
		printf("Node %d: Failed to count neighbor nodes in layout.c! Found %d, but theoretical max is %d.\n", p->rank, nn, MAXNODES);
		die(11);
	}

	// realloc receive structs
	comlist->neighbors = nn;
	comlist->from_node = realloc(comlist->from_node, nn * sizeof(*(comlist->from_node)));
	for (k=0; k<nn; k++) {
		comlist->from_node[k].sitelist = realloc(comlist->from_node[k].sitelist, 
																	comlist->from_node[k].sites * sizeof(comlist->from_node[k].sitelist));
	}
	
	
	// Where to send? These should be the same nodes where we receive from.
	// To fill in sitelist in sendrecv_structs, we need xphys on their node, so send this with MPI.
	comlist->to_node = malloc(nn * sizeof(*(comlist->to_node)) );

	// receive buffer:
	long** xphys_nn[nn];
	// send buffer (can send same buffer to only one node at a time, so need an array):
	long** buf[nn];
	MPI_Request req[nn];
	MPI_Request req2[nn];
	MPI_Barrier(MPI_COMM_WORLD);

	// shorthands
	sendrecv_struct* recv;
	sendrecv_struct* send;

	for (k=0; k<nn; k++) {
		recv = &(comlist->from_node[k]);
		
		// alloc list of sites in sendrecv_struct, realloc later (best to do this before MPI comms).
		comlist->to_node[k].sitelist = malloc(p->sites * sizeof(comlist->to_node[k].sitelist));
		
		buf[k] = alloc_latticetable(p->dim, maxindex); 
		// now buf[k] points to a contiguous memory address that can hold xphys
		xphys_nn[k] = alloc_latticetable(p->dim, maxindex);
		

		memcpy(&buf[k][0][0], &xphys[0][0], maxindex * p->dim * sizeof(xphys[0][0]));
		
		if (recv->node != p->rank) {		
			printf("Node %d: Sending xphys to node %d\n", p->rank, recv->node);
			MPI_Isend(&buf[k][0][0], size, MPI_LONG, recv->node, tag, MPI_COMM_WORLD, &req[k]);		
			
			
			// request xphys from the receiving node
			printf("Node %d: Requesting xphys from node %d\n", p->rank, recv->node);
			// blocking receive here
			
			start = clock();
			
			MPI_Irecv(&xphys_nn[k][0][0], size, MPI_LONG, recv->node, tag, MPI_COMM_WORLD, &req2[k]);
			//printf("Node %d: Received xphys from node %d\n", p->rank, recv->node);
			
			end = clock();
			
		}
	}
	// MPI_Isend is non-blocking, so we move on.

	// Now receive their xphys_nn and use that to construct sendrecv_struct
	for (k=0; k<nn; k++) {

		recv = &(comlist->from_node[k]);
		send = &(comlist->to_node[k]);
		// initialize sendrecv_struct for this neighbor
		send->node = recv->node;
		send->sites = 0;
		send->even = 0;
		send->odd = 0;
		
		MPI_Wait(&req2[k], MPI_STATUS_IGNORE);

		int match;
		// find matching site in both nodes, and their indices
		for (long j=0; j<maxindex; j++) {
			// j = their index	
			for (i=0; i<p->sites; i++) {
				// i = my index
				match = 1;
				for (dir=0; dir<p->dim; dir++) {
					if (xphys[i][dir] != xphys_nn[k][j][dir]) {
						match = 0;
						break;
					}
				}
				
				if (match) {
					// Found matching site. Quick sanity check:
					if (j < p->sites) {
						// but now it looks like my physical non-halo site in the other node is also a non-halo site there!
						printf("Nodes %d and %d: Error in make_comlists()! Check halo indexing?\n", p->rank, recv->node);
						die(14);
					}
					// add to sendrecv_struct 
					send->even += p->parity[i];
					send->odd += p->parity[i];
					
					/* store site index i in my node to my to_node. 
					* Note that since we loop over their xphys (j loop) in the SAME order as 
					* when we constructed sitelist for THEIR from_node, we automatically get the 
					* sends in the same order as their receives!
					* (apart from parity ordering, which is done later and does not reorder even/odd sites among themselves)
					*/ 
					send->sitelist[send->sites] = i;
					send->sites++;
				}

			} // end i

		} // end j

	} // end k
	
	comms_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Node %d: Receiving took %lf seconds.\n", p->rank, comms_time);
	
	// wait for sends to finish and free buffer (MPI_Wait also frees request)
	for (k=0; k<nn; k++) {
		start = clock();
		
		MPI_Wait(&req[k], MPI_STATUS_IGNORE);
		
		end = clock();
		
		free_latticetable(maxindex, buf[k]);
		free_latticetable(maxindex, xphys_nn[k]);
		
		// Finally, realloc and rearrange sitelist so that sites with parity = 0 come first, then parity = 1.
		comlist->to_node[k].sitelist = realloc(comlist->to_node[k].sitelist, 
																	comlist->to_node[k].sites * sizeof(comlist->to_node[k].sitelist));
		reorder_sitelist(p, &comlist->from_node[k]);
		reorder_sitelist(p, &comlist->to_node[k]);
	}

	comms_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Node %d: Waiting took %lf seconds.\n", p->rank, comms_time);
	
	
	
	free(neighbornodes);
	free(siterank);

}

/* Quick routine for ordering sitelist in sendrecv_struct
* so that sites with parity = 0 come before parity = 1.
*/
void reorder_sitelist(params* p, sendrecv_struct* sr) {
	
	long* temp = malloc(sr->sites * sizeof(sr->sitelist));
	memcpy(temp, sr->sitelist, sr->sites * sizeof(sr->sitelist));
	
	long site;
	long even = 0;
	long odd = 0;	
	for (long i=0; i<sr->sites; i++) {
		site = temp[i];
		if (p->parity[site] == 0) {
			sr->sitelist[odd] = site;
			odd++;
		} else {
			sr->sitelist[sr->odd + even] = site;
			even++;
		}
	}
	
	free(temp);
}
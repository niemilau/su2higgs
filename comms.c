/** @file comms.c
*
*	Routines for constructing communication lookup tables,
* and for halo communications and reducing numbers between the nodes.
*
* Performance goal: communications take max ~10% of all time (Kari has ~6% - 17% in the simulations I did)
*
* TODO optimize make_comlists. This routine is very very slow if the lattice is large AND only a few nodes are used.
* Slowness is probably because of the xphys lookup when constructing send_to structure...
* There is no problem when the volume/node ratio is good
*
*/

#include "su2.h"
#include "comms.h"


#ifdef MPI

/*
*
*/
void make_comlists(params *p, comlist_struct *comlist, long** xphys) {

	long i;
	int dir, k;

	// work arrays
	int* neighbornodes;
	int* siterank;
	long** xphys_nn;

	/* N-dimensional hypercube has N^3 - 1 neighboring hypercubes on a lattice (imagine a 3x3x3x... grid)
	* However, in some directions the neighbor may be the node itself, due to periodicity (small lattice/few nodes).
	* Here we assume that these have been dealt with already in layout.c, (references to "self halos" changed to point to real sites instead),
	* so we don't have to worry about including the these in comlists.
	*/

	int MAXNODES = p->dim * p->dim * p->dim - 1;
	neighbornodes = malloc(MAXNODES * sizeof(*neighbornodes));
	for (k=0; k<MAXNODES; k++) {
		neighbornodes[k] = -1;
	}

	long maxindex = p->sites_total;

	siterank = malloc(maxindex * sizeof(*siterank)); // rank of the node where site i resides in
	for (i=0; i<maxindex; i++) {
		siterank[i] = coordsToRank(*p, xphys[i]);
	}

	// allocate enough memory to host all needed receive structs, realloc later
	// note that comlist->recv_from itself is a pointer to sendrecv_struct
	comlist->recv_from = malloc(MAXNODES * sizeof(*(comlist->recv_from)));

	// message size and tag for sending xphys
	long size = maxindex * p->dim;
	int tag = 0;

	// Where to receive from?
	int nn = 0; // how many distinct neighbors have we found (nn = neighboring node)
	int newnode;
	for (i=p->sites; i<maxindex; i++) {

		if (siterank[i] == p->rank) {
			// shouldn't get here!
			printf("Node %d: Self halos not removed! in comms.c\n", p->rank);
			continue;
		}

		newnode = 1;
		// new neighbor node?
		for (k=0; k<MAXNODES; k++) {
			if (siterank[i] == neighbornodes[k]) {
				// node already in neighbors, so add the site in its comlists
				newnode = 0;
				if (p->parity[i] == EVEN) {
					comlist->recv_from[k].even++;
				} else {
					comlist->recv_from[k].odd++;
				}
				comlist->recv_from[k].sitelist[ comlist->recv_from[k].sites ] = i;
				comlist->recv_from[k].sites++;
			}
		}
		if (newnode == 1) {
			// found new neighbor node, so initialize sendrecv_struct!
			comlist->recv_from[nn].node = siterank[i];
			neighbornodes[nn] = siterank[i];
			comlist->recv_from[nn].sites = 1;
			if (p->parity[i] == EVEN) {
				comlist->recv_from[nn].even = 1;
				comlist->recv_from[nn].odd = 0;
			} else {
				comlist->recv_from[nn].odd = 1;
				comlist->recv_from[nn].even = 0;
			}
			// list of sites to receive from this node. Realloc later.
			comlist->recv_from[nn].sitelist = malloc(p->halos * sizeof(*(comlist->recv_from[nn].sitelist)));
			comlist->recv_from[nn].sitelist[0] = i;
			nn++;
		}
	}

	if (nn > MAXNODES) {
		printf("Node %d: Failed to count neighbor nodes in layout.c! Found %d, but theoretical max is %d.\n", p->rank, nn, MAXNODES);
		die(11);
	}

	// realloc receive structs (TODO: add error handling here)
	comlist->neighbors = nn;
	comlist->recv_from = realloc(comlist->recv_from, nn * sizeof(*(comlist->recv_from)));
	for (k=0; k<nn; k++) {
		comlist->recv_from[k].sitelist = realloc(comlist->recv_from[k].sitelist,
																	comlist->recv_from[k].sites * sizeof(*(comlist->recv_from[k].sitelist)));
	}


	// Where to send? These should be the same nodes where we receive from.
	// To fill in sitelist in sendrecv_structs, we need xphys on their node, so send this with MPI.
	comlist->send_to = malloc(nn * sizeof(*(comlist->send_to)) );

	// receive buffer:
	xphys_nn = alloc_latticetable(p->dim, maxindex);

	// send buffer (can send same buffer to only one node at a time, so need an array):
	long** buf[nn];
	MPI_Request req[nn];
	MPI_Barrier(MPI_COMM_WORLD);

	// shorthands
	sendrecv_struct* recv;
	sendrecv_struct* send;

	for (k=0; k<nn; k++) {
		recv = &(comlist->recv_from[k]);

		// alloc list of sites in sendrecv_struct, realloc later
		comlist->send_to[k].sitelist = malloc(p->halos * sizeof(*(comlist->send_to[k].sitelist)));

		buf[k] = alloc_latticetable(p->dim, maxindex);
		// now buf[k] points to a contiguous memory address that can hold xphys

		memcpy(&buf[k][0][0], &xphys[0][0], maxindex * p->dim * sizeof(xphys[0][0]));

		if (recv->node != p->rank) {
			//printf("Node %d: Sending xphys to node %d\n", p->rank, recv->node);
			MPI_Isend(&buf[k][0][0], size, MPI_LONG, recv->node, tag, MPI_COMM_WORLD, &req[k]);
			}

	}
	// MPI_Isend is non-blocking, so we move on.

	// Now receive their xphys_nn and use that to construct sendrecv_struct
	for (k=0; k<nn; k++) {

		recv = &(comlist->recv_from[k]);
		send = &(comlist->send_to[k]);
		// initialize sendrecv_struct for this neighbor
		send->node = recv->node;
		send->sites = 0;
		send->even = 0;
		send->odd = 0;

		// request xphys from the receiving node

		// blocking receive here
		MPI_Recv(&(xphys_nn[0][0]), size, MPI_LONG, recv->node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int match;
		// find matching site in both nodes, and their indices.
		for (long j=p->sites; j<maxindex; j++) {
			// j = their index. Halos start from j = p->sites
			for (i=0; i<p->sites; i++) {
				// i = my index
				match = 1;
				for (dir=0; dir<p->dim; dir++) {
					if (xphys[i][dir] != xphys_nn[j][dir]) {
						match = 0;
						break;
					}
				}

				if (match) {
					// Found matching site, so add to sendrecv_struct
					if (p->parity[i] == EVEN) {
						send->even++;
					} else {
						send->odd++;
					}
					/* store site index i in my node to my send_to.
					* Note that since we loop over their xphys (j loop) in the SAME order as
					* when we constructed sitelist for THEIR recv_from, we automatically get the
					* sends in the same order as their receives!
					* (apart from parity ordering, which does not reorder even/odd sites among themselves).
					* NB! Since the same site in our node may map onto many halo sites in their node, we may end up sending the same
					* site multiple times. For optimization, we could implement this in receive structures instead...
					*/
					send->sitelist[send->sites] = i;
					send->sites++;
				}

			} // end i

		} // end j

	} // end k


	// wait for sends to finish and free buffer (MPI_Wait also frees request)
	for (k=0; k<nn; k++) {

		MPI_Wait(&req[k], MPI_STATUS_IGNORE);

		free_latticetable(buf[k]);

		// Finally, realloc and rearrange sitelist so that sites with parity = EVEN come first, then parity = ODD.
		comlist->send_to[k].sitelist = realloc(comlist->send_to[k].sitelist,
																	comlist->send_to[k].sites * sizeof(*(comlist->send_to[k].sitelist)));
		reorder_sitelist(p, &comlist->recv_from[k]);
		reorder_sitelist(p, &comlist->send_to[k]);
	}
	reorder_comlist(p, comlist); // arranges send_to/recv_from by node ranks


	free_latticetable(xphys_nn);
	free(neighbornodes);
	free(siterank);

}


/* Update all halos for a gauge link in direction dir on sites with given parity.
* Uses nonblocking sends to first send all required data to all neighbors,
* and then blocking receive to update my own halos.
* The advantage over blocking sendrecv is that since we just send everything
* before starting receives, waiting for neighbors to send is significantly reduced.
*
* Return value is the time in seconds spent in the routine call.
*/
double update_gaugehalo(comlist_struct* comlist, char parity, double*** field, int dofs, int dir) {

	double start, end, time = 0.0;
	start = clock();

	int neighbors = comlist->neighbors;
	// first do a nonblocking send to all neighbors
	MPI_Request send_req[neighbors];

	for (int k=0; k<neighbors; k++) {
		// send buffer is allocated here, so need to free it later.
		send_gaugefield(&comlist->send_to[k], &send_req[k], parity, field, dofs, dir);
	}

	// all sends done, now receive from all neighbors one by one
	for (int k=0; k<neighbors; k++) {
		recv_gaugefield(&comlist->recv_from[k], parity, field, dofs, dir);
	}

	// finally, wait until my sends have been received and free the send buffers
	double s, e;
	s = clock();
	for (int k=0; k<neighbors; k++) {
		MPI_Wait(&send_req[k], MPI_STATUS_IGNORE);
		free(comlist->send_to[k].buf);
	}
	e = clock();
	waittime += (e - s) /CLOCKS_PER_SEC;


	end = clock();

	time += (double)(end - start) / CLOCKS_PER_SEC;
	return time;
}


/* Nonblocking send to a given neighbor for a gauge link.
* Send buffer is allocated here but not freed; freeing is performed
* in update_gaugehalo() after all receives are complete.
*/
void send_gaugefield(sendrecv_struct* send, MPI_Request* req, char parity, double*** field, int dofs, int dir) {

	int dest = send->node; // rank of the receiving node

	long send_offset, send_max;
	long send_count; // how many MPI_DOUBLEs do we send
	int tag=0;

	// sites with parity = EVEN come first
	if (parity == EVEN) {
		// update even
		send_offset = 0;
		send_count = send->even;
		send_max = send->even;
	} else {
		// update odd
		send_offset = send->even;
		send_count = send->odd;
		send_max = send->sites;
	}
	send_count *= dofs;

	// alloc buffer and copy field values there
	send->buf = malloc(send_count * sizeof(*(send->buf)));

	long j = 0, index;
	for (long i = send_offset; i<send_max; i++) {
		index = send->sitelist[i];
		for (int d=0; d<dofs; d++) {
			// this line is the only difference to send_field()
			send->buf[j + d] = field[index][dir][d];
		}
		j += dofs;
	}

	// nonblocking send
	MPI_Isend(send->buf, send_count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);

}


/* Blocking receive from a given neighbor for a gauge link.
* Receive buffer is both allocated and freed here.
*/
void recv_gaugefield(sendrecv_struct* recv, char parity, double*** field, int dofs, int dir) {

	int source = recv->node; // rank of the sending node

	long recv_offset, recv_max;
	long recv_count; // how many MPI_DOUBLEs do we receive
	int tag=0;

	// sites with parity = EVEN come first
	if (parity == EVEN) {
		// update even
		recv_offset = 0;
		recv_count = recv->even;
		recv_max = recv->even;
	} else {
		// update odd
		recv_offset = recv->even;
		recv_count = recv->odd;
		recv_max = recv->sites;
	}
	recv_count *= dofs;

	// alloc buffer and receive
	recv->buf = malloc(recv_count * sizeof(*(recv->buf)));
	MPI_Recv(recv->buf, recv_count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// copy values from buffer to the field array
	long j = 0, index;
	for (long i = recv_offset; i<recv_max; i++) {
		index = recv->sitelist[i];
		for (int d=0; d<dofs; d++) {
			// this bit is slightly different from sendrecv_field(){
			field[index][dir][d] = recv->buf[j + d];
		}
		//printf("Received site %ld, value now: %lf, %lf, %lf \n", index, field[index][dir][0], field[index][dir][1], field[index][dir][2]);
		j += dofs;
	}

	free(recv->buf);

}

/* Same as update_gaugehalo(), but for a normal field with dof components.
* TODO: generalize this to also work for gauge fields, so that separate routines are not needed (see write_field() in checkpoint.c)
*/
double update_halo(comlist_struct* comlist, char parity, double** field, int dofs) {

	double start, end, time = 0.0;
	start = clock();

	int neighbors = comlist->neighbors;
	// first do a nonblocking send to all neighbors
	MPI_Request send_req[neighbors];

	for (int k=0; k<neighbors; k++) {
		// send buffer is allocated here, so need to free it later.
		send_field(&comlist->send_to[k], &send_req[k], parity, field, dofs);
	}

	// all sends done, now receive from all neighbors one by one
	for (int k=0; k<neighbors; k++) {
		recv_field(&comlist->recv_from[k], parity, field, dofs);
	}

	// finally, wait until my sends have been received and free the send buffers
	double s, e;
	s = clock();
	for (int k=0; k<neighbors; k++) {
		MPI_Wait(&send_req[k], MPI_STATUS_IGNORE);
		free(comlist->send_to[k].buf);
	}
	e = clock();
	waittime += (e - s) /CLOCKS_PER_SEC;


	end = clock();

	time += (double)(end - start) / CLOCKS_PER_SEC;
	return time;
}


/* Same as send_gaugehalo(), but for a normal field with dof components.
*/
void send_field(sendrecv_struct* send, MPI_Request* req, char parity, double** field, int dofs) {

	int dest = send->node; // rank of the receiving node

	long send_offset, send_max;
	long send_count; // how many MPI_DOUBLEs do we send
	int tag=0;

	// sites with parity = EVEN come first
	if (parity == EVEN) {
		// update even
		send_offset = 0;
		send_count = send->even;
		send_max = send->even;
	} else {
		// update odd
		send_offset = send->even;
		send_count = send->odd;
		send_max = send->sites;
	}
	send_count *= dofs;

	// alloc buffer and copy field values there
	send->buf = malloc(send_count * sizeof(*(send->buf)));

	long j = 0, index;
	for (long i = send_offset; i<send_max; i++) {
		index = send->sitelist[i];
		for (int d=0; d<dofs; d++) {
			// this line is the only difference to send_gaugefield()
			send->buf[j + d] = field[index][d];
		}
		j += dofs;
	}

	// nonblocking send
	MPI_Isend(send->buf, send_count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req);

}


/* Same as recv_gaugefield, but for a normal field with dofs components.
*/
void recv_field(sendrecv_struct* recv, char parity, double** field, int dofs) {

	int source = recv->node; // rank of the sending node

	long recv_offset, recv_max;
	long recv_count; // how many MPI_DOUBLEs do we receive
	int tag=0;

	// sites with parity = EVEN come first
	if (parity == EVEN) {
		// update even
		recv_offset = 0;
		recv_count = recv->even;
		recv_max = recv->even;
	} else {
		// update odd
		recv_offset = recv->even;
		recv_count = recv->odd;
		recv_max = recv->sites;
	}
	recv_count *= dofs;

	// alloc buffer and receive
	recv->buf = malloc(recv_count * sizeof(*(recv->buf)));
	MPI_Recv(recv->buf, recv_count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// copy values from buffer to the field array
	long j = 0, index;
	for (long i = recv_offset; i<recv_max; i++) {
		index = recv->sitelist[i];
		for (int d=0; d<dofs; d++) {
			// this line is the only difference to recv_gaugefield()
			field[index][d] = recv->buf[j + d];
		}
		//printf("Received site %ld, value now: %lf, %lf, %lf \n", index, field[index][0], field[index][1], field[index][2]);
		j += dofs;
	}

	free(recv->buf);

}


/* Quick routine for ordering sitelist in sendrecv_struct
* so that sites with parity = 0 come before parity = 1.
* NB! This routine fails if some lattice side L is an odd number,
* in which case parity is not even meaningful really.
* Can result in memory errors in such cases!!
*/
void reorder_sitelist(params* p, sendrecv_struct* sr) {

	long* temp = malloc(sr->sites * sizeof(*(sr->sitelist)));
	memcpy(temp, sr->sitelist, sr->sites * sizeof(*(sr->sitelist)));

	long site;
	long even = 0;
	long odd = 0;

	for (long i=0; i<sr->sites; i++) {
		site = temp[i];
		if (p->parity[site] == EVEN) {
			sr->sitelist[even] = site;
			even++;
		} else {
			sr->sitelist[sr->even + odd] = site;
			odd++;
		}
	}

	free(temp);
}

/* Routine for ordering send/recv structs in comlist
* by their MPI rank.
*/
void reorder_comlist(params* p, comlist_struct* comlist) {

	int rank;

	sendrecv_struct* send;
	sendrecv_struct* recv;

	send = malloc(comlist->neighbors * sizeof(*send));
	// copy contents to temp struct. this also copies pointer to address of sitelist, no need to alloc/free it
	memcpy(send, comlist->send_to, comlist->neighbors * sizeof(*send));

	int n = 0;
	for (rank=0; rank<p->size; rank++) {
		for (int k=0; k<comlist->neighbors; k++) {
			if (send[k].node == rank) {
				comlist->send_to[n].sitelist = send[k].sitelist;
				comlist->send_to[n].node = rank;
				comlist->send_to[n].even = send[k].even;
				comlist->send_to[n].odd = send[k].odd;
				comlist->send_to[n].sites = send[k].sites;
				n++;
				break;
			}
		}
	}

	free(send);

	// same for recv_from
	recv = malloc(comlist->neighbors * sizeof(*recv));
	memcpy(recv, comlist->recv_from, comlist->neighbors * sizeof(*recv));

	n = 0;
	for (rank=0; rank<p->size; rank++) {
		for (int k=0; k<comlist->neighbors; k++) {
			if (recv[k].node == rank) {
				comlist->recv_from[n].sitelist = recv[k].sitelist;
				comlist->recv_from[n].node = rank;
				comlist->recv_from[n].even = recv[k].even;
				comlist->recv_from[n].odd = recv[k].odd;
				comlist->recv_from[n].sites = recv[k].sites;
				n++;
				break;
			}
		}
	}

	free(recv);

}


/* Test halo send and receive (specifically update_gaugehalo())
* This routine first checks that individual parity EVEN or ODD updates
* do not change sites of the opposite parity,
* and that all halos are eventually updated while real sites are not.
*
* See test_comms_individual() for another test routine.
*/
void test_comms(params p, comlist_struct comlist) {

	long i;
	int dir, dof;

	// field for testing purposes
	int maxdof = 3;
	double*** field = make_gaugefield(p.sites_total, p.dim, maxdof);
	// give some values that are easily tracked (0.0 for halos)
	for (i=0; i<p.sites_total; i++) {
		for (dir=0; dir<p.dim; dir++) {
			for (dof=0; dof<maxdof; dof++) {
				field[i][dir][dof] = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
			}
		}
	}

	// update EVEN links in some direction
	int testdir = p.dim - 1;
	update_gaugehalo(&comlist, EVEN, field, maxdof, testdir);

	// check that all EVEN halos changed in testdir, and that nothing else changed.
	// since we work with doubles, need to be careful when comparing values.
	// precision 10^(-4) should be sufficient for small maxdof
	for (i=0; i<p.sites_total; i++) {
		for (dir=0; dir<p.dim; dir++) {
			if (i >= p.sites && p.parity[i] == EVEN && dir == testdir) {
				// should have changed
				for (dof=0; dof<maxdof; dof++) {
					double oldval = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld with EVEN parity was not updated\n", p.rank, i);
						die(-120);
					}
				}
			} else {
				// should have no change
				for (dof=0; dof<maxdof; dof++) {
					double oldval = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) > 0.0001 ) {
						printf("Node %d: Error in test_comms! Site %ld was updated in EVEN sweep, when it should not have been\n", p.rank, i);
						die(-121);
					}
				}
			}
		}
	}

	// reset the field and repeat for ODD sites
	for (i=0; i<p.sites_total; i++) {
		for (dir=0; dir<p.dim; dir++) {
			for (dof=0; dof<maxdof; dof++) {
				field[i][dir][dof] = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
			}
		}
	}

	update_gaugehalo(&comlist, ODD, field, maxdof, testdir);

	for (i=0; i<p.sites_total; i++) {
		for (dir=0; dir<p.dim; dir++) {
			if (i >= p.sites && p.parity[i] == ODD && dir == testdir) {
				// should have changed
				for (dof=0; dof<maxdof; dof++) {
					double oldval = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld with ODD parity was not updated \n", p.rank, i);
						die(-122);
					}
				}
			} else {
				// should have no change
				for (dof=0; dof<maxdof; dof++) {
					double oldval = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) > 0.0001 ) {
						printf("Node %d: Error in test_comms! Site %ld was updated in ODD sweep, when it should not have been \n", p.rank, i);
						die(-123);
					}
				}
			}
		}
	}

	// update EVEN again and check that all halos have changed in testdir
	update_gaugehalo(&comlist, EVEN, field, maxdof, testdir);

	for (i=p.sites; i<p.sites_total; i++) {
		for (dir=0; dir<p.dim; dir++) {
			if (dir == testdir) {
				for (dof=0; dof<maxdof; dof++) {
					double oldval = p.rank * p.sites_total * p.dim + i * p.dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld was not updated in either EVEN nor ODD sweep \n", p.rank, i);
						die(-124);
					}
				}
			}
		}
	}

	// tests done, can free test field
	free_gaugefield(p.sites_total, field);
}


/*
* Perform a complete check on the comlist. By complete we mean that
* for each lattice site, we check that the field we send from our site actually ends up in their correct site.
* Such check requires either:
* 	1) Exhaustive communications to know the physical coords of my real site in their halo
*			and values of the sent/received field at each site
* or
*		2) A way to assign an unique value to the field at each site, based on the physical coordinates of the site.
*
* We perform 2), using a recursive generalization of the Cantor pairing function p(x,y) = 0.5*(x+y)*(x+y+1) + y.
* Speficially, if a site has coords (x1, x2, x3, ...), we calculate y1 = p(x1,x2),
* then y2 = p(x3, y) etc. This gives an unique number for each site that we assign to the field,
* plus a decimal number for each component. It is then easy to use xphys to predict what the received value should be.
*
* Uses update_field() instead of update_gaugefield() for simplicity.
*/
void test_comms_individual(params p, comlist_struct comlist, long** xphys) {
	long i, x, y;
	int dir, dof;

	// field for testing purposes
	int maxdof = 4;
	double** field = make_field(p.sites_total, maxdof);
	// give values according to the Cantor pairing function, but set halos to 0
	for (i=0; i<p.sites_total; i++) {

		y = xphys[i][0];
		for (dir=1; dir<p.dim; dir++) {
			x = xphys[i][dir];
			y = y + 0.5 * (x + y + 1) * (x + y);
		}

		for (dof=0; dof<maxdof; dof++) {
			if (i >= p.sites) {
				field[i][dof] = 0;
			} else {
				// Obtain base number from Cantor, set different decimals for different components
				field[i][dof] = y + (double) dof / maxdof;
			}
		}
	}

	// now receive my halo fields from neighbors
	update_halo(&comlist, EVEN, field, maxdof);
	update_halo(&comlist, ODD, field, maxdof);

	/* use xphys to calculate what the field value should be in my halos,
	* according to the Cantor pairing. This is a strong check on sitelist in sendrecv_structs,
	* because the (x,y,z,...) coords of my halo site should match those of the real site
	* in the node where we received the field value from.
	*/
	for (i=p.sites; i<p.sites_total; i++) {

		y = xphys[i][0];
		for (dir=1; dir<p.dim; dir++) {
			x = xphys[i][dir];
			y = y + 0.5 * (x + y + 1) * (x + y);
		}

		for (dof=0; dof<maxdof; dof++) {
			long val = y + (double) dof / maxdof;
			if (abs(val - field[i][dof]) > 0.001) {
				// predicted value does not match what was sent...
				// print error, but don't die
				printf("Node %d: Error in test_comms! Halo site %ld did not receive correct value from update_halo() \n", p.rank, i);
			}
		}

	}

	//printf("p.sites_total = %ld \n", p.sites_total);
	//printf("Node %d: site i = 1171, x = %ld, y = %ld, z = %ld \n", p.rank, xphys[1171][0], xphys[1171][1], xphys[1171][2] );

	// done, free test field
	free_field(field);
}


/* Add together doubles from each node
* and collect the result in master node (rank = 0).
* res is the variable that is to be collected and added.
*/
double reduce_sum(double res) {

  double total = 0.0;
  MPI_Reduce(&res, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return total;
}

// Same as reduce_sum(), but the result is distributed to all nodes instead of just master node
double allreduce(double res) {
	double total = 0.0;
	MPI_Allreduce(&res, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return total;
}

// Broadcast integer from root node (rank = 0) to all other nodes.
void bcast_int(int *res) {
  MPI_Bcast(res, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
}

// Broadcast long integer from root node (rank = 0) to all other nodes.
void bcast_long(long *res) {
  MPI_Bcast(res, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}

// Broadcast double from root node (rank = 0) to all other nodes.
void bcast_double(double *res) {
  MPI_Bcast(res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Barrier MPI_COMM_WORLD
void barrier() {
	MPI_Barrier(MPI_COMM_WORLD);
}

#else // No MPI, dummy routines. comlist is not even needed in this case

void barrier() {}

double update_gaugehalo(comlist_struct* comlist, char parity, double*** field, int dofs, int dir) {
	return 0;
}

double update_halo(comlist_struct* comlist, char parity, double** field, int dofs) {
	return 0;
}

double reduce_sum(double res) {
  return res;
}

double allreduce(double res) {
	return res;
}

void bcast_int(int *res) {
	return;
}

void bcast_long(long *res) {
	return;
}

#endif

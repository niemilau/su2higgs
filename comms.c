/** @file comms.c
*
*	Routines for constructing communication lookup tables,
* and for halo communications and reducing numbers between the nodes.
*
* TODO optimize make_comlists. This routine is quite slow if the lattice is large AND only a few nodes are used.
* Slowness is probably because of the l->coords lookup when constructing send_to structure.
* There is no problem when the volume/node ratio is good
*
*/

#include "su2.h"
#include "comms.h"
#include "generic/stddefs.h"


/* Formatted printing to stdout from root node only. */
void printf0(char *msg, ...) {
  va_list fmtargs;
  va_start(fmtargs, msg);

  if(!myRank) {
    vprintf(msg, fmtargs);
		fflush(stdout);
	}
}

/* Add a site to the sitelist in comlist.send or comlist.recv.
* First checks if the given node already has an associated sendrecv structure;
* if not, allocs sitelist (of size init_max) in a new sendrecv struct and initializes it.
* sendrecv = SEND -> add to send struct,
* sendrecv = RECV -> add to recv struct.
* Return 1 if a new node was found, 0 otherwise. */
int addto_comlist(comlist_struct* comlist, int rank, long i, int sendrecv, char evenodd, long init_max) {

  sendrecv_struct* sr;
  int* sr_nodes;
  if (sendrecv == SEND) {
    sr_nodes = &comlist->sends;
  } else if (sendrecv == RECV) {
    sr_nodes = &comlist->recvs;
  } else {
    printf("Error in comlist.c! invalid input\n");
  }

  // check if the node already has an associated sendrecv struct
  for (int n=0; n < *sr_nodes; n++) {

    if (sendrecv == SEND) {
      sr = &comlist->send_to[n];
    } else {
      sr = &comlist->recv_from[n];
    }

    if (sr->node == rank) {
      // found existing sendrecv struct

      sr->sitelist[ sr->sites ] = i;

      sr->sites++;
      if (evenodd == EVEN) {
        sr->even++;
      } else {
        sr->odd++;
      }

      return 0;
    }
  }

  // new node to send/recv, so alloc and initialize the sendrecv struct
  int k = *sr_nodes;
  *sr_nodes = *sr_nodes + 1;

  if (sendrecv == SEND) {
    if (k == 0) {
      comlist->send_to = malloc(1 * sizeof(*(comlist->send_to))); // alloc just one first
    } else {
      // need more, so realloc
      comlist->send_to = realloc(comlist->send_to, (k+1)* sizeof(*(comlist->send_to)));
    }
    sr = &comlist->send_to[k];
  } else {
    if (k == 0) {
      comlist->recv_from = malloc(1 * sizeof(*(comlist->recv_from)));
    } else {
      comlist->recv_from = realloc(comlist->recv_from, (k+1)* sizeof(*(comlist->recv_from)));
    }
    sr = &comlist->recv_from[k];
  }

  // then the sitelist
  sr->sitelist = malloc(init_max * sizeof(*sr->sitelist));

	if (sr->sitelist == NULL) {
		printf("malloc error in comms.c!!\n");
		die(-111);
	}

  sr->sitelist[0] = i;
  sr->node = rank;
  sr->sites = 1;
  if (evenodd == EVEN) {
    sr->even = 1;
    sr->odd = 0;
  } else {
    sr->even = 0;
    sr->odd = 1;
  }

	return 1;
}


#ifdef MPI

/* Huge routine for preparing the comlists. l->coords is assumed to be the table of (x,y,z,...)
* coordinates on the full lattice, calculated in layout().
*/
void make_comlists(lattice *l, comlist_struct *comlist) {

	long i;
	int dir, k;

	// work arrays
	int* whichnode;
	long** coords_nn;

	/* N-dimensional hypercube has N^3 - 1 neighboring hypercubes on a lattice (imagine a 3x3x3x... grid)
	* However, in some directions the neighbor may be the node itself, due to periodicity (small lattice/few nodes).
	* Here we assume that these have been dealt with already in layout.c,
	* (references to "self halos" changed to point to real sites instead),
	* so we don't have to worry about including the these in comlists. */
	int MAXNODES = l->dim * l->dim * l->dim - 1;
  // if using larger halos, this may be more complicated (below)


	long maxindex = l->sites_total;

	whichnode = malloc(maxindex * sizeof(*whichnode)); // rank of the node where site i resides in
	for (i=0; i<maxindex; i++) {
		whichnode[i] = coordsToRank(l, l->coords[i]);
	}

	// allocating of send/recv structs is done by addto_comlist, just initialize here
  comlist->sends = 0; comlist->recvs = 0;

	// message size and tag for sending l->coords
	long size = maxindex * l->dim;
	int tag = 0;

	/* Where to receive from? Loop over all halos here and figure out the node
	* where the halo site lives in, and store it in comlist->recv_from[node].sitelist*/
	int nn = 0; // how many distinct neighbors have we found (nn = neighboring node)
	int newnode;
	for (i=l->sites; i<maxindex; i++) {

		if (whichnode[i] == l->rank) {
			// shouldn't get here!
			printf("WARNING in node %d: Self halos not removed! in comms.c\n", l->rank);
			continue;
		}

		nn += addto_comlist(comlist, whichnode[i], i, RECV, l->parity[i], l->halos);
	}
	// receives done and comlist->recvs up to date (together with nn)

	if (nn > MAXNODES && HALOWIDTH == 1) {
		printf("Node %d: Failed to count neighbor nodes in layout.c! Found %d, but theoretical max is %d.\n", l->rank, nn, MAXNODES);
		die(11);
	}

	// Where to send? These should be the same nodes where we receive from.
	// To fill in sitelist in sendrecv_structs, we need l->coords on their node, so send this with MPI.

	// receive buffer:
	coords_nn = alloc_latticetable(l->dim, maxindex);

	// send buffer (can send same buffer to only one node at a time, so need an array):
	long** buf[nn];
	MPI_Request req[nn];
	barrier(l->comm);

	// shorthands
	sendrecv_struct* recv;
	sendrecv_struct* send;

	for (k=0; k<nn; k++) {
		recv = &(comlist->recv_from[k]);

		buf[k] = alloc_latticetable(l->dim, maxindex);
		// now buf[k] points to a contiguous memory address that can hold l->coords

		memcpy(&buf[k][0][0], &l->coords[0][0], maxindex * l->dim * sizeof(l->coords[0][0]));

		if (recv->node != l->rank) {
			//printf("Node %d: Sending l->coords to node %d\n", l->rank, recv->node);
			MPI_Isend(&buf[k][0][0], size, MPI_LONG, recv->node, tag, l->comm, &req[k]);
		}

	}
	// MPI_Isend is non-blocking, so we move on.

	// Now receive their coords_nn and use that to construct sendrecv_struct
	for (k=0; k<nn; k++) {

		recv = &(comlist->recv_from[k]);

		// request l->coords from the receiving node
		MPI_Recv(&(coords_nn[0][0]), size, MPI_LONG, recv->node, tag, l->comm, MPI_STATUS_IGNORE);

		int match;
		// find matching site in both nodes, and their indices.
		for (long j=l->sites; j<maxindex; j++) {
			// j = their index. Halos start from j = l->sites
			for (i=0; i<l->sites; i++) {
				// i = my index
				match = 1;
				for (dir=0; dir<l->dim; dir++) {
					if (l->coords[i][dir] != coords_nn[j][dir]) {
						match = 0;
						break;
					}
				}

				if (match) {
					// Found matching site, so add to my send_to
					addto_comlist(comlist, recv->node, i, SEND, l->parity[i], l->sites_total);
					/* need l->sites_total here as the maximum amount of sites to be sent,
					* instead of l->sites, because some of my sites may map onto multiple
					* halo sites in the receiving node. To be absolutely sure, we should
					* actually use THEIR l->halos, but here I assume that all nodes have the
					* same amount of sites. */

					/* Note that since we loop over their l->coords (j loop) in the SAME order as
					* when we constructed sitelist for THEIR recv_from, we automatically get the
					* sends in the same order as their receives!
					* (apart from parity ordering, which does not reorder even/odd sites among themselves). */
				}

			} // end i

		} // end j

	} // end k


	// wait for sends to finish and free buffer (MPI_Wait also frees request)
	for (k=0; k<nn; k++) {

		MPI_Wait(&req[k], MPI_STATUS_IGNORE);

		free_latticetable(buf[k]);

		// Finally, realloc and rearrange sitelist so that sites with parity = EVEN come first, then parity = ODD.
		realloc_comlist(comlist, SEND);
		realloc_comlist(comlist, RECV);
		reorder_sitelist(l, &comlist->recv_from[k]);
		reorder_sitelist(l, &comlist->send_to[k]);
	}
	reorder_comlist(l, comlist); // arranges send_to/recv_from by node ranks


	free_latticetable(coords_nn);
	free(whichnode);

}


/* Update all halos for a gauge link in direction dir on sites with given parity.
* Uses nonblocking sends to first send all required data to all neighbors,
* and then blocking receive to update my own halos.
* The advantage over blocking sendrecv is that since we just send everything
* before starting receives, waiting for neighbors to send is significantly reduced. */
void update_gaugehalo(lattice* l, char parity, double*** field, int dofs, int dir) {

	double start, end, time = 0.0;
	start = clock();
	comlist_struct* comlist = &l->comlist;

	int neighbors = comlist->sends;
	if (neighbors <= 0) {
		return;
	}

	// first do a nonblocking send to all neighbors
	MPI_Request send_req[neighbors];

	for (int k=0; k<neighbors; k++) {
		// send buffer is allocated here, so need to free it later.
		send_gaugefield(&comlist->send_to[k], l->comm, &send_req[k], parity, field, dofs, dir);
	}

	// all sends done, now receive from all neighbors one by one
	for (int k=0; k<neighbors; k++) {
		recv_gaugefield(&comlist->recv_from[k], l->comm, parity, field, dofs, dir);
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
	Global_comms_time += time;

}


/* Nonblocking send to a given neighbor for a gauge link.
* Send buffer is allocated here but not freed; freeing is performed
* in update_gaugehalo() after all receives are complete.
*/
void send_gaugefield(sendrecv_struct* send, MPI_Comm comm, MPI_Request* req, char parity,
					double*** field, int dofs, int dir) {

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
	} else if (parity == ODD) {
		// update odd
		send_offset = send->even;
		send_count = send->odd;
		send_max = send->sites;
	} else if (parity == EVENODD) {
    // send both EVEN and ODD
    send_offset = 0;
    send_count = send->sites;
    send_max = send->sites;
  } else {
    printf("invalid parity in send_gaugefield!! (comms.c)\n");
    return;
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
	MPI_Isend(send->buf, send_count, MPI_DOUBLE, dest, tag, comm, req);

}


/* Blocking receive from a given neighbor for a gauge link.
* Receive buffer is both allocated and freed here. */
void recv_gaugefield(sendrecv_struct* recv, MPI_Comm comm, char parity, double*** field, int dofs, int dir) {

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
	} else if (parity == ODD) {
		// update odd
		recv_offset = recv->even;
		recv_count = recv->odd;
		recv_max = recv->sites;
	} else if (parity == EVENODD) {
    // send both EVEN and ODD
    recv_offset = 0;
    recv_count = recv->sites;
    recv_max = recv->sites;
  } else {
    printf("invalid parity in recv_gaugefield!! (comms.c)\n");
    return;
  }
	recv_count *= dofs;

	// alloc buffer and receive
	recv->buf = malloc(recv_count * sizeof(*(recv->buf)));
	MPI_Recv(recv->buf, recv_count, MPI_DOUBLE, source, tag, comm, MPI_STATUS_IGNORE);

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
void update_halo(lattice* l, char parity, double** field, int dofs) {

	double start, end, time = 0.0;
	start = clock();
	comlist_struct* comlist = &l->comlist;

	int neighbors = comlist->sends;
	if (neighbors <= 0) {
		return;
	}

	// first do a nonblocking send to all neighbors
	MPI_Request send_req[neighbors];

	for (int k=0; k<neighbors; k++) {
		// send buffer is allocated here, so need to free it later.
		send_field(&comlist->send_to[k], l->comm, &send_req[k], parity, field, dofs);
	}

	// all sends done, now receive from all neighbors one by one
	for (int k=0; k<neighbors; k++) {
		recv_field(&comlist->recv_from[k], l->comm, parity, field, dofs);
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
	Global_comms_time += time;
}


/* Same as send_gaugehalo(), but for a normal field with dof components. */
void send_field(sendrecv_struct* send, MPI_Comm comm, MPI_Request* req, char parity, double** field, int dofs) {

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
	} else if (parity == ODD) {
		// update odd
		send_offset = send->even;
		send_count = send->odd;
		send_max = send->sites;
	} else if (parity == EVENODD) {
    // send both EVEN and ODD
    send_offset = 0;
    send_count = send->sites;
    send_max = send->sites;
  } else {
    printf("invalid parity in send_field!! (comms.c)\n");
    return;
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
	MPI_Isend(send->buf, send_count, MPI_DOUBLE, dest, tag, comm, req);
}


/* Same as recv_gaugefield, but for a normal field with dofs components. */
void recv_field(sendrecv_struct* recv, MPI_Comm comm, char parity, double** field, int dofs) {

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
	} else if (parity == ODD) {
		// update odd
		recv_offset = recv->even;
		recv_count = recv->odd;
		recv_max = recv->sites;
	} else if (parity == EVENODD) {
    // send both EVEN and ODD
    recv_offset = 0;
    recv_count = recv->sites;
    recv_max = recv->sites;
  } else {
    printf("invalid parity in recv_field!! (comms.c)\n");
    return;
  }
	recv_count *= dofs;

	// alloc buffer and receive
	recv->buf = malloc(recv_count * sizeof(*(recv->buf)));
	MPI_Recv(recv->buf, recv_count, MPI_DOUBLE, source, tag, comm, MPI_STATUS_IGNORE);

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
* so that sites with parity = EVEN come before parity = ODD. */
void reorder_sitelist(lattice* l, sendrecv_struct* sr) {

	long* temp = malloc(sr->sites * sizeof(*(sr->sitelist)));
	memcpy(temp, sr->sitelist, sr->sites * sizeof(*(sr->sitelist)));

	long site;
	long even = 0;
	long odd = 0;

	for (long i=0; i<sr->sites; i++) {
		site = temp[i];
		if (l->parity[site] == EVEN) {
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
* by their MPI rank. */
void reorder_comlist(lattice* l, comlist_struct* comlist) {

	int rank;

	sendrecv_struct* send;
	sendrecv_struct* recv;

  if (comlist->sends > 0) {
	  send = malloc(comlist->sends * sizeof(*send));
  	// copy contents to temp struct. this also copies pointer to address of sitelist, no need to alloc/free it
  	memcpy(send, comlist->send_to, comlist->sends * sizeof(*send));

  	int n = 0;
  	for (rank=0; rank<l->size; rank++) {
  		for (int k=0; k<comlist->sends; k++) {
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
  }

	// same for recv_from
  if (comlist->recvs > 0) {
  	recv = malloc(comlist->recvs * sizeof(*recv));
  	memcpy(recv, comlist->recv_from, comlist->recvs * sizeof(*recv));

  	int n = 0;
  	for (rank=0; rank<l->size; rank++) {
  		for (int k=0; k<comlist->recvs; k++) {
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

}


/* Test halo send and receive (specifically update_gaugehalo())
* This routine first checks that individual parity EVEN or ODD updates
* do not change sites of the opposite parity,
* and that all halos are eventually updated while real sites are not.
*
* See test_comms_individual() for another test routine. */
void test_comms(lattice* l) {

	long i;
	int dir, dof;

	// field for testing purposes
	int maxdof = 3;
	double*** field = make_gaugefield(l->sites_total, l->dim, maxdof);
	// give some values that are easily tracked (0.0 for halos)
	for (i=0; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			for (dof=0; dof<maxdof; dof++) {
				field[i][dir][dof] = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
			}
		}
	}

	// update EVEN links in some direction
	int testdir = l->dim - 1;
	update_gaugehalo(l, EVEN, field, maxdof, testdir);

	// check that all EVEN halos changed in testdir, and that nothing else changed.
	// since we work with doubles, need to be careful when comparing values.
	// precision 10^(-4) should be sufficient for small maxdof
	for (i=0; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			if (i >= l->sites && l->parity[i] == EVEN && dir == testdir) {
				// should have changed
				for (dof=0; dof<maxdof; dof++) {
					double oldval = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld with EVEN parity was not updated\n", l->rank, i);
						die(-120);
					}
				}
			} else {
				// should have no change
				for (dof=0; dof<maxdof; dof++) {
					double oldval = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) > 0.0001 ) {
						printf("Node %d: Error in test_comms! Site %ld was updated in EVEN sweep, when it should not have been\n", l->rank, i);
						die(-121);
					}
				}
			}
		}
	}

	// reset the field and repeat for ODD sites
	for (i=0; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			for (dof=0; dof<maxdof; dof++) {
				field[i][dir][dof] = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
			}
		}
	}

	update_gaugehalo(l, ODD, field, maxdof, testdir);

	for (i=0; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			if (i >= l->sites && l->parity[i] == ODD && dir == testdir) {
				// should have changed
				for (dof=0; dof<maxdof; dof++) {
					double oldval = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld with ODD parity was not updated \n", l->rank, i);
						die(-122);
					}
				}
			} else {
				// should have no change
				for (dof=0; dof<maxdof; dof++) {
					double oldval = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) > 0.0001 ) {
						printf("Node %d: Error in test_comms! Site %ld was updated in ODD sweep, when it should not have been \n", l->rank, i);
						die(-123);
					}
				}
			}
		}
	}

	// update EVEN again and check that all halos have changed in testdir
	update_gaugehalo(l, EVEN, field, maxdof, testdir);

	for (i=l->sites; i<l->sites_total; i++) {
		for (dir=0; dir<l->dim; dir++) {
			if (dir == testdir) {
				for (dof=0; dof<maxdof; dof++) {
					double oldval = l->rank * l->sites_total * l->dim + i * l->dim + dir + (double) dof / maxdof;
					if (fabs(field[i][dir][dof] - oldval) < 0.0001 ) {
						printf("Node %d: Error in test_comms! Halo site %ld was not updated in either EVEN nor ODD sweep \n", l->rank, i);
						die(-124);
					}
				}
			}
		}
	}

	// tests done, can free test field
	free_gaugefield(l->sites_total, field);
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
* I perform 2), using a recursive generalization of the Cantor pairing function p(x,y) = 0.5*(x+y)*(x+y+1) + y.
* Speficially, if a site has coords (x1, x2, x3, ...), we calculate y1 = p(x1,x2),
* then y2 = p(x3, y) etc. This gives an unique number for each site that we assign to the field,
* plus a decimal number for each component. It is then easy to use l->coords to predict what the received value should be.
*
* Uses update_field() instead of update_gaugefield() for simplicity. */
void test_comms_individual(lattice* l) {
	long i, x, y;
	int dir, dof;
	long n_err = 0;

	// field for testing purposes
	int maxdof = 4;
	double** field = make_field(l->sites_total, maxdof);
	// give values according to the Cantor pairing function, but set halos to 0
	for (i=0; i<l->sites_total; i++) {

		y = l->coords[i][0];
		for (dir=1; dir<l->dim; dir++) {
			x = l->coords[i][dir];
			y = y + 0.5 * (x + y + 1) * (x + y);
		}

		for (dof=0; dof<maxdof; dof++) {
			if (i >= l->sites) {
				field[i][dof] = 0;
			} else {
				// Obtain base number from Cantor, set different decimals for different components
				field[i][dof] = y + (double) dof / maxdof;
			}
		}
	}

	// now receive my halo fields from neighbors
	update_halo(l, EVEN, field, maxdof);
	update_halo(l, ODD, field, maxdof);

	/* use l->coords to calculate what the field value should be in my halos,
	* according to the Cantor pairing. This is a strong check on sitelist in sendrecv_structs,
	* because the (x,y,z,...) coords of my halo site should match those of the real site
	* in the node where we received the field value from.
	*/
	for (i=l->sites; i<l->sites_total; i++) {

		y = l->coords[i][0];
		for (dir=1; dir<l->dim; dir++) {
			x = l->coords[i][dir];
			y = y + 0.5 * (x + y + 1) * (x + y);
		}

		for (dof=0; dof<maxdof; dof++) {
			long val = y + (double) dof / maxdof;
			if (abs(val - field[i][dof]) > 0.001) {
				// predicted value does not match what was sent...
				// print error, but don't die
				n_err++;
				printf("Node %d: Error in test_comms! Halo site %ld did not receive correct value from update_halo() \n", l->rank, i);
			}
		}

	}

	// done, free test field
	free_field(field);
	if (n_err > 0) {
		printf("Node %d: test_comms_individual() failed with %ld errors!!\n", l->rank, n_err);
		die(-4222);
	}
}


/* Add together doubles from each node
* and collect the result in master node (rank = 0).
* res is the variable that is to be collected and added. */
double reduce_sum(double res, MPI_Comm comm) {

  double total = 0.0;
  MPI_Reduce(&res, &total, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  return total;
}

// Same as reduce_sum(), but the result is distributed to all nodes instead of just master node
double allreduce(double res, MPI_Comm comm) {
	double total = 0.0;
	MPI_Allreduce(&res, &total, 1, MPI_DOUBLE, MPI_SUM, comm);
	return total;
}

// Same as reduce_sum, but for a long type
long reduce_sum_long(long res, MPI_Comm comm) {
	long total = 0;
	MPI_Reduce(&res, &total, 1, MPI_LONG, MPI_SUM, 0, comm);
	return total;
}

// Broadcast integer from root node (rank = 0) to all other nodes.
void bcast_int(int *res, MPI_Comm comm) {
  MPI_Bcast(res, 1, MPI_INTEGER, 0, comm);
}

// Broadcast long integer from root node (rank = 0) to all other nodes.
void bcast_long(long *res, MPI_Comm comm) {
  MPI_Bcast(res, 1, MPI_LONG, 0, comm);
}

// Broadcast double from root node (rank = 0) to all other nodes.
void bcast_double(double *res, MPI_Comm comm) {
  MPI_Bcast(res, 1, MPI_DOUBLE, 0, comm);
}

// Broadcast an array of integers
void bcast_int_array(int *arr, int size, MPI_Comm comm) {
  MPI_Bcast(arr, size, MPI_INTEGER, 0, comm);
}

// Broadcast an array of long integers
void bcast_long_array(long *arr, int size, MPI_Comm comm) {
  MPI_Bcast(arr, size, MPI_LONG, 0, comm);
}

void bcast_double_array(double *arr, int size, MPI_Comm comm) {
	MPI_Bcast(arr, size, MPI_DOUBLE, 0, comm);
}

// Broadcast a string (=array of chars)
void bcast_string(char *str, int len, MPI_Comm comm) {
  MPI_Bcast(str, len, MPI_CHAR, 0, comm);
}

// barrier a given MPI commutator
void barrier(MPI_Comm comm) {
	MPI_Barrier(comm);
}

#else // No MPI, dummy routines. comlist is not even needed in this case

void barrier(MPI_Comm comm) {}

void update_gaugehalo(lattice* l, char parity, double*** field, int dofs, int dir) {
	return;
}

void update_halo(lattice* l, char parity, double** field, int dofs) {
	return;
}

double reduce_sum(double res, MPI_Comm comm) {
  return res;
}

double allreduce(double res, MPI_Comm comm) {
	return res;
}

long reduce_sum_long(long res, MPI_Comm comm) {
	return res;
}

void bcast_int(int *res, MPI_Comm comm) {
	return;
}

void bcast_long(long *res, MPI_Comm comm) {
	return;
}

void bcast_int_array(int *arr, int size, MPI_Comm comm) {
  return;
}

void bcast_long_array(long *arr, int size, MPI_Comm comm) {
  return;
}

void bcast_double_array(double *arr, int size, MPI_Comm comm) {
	return;
}

void bcast_string(char *str, int len, MPI_Comm comm) {
  return;
}

#endif

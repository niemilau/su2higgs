/** @file blocking.c
*
* Routines for creating "blocked", i.e. coarser, lattices from the original one
* and mapping fields there. The blocked lattice struct can be used just as the
* original one, particularly for calculating correlation lengths.
*
* //TODO
*/

#ifdef BLOCKING

#include "su2.h"
#include "comms.h"

/* Calculate how many blocking levels we can create for a given lattice
* Does not accept blocking that shrinks a dimension down to a point!
* The current routines also do not behave well in situations where
* the length shrinks 2 -> 1, see test_blocking(). */
int max_block_level(lattice const* l, int const* block_dir) {
  int count = 0;
  int ok = 1;
  int len[l->dim];
  memcpy(len, l->L, l->dim * sizeof(*(l->L)) );

  while (ok) {
    for (int dir=0; dir<l->dim; dir++) {
      if (block_dir[dir] && (len[dir] % 2 != 0 || len[dir] <= 2) ) {
        // cannot block this direction anymore
        ok = 0;
        return count;

      } else if (block_dir[dir]) {
        len[dir] /= 2;
      }
    }
    count++;
  }

}

/* Create a new, coarser lattice by effectively doubling the lattice
* spacing in given directions.
* block_dir[dir] = 1 if we do blocking in this direction, 0 otherwise.
*/
void block_lattice(lattice* l, lattice* b, int const* block_dir) {

  /* test that number of sites is even in each block dir,
  * and not exactly 2 because then the blocking would remove that dimension altogether */
  for (int dir=0; dir<l->dim; dir++) {
    if (block_dir[dir] && (l->L[dir] % 2 != 0 || l->L[dir] <= 2) ) {
      printf0(*l, "Can't perform any more blocking in direction %d !!\n", dir+1);
      return;
    }
  }

  // inherit dimension and MPI rank/size
  b->dim = l->dim;
  b->rank = l->rank;
  b->size = l->size;

  b->L = malloc(b->dim * sizeof(*b->L));

  // calculate number of lattice sites on the blocked lattice
  b->vol = 1;
  for (int dir=0; dir<b->dim; dir++) {
    if (block_dir[dir]) {
      b->L[dir] = l->L[dir] / 2;
    } else {
      // no blocking in this direction
      b->L[dir] = l->L[dir];
    }
    b->vol *= b->L[dir];
  }


  /* need to figure out how to layout the blocked lattice
  * on MPI nodes. Because the new lattice has less sites than the original,
  * we may have to use LESS MPI nodes to get an even slicing. The 'extra' nodes
  * will have to standby until the original lattice is used again */
  while (b->vol % b->size != 0) {
    b->size--;
  }

  if (b->rank < b->size) {
    b->standby = 0;
  } else {
    b->standby = 1;
  }

  if (!l->rank) {
    printf("Initializing blocked lattice of volume: ");
    for (int dir=0; dir<b->dim; dir++) {
      if (dir>0) printf(" x ");
      printf("%d", b->L[dir]);
    }
    printf(" (MPI size = %d)\n", b->size);
  }

  #ifdef MPI
    // create a new communicator for the blocked lattice
    int color = 1;
    if (b->standby) {
      // do not include the 'extra' nodes in the new communicator
      color = MPI_UNDEFINED;
    }

    // same ranks as in the original comm
    MPI_Comm_split(l->comm, color, l->rank, &b->comm);

  #else
    // just give some value to l->comm
  #endif

  int do_prints = 0;
  int run_checks = 1;
  if (!b->standby) {
    layout(b, do_prints, run_checks);
  } else {
    // make dummy lattice that doesn't communicate with anyone.
    standby_layout(b);
  }


  /* Now we need a mapping that tells us how to relate fields on the original lattice
  * to those on the blocked lattice. This is stored in the "blocklist", which
  * tells each node where they need to send which sites, and tells nodes on the
  * blocked lattice what they need to receive */

  make_blocklists(l, b, block_dir);
  /* all nodes need to call this, because "standby" nodes still need to
  * send their fields to non-standby nodes*/

  // finally, test
  test_blocking(l, b, block_dir);

  barrier(l->comm);

}


/* Construct comlists for communication between the original lattice
* and the blocked lattice, so that all nodes know where to send/receive
* their blocked fields.
* l = original lattice, b = blocked lattice
*
* Specifically, this routine creates and fills in send_to structure
* in 'l', and recv_from structure in 'b', and does not touch the other
* sendrecv struct. Otherwise this is quite similar to the routine in comms.c.
*
* Assumes that all blocked nodes have the same
* number of sites and that halos come after real sites!
*/
void make_blocklists(lattice* l, lattice* b, int const* block_dir) {

  long block_sites = b->sites;

  // keep track of how many sites the blocked lattice has received
  long sites_received = 0;

  #ifdef MPI

    /* broadcast the number of sites in the blocked nodes.
    * this is necessary because the nodes on standby really have 0 sites,
    * but routines below are easier if we treat everyone more or less equally */

    bcast_long(&block_sites, l->comm); // root node always operates on the new lattice

    if (block_sites != b->sites && !b->standby) {
      printf("Error in blocklist.c!! nodes have different number of sites \n");
    }

    /* For nodes NOT on standby on the blocked lattice,
    * send the coordinate tables to all nodes that work on the original lattice */
    long** coord_buf[l->size];
    MPI_Request coord_req[l->size];
    int coord_tag = 0;

    // send blocked nodes a list of sites that they will receive from me
    long* sites_buf[l->size];
    MPI_Request sites_req[l->size];
    int sites_tag = 1;

    /* Initialization: allocate everything here even if some nodes do not need
    * all of the buffers. This makes it easier to keep track of freeing.
    * Also initialize the requests to MPI_REQUEST_NULL;
    * this ensures that MPI_Wait works even if we never sent anything
    */
    for (int r=0; r<l->size; r++) {
      coord_buf[r] = alloc_latticetable(b->dim, block_sites);
      coord_req[r] = MPI_REQUEST_NULL;
      sites_buf[r] = malloc(l->sites * sizeof(*sites_buf[r])); // will realloc later
      sites_req[r] = MPI_REQUEST_NULL;
    }

    // if not on standby, send the coordinate table to everyone
    if (!b->standby) {
      for (int r=0; r<l->size; r++) {

        if (r == l->rank) {
          // don't send own stuff
          continue;
        }

        // copy site coordinates (excl. halos)
        for (long j=0; j<block_sites; j++) {
          for (int dir=0; dir<b->dim; dir++) {
            coord_buf[r][j][dir] = b->coords[j][dir];
          }
        }

        int size = b->dim * block_sites;
        // send to everyone in l->comm
        MPI_Isend(&coord_buf[r][0][0], size, MPI_LONG, r, coord_tag, l->comm, &coord_req[r]);
        //printf("node %d: sent blocked coordinates to %d\n", l->rank, r);
      }
    }

  #endif // end MPI

  long tot = 0; // keep track of how many sites have been mapped

  /* For each node on the new lattice, collect their coordinate tables */
  long** blocked_coords = alloc_latticetable(b->dim, b->sites);

  // loop over MPI ranks on the blocked lattice
  for (int r=0; r<b->size; r++) {

    long sends = 0;
    long recvlist[l->sites];

    if (r == b->rank) {
      // own node and not on standby, just copy the coords table
      for (long j=0; j<block_sites; j++) {
        for (int dir=0; dir<b->dim; dir++) {
          blocked_coords[j][dir] = b->coords[j][dir];
        }
      }
    } else {
      #ifdef MPI
        // other node, so receive their coords
        int size;
        MPI_Status status;
        MPI_Probe(r, coord_tag, l->comm, &status);
        MPI_Get_count(&status, MPI_LONG, &size);

        MPI_Recv(&blocked_coords[0][0], size, MPI_LONG, r, coord_tag, l->comm, MPI_STATUS_IGNORE);
      #else
        printf("Should not get here!! in blocking.c (2)\n");
        die(941);
      #endif
    }

    /* For each site on the original lattice, find the corresponding coordinates
    on the blocked lattice and check if that site exists in the node with rank = r */
    for (long i=0; i<l->sites; i++) {

      int done = 0;
      int coords[l->dim];
      // find new coords
      for (int dir=0; dir<l->dim; dir++) {
        coords[dir] = l->coords[i][dir];
        // odd coordinate in a blocked direction?
        if (block_dir[dir] && coords[dir] % 2 != 0) {
          // no need to send these
          done = 1;
          break;
        } else if (block_dir[dir] && coords[dir] % 2 == 0) {
          // even site in blocked direction, coordinate gets halved
          coords[dir] /= 2;
        } else {
          // no blocking in this direction, so coordinate is unchanged
        }

      } // end dir

      // now check if this site exists on the new lattice
      if (!done) {
        for (long j=0; j<block_sites; j++) {
          // j loop: all sites on the blocked node
          int match = 1;
          for (int dir=0; dir<l->dim; dir++) {
            if (coords[dir] != blocked_coords[j][dir]) {
              match = 0;
              break;
            }
          }

          if (match) {
            // site found
            tot++;
            /* add to send blocklist in my node, store their site index */
            addto_comlist(&l->blocklist, r, i, SEND, l->parity[i], b->sites);
            recvlist[sends] = j;
            #ifdef MPI
              sites_buf[r][sends] = j;
            #endif
            sends++;
          }

        } // end j
      }

    } // end i

    /* all sites checked, now let the other node know what they need to receive.
    * If own node, add to blocklist directly. Note that the receives are
    * in the same order as my sends. */


    // list of sites in the receiving node are now stored in recvlist; send it
    if (sends > 0) {

      if (r == l->rank) {
        for (long s=0; s<sends; s++) {
          long mysite = recvlist[s];
          addto_comlist(&b->blocklist, r, mysite, RECV, b->parity[mysite], block_sites);
          sites_received++;
        }

      } else {
        #ifdef MPI
          // not my node, so send the site list
          memcpy(&sites_buf[r][0], &recvlist[0], sends * sizeof(*sites_buf[r]));
          sites_buf[r] = realloc(sites_buf[r], sends * sizeof(*sites_buf[r]));
          if (sites_buf[r] == NULL) {
            printf("Failed to realloc buffer! in blocking.c\n");
            die(4221);
          }
          MPI_Isend(&sites_buf[r][0], sends, MPI_LONG, r, sites_tag, l->comm, &sites_req[r]);
          //printf("Node %d: sent blocked sites to %d\n", l->rank, r);
        #else
          printf("Should not get here!! in blocking.c\n");
        #endif
      }
    }

  } // end rank loop

  // send lists done, now fill in the receives for nodes using the blocked lattice
  #ifdef MPI

  if (!b->standby) {

    while (sites_received < b->sites) {

      // receive the site list from ANY node
      MPI_Status status;
      int count;
      MPI_Probe(MPI_ANY_SOURCE, sites_tag, l->comm, &status);
      MPI_Get_count(&status, MPI_LONG, &count);

      long sitelist[count];
      MPI_Recv(&sitelist[0], count, MPI_LONG, status.MPI_SOURCE, sites_tag, l->comm, MPI_STATUS_IGNORE);

      for (int j=0; j<count; j++) {
        addto_comlist(&b->blocklist, status.MPI_SOURCE, sitelist[j], RECV, b->parity[j], block_sites);
      }
      sites_received += count;

    } // end r
  }
  #endif

  if (sites_received != b->sites && !b->standby) {
    printf("WARNING from node %d: did not map all blocked sites (in blocking.c)\n", l->rank);
  }


  // receives done, now free the buffers and temp arrays
  #ifdef MPI
    for (int r=0; r<l->size; r++) {
      MPI_Wait(&coord_req[r], MPI_STATUS_IGNORE);
      MPI_Wait(&sites_req[r], MPI_STATUS_IGNORE);
      free_latticetable(coord_buf[r]);
      free(sites_buf[r]);
    }
  #endif

  free_latticetable(blocked_coords);


  // check that we didn't miss any sites on the blocked lattice
  tot = reduce_sum_long(tot, l->comm);
  if (l->rank == 0 && tot != b->vol) {
    printf("blocking.c: site mapping failed!! \n");
    die(-1921);
  }

  // finally, realloc blocklists
  realloc_comlist(&l->blocklist, SEND);
  realloc_comlist(&b->blocklist, RECV);

}

/* Make a dummy lattice struct for nodes on standby. This allocs everything that
* normal layout() does, so the lattice can be freed using free_lattice(),
* but the resulting lattice is essentially empty.
*/
void standby_layout(lattice* l) {

  if (!l->standby) {
    return;
  }

  l->sliceL = malloc(l->dim * sizeof(*(l->sliceL)));
  l->nslices = malloc(l->dim * sizeof(*(l->nslices)));
  l->offset = malloc(l->dim * sizeof(*(l->offset)));

  l->sites = 1;
	l->halos = 0;
	l->sites_total = 1;
  l->evensites = 1; l->oddsites = 0;
  l->evenhalos = 1; l->oddhalos = 0;

  alloc_lattice_arrays(l, l->sites_total);
  l->parity[0] = EVEN;
  for (int dir=0; dir<l->dim; dir++) {
    l->coords[0][dir] = -1;
    l->sliceL[dir] = 1;
    l->nslices[dir] = 1;
    l->offset[dir] = 0;
  }

  l->longest_dir = 0;
  l->sites_per_coord = malloc(1 * sizeof(*l->sites_per_coord));
  l->sites_per_coord[0] = 1;
  l->sites_at_coord = malloc(1 * sizeof(*l->sites_at_coord));
  l->sites_at_coord[0] = alloc_latticetable(l->sites_per_coord[0], l->sliceL[0]);

  #ifdef MEASURE_Z
    l->sites_per_z = 1;
    l->z_dir = 0;
    l->site_at_z = alloc_latticetable(l->sites_per_z, l->sliceL[0]);
  #endif

  alloc_comlist(&l->comlist, 1); // sets sends,recvs = 0
}


/* Copy a smeared field on the original lattice 'l' to a blocked field
* on blocked lattice 'b'. Only does communications between distinct nodes,
* for transfering blocked fields between the same node use block_fields_ownnode().
*/
void transfer_blocked_field(lattice* l, lattice* b, double** field, double** field_b, int dofs) {

  #ifdef MPI

    int sends = l->blocklist.sends;
    if (sends <= 0) {
      return;
    }

    double start, end;
	  start = clock();

    // first nonblocking sends to everyone
    MPI_Request req[sends];
    for (int k=0; k<sends; k++) {
      req[k] = MPI_REQUEST_NULL;

      sendrecv_struct* send = &l->blocklist.send_to[k];
      int rank = send->node;
      if (rank == l->rank) {
       // don't send to self. these are all done by block_fields_ownnode(),
       // which is assumed to be called separately.
       continue;
      }

      // allocates send buffer and does nonblocking send:
      send_field(send, l->comm, &req[k], EVENODD, field, dofs);
    }

    // sends done, then blocking receives in the blocked nodes
    if (!b->standby) {
      int recvs = b->blocklist.recvs;
      for (int m=0; m<recvs; m++) {
        sendrecv_struct* recv = &b->blocklist.recv_from[m];
        if (recv->node == b->rank) {
          continue; // don't receive from self, call block_fields_ownnode() elsewhere
        }
        recv_field(recv, l->comm, EVENODD, field_b, dofs);
      }
    }

    // recvs done, now wait and free send buffers
    for (int k=0; k<sends; k++) {
      sendrecv_struct* send = &l->blocklist.send_to[k];
      if (send->node == l->rank) continue; // buffer for own node never allocated

      MPI_Wait(&req[k], MPI_STATUS_IGNORE);
      free(send->buf);
    }

    end = clock();
    Global_comms_time += (double)(end - start) / CLOCKS_PER_SEC;

    #endif
    // do nothing in serial. Use block_fields_ownnode() instead.
}


/* Same as transfer_blocked_field(), but for a gauge field of type double***. */
void transfer_blocked_gaugefield(lattice* l, lattice* b, double*** field, double*** field_b, int dofs, int dir) {

  #ifdef MPI

    int sends = l->blocklist.sends;
    if (sends <= 0) {
      return;
    }

    double start, end;
	  start = clock();

    // first nonblocking sends to everyone
    MPI_Request req[sends];
    for (int k=0; k<sends; k++) {
      req[k] = MPI_REQUEST_NULL;

      sendrecv_struct* send = &l->blocklist.send_to[k];
      int rank = send->node;
      if (rank == l->rank) {
       // don't send to self. these are all done by block_fields_ownnode(),
       // which is assumed to be called separately.
       continue;
      }

      // allocates send buffer and does nonblocking send:
      send_gaugefield(send, l->comm, &req[k], EVENODD, field, dofs, dir);
    }

    // sends done, then blocking receives in the blocked nodes
    if (!b->standby) {
      int recvs = b->blocklist.recvs;
      for (int m=0; m<recvs; m++) {
        sendrecv_struct* recv = &b->blocklist.recv_from[m];
        if (recv->node == b->rank) {
          continue; // don't receive from self, call block_fields_ownnode() elsewhere
        }

        recv_gaugefield(recv, l->comm, EVENODD, field_b, dofs, dir);
      }
    }

    // recvs done, now wait and free send buffers
    for (int k=0; k<sends; k++) {
      sendrecv_struct* send = &l->blocklist.send_to[k];
      if (send->node == l->rank) continue;

      MPI_Wait(&req[k], MPI_STATUS_IGNORE);
      free(send->buf);
    }

    end = clock();
    Global_comms_time += (double)(end - start) / CLOCKS_PER_SEC;

    #endif
    // do nothing in serial. Use block_fields_ownnode() instead.
}

/* Transfer smeared fields on the full lattice (f_smeared) to
* fields on the blocked lattice with less sites (f_blocked).
*/
void make_blocked_fields(lattice* l, lattice* b, fields const* f_smeared, fields* f_blocked) {

  // first do all "sends" within my own node
  block_fields_ownnode(l, b, f_smeared, f_blocked);

  #ifdef MPI
    // then contributions from other nodes

    for (int dir=0; dir<l->dim; dir++) {
      transfer_blocked_gaugefield(l, b, f_smeared->su2link, f_blocked->su2link, SU2LINK, dir);
    }
    #ifdef U1
      transfer_blocked_field(l, b, f_smeared->u1link, f_blocked->u1link, l->dim);
    #endif

    #ifdef HIGGS
      transfer_blocked_field(l, b, f_smeared->su2doublet, f_blocked->su2doublet, SU2DB);
    #endif
    #ifdef TRIPLET
      transfer_blocked_field(l, b, f_smeared->su2triplet, f_blocked->su2triplet, SU2TRIP);
    #endif

  #endif // end MPI

  // done, now just sync halos on the blocked lattice
  sync_halos(b, f_blocked);
}

/* Copy smeared fields to the corresponding location on a blocked lattice,
* but for own node only, without communications. Used by make_blocked_fields()
* to avoid communicating with self (and in serial).
*/
void block_fields_ownnode(lattice const* l, lattice const* b, fields const* f_smeared, fields* f_blocked) {

  if (b->standby) {
    return;
  }

  // first find own node in the comlists
  sendrecv_struct* send = NULL;
  sendrecv_struct* recv = NULL;
  int sends = l->blocklist.sends;
  for (int k=0; k<sends; k++) {
    if (l->blocklist.send_to[k].node == l->rank) {
      send = &l->blocklist.send_to[k];
      break;
    }
  }
  if (send == NULL) {
    return; // nothing to send to self
  }

  // same for receives
  int recvs = b->blocklist.recvs;
  for (int k=0; k<recvs; k++) {
    if (b->blocklist.recv_from[k].node == b->rank) {
      recv = &b->blocklist.recv_from[k];
      break;
    }
  }
  if (recv == NULL) {
    // should not get here, because above we found stuff to send to self
    printf("Error (1) in block_fields_ownnode()! (blocking.c)\n");
    die(1002);
  }

  if (send->sites != recv->sites) {
    printf("Error (2) in block_fields_ownnode()! (blocking.c)\n");
    die(1003);
  }

  for (long i=0; i<send->sites; i++) {
    long site_l = send->sitelist[i]; // site on the original lattice
    long site_b = recv->sitelist[i]; // site on the blocked lattice

    // now copy the smeared fields into the right places in f_blocked
    for (int dir=0; dir<l->dim; dir++) {
      memcpy(f_blocked->su2link[site_b][dir], f_smeared->su2link[site_l][dir], SU2LINK * sizeof(f_blocked->su2link[site_b][dir][0]));
      #ifdef U1
        f_blocked->u1link[site_b][dir] = f_smeared->u1link[site_l][dir];
      #endif
    }
    #ifdef HIGGS
      memcpy(f_blocked->su2doublet[site_b], f_smeared->su2doublet[site_l], SU2DB * sizeof(f_blocked->su2doublet[site_b][0]));
    #endif
    #ifdef TRIPLET
      memcpy(f_blocked->su2triplet[site_b], f_smeared->su2triplet[site_l], SU2TRIP * sizeof(f_blocked->su2triplet[site_b][0]));
    #endif
  }
  // self blocking done
}

/* Test routine: checks that fields on the original lattice end up
* at the correct sites on the blocked lattice. I do this by assigning the field
* a unique number at each site (using Cantor pairing). See test_comms_individual()
* for a similar routine. */
void test_blocking(lattice* l, lattice* b, int const* block_dir) {
  fields f;
  fields f_b;

  alloc_fields(l, &f);
  alloc_fields(b, &f_b);

  long x, y;
  // use just the SU(2) links for testing. no distinction between different dirs though
  for (long i=0; i<l->sites; i++) {
    // calculate a unique value using Cantor pairing
    y = l->coords[i][0];
		for (int dir=1; dir<l->dim; dir++) {
			x = l->coords[i][dir];
			y = y + 0.5 * (x + y + 1) * (x + y);
		}

    for (int dir=0; dir<l->dim; dir++) {
      for (int k=0; k<SU2LINK; k++) {
        // Obtain base number from Cantor, set different decimals for different components
        f.su2link[i][dir][k] = y + (double) k / SU2LINK;
        // initialize the blocked fields to something easy
        if (i < b->sites) {
          f_b.su2link[i][dir][k] = -1.0;
        }
      }
    }
  }

  make_blocked_fields(l, b, &f, &f_b);

  if (!b->standby) {
    // check the new values in f_b
    for (long i=0; i<b->sites; i++) {

      // what are the coords BEFORE blocking?
      long coords[b->dim];
      for (int dir=0; dir<b->dim; dir++) {
        if (block_dir[dir]) {
          coords[dir] = b->coords[i][dir] * 2;
        } else {
          coords[dir] = b->coords[i][dir];
        }
      }

      // repeat the Cantor pairing
      y = coords[0];
  		for (int dir=1; dir<b->dim; dir++) {
  			x = coords[dir];
  			y = y + 0.5 * (x + y + 1) * (x + y);
  		}

      for (int dir=0; dir<b->dim; dir++) {
        for (int k=0; k<SU2LINK; k++) {
          // predicted value:
          double val = y + (double) k / SU2LINK;

          if (fabs(f_b.su2link[i][dir][k] - val) > 0.0001) {
            printf("Node %d: error in test_blocking at site %ld, dir %d!! field val is %lf; was supposed to be %lf (blocking.c)\n"
                , b->rank, i, dir, f_b.su2link[i][dir][k], val);
          }
        }
      }

    } // end i

  }

  free_fields(l, &f);
  free_fields(b, &f_b);
}

#endif

/* Send and receive a field with dofs components to/from all my neighbors.
* Uses comlist.send_to struct to figure out which sites need to be sent where,
* and inludes only sites with parity as defined in the argument.
* NB! Can result in a deadlock if the sendrecv_structs in comlist are not ordered
* by node rank!
*
* Return value is the time spent on MPI_Sendrecv.
*
* We end up malloc'ing and freeing the buffer many times in one update.
* This should still be faster than realloc, since we don't need to copy 
* data from the old buffer to the new one...
*/
double sendrecv_field(comlist_struct* comlist, char parity, double** field, int dofs) {
	
	double start, end, time = 0.0;
	
	long send_offset, recv_offset, send_max, recv_max,index;
	long send_count, recv_count; // how many MPI_DOUBLEs do we send/receive
	int tag =0;

	double* send_buf;
	double* recv_buf;
	sendrecv_struct* send;
	sendrecv_struct* recv;

	for (int k=0; k<comlist->neighbors; k++) {

		send = &(comlist->send_to[k]);
		recv = &(comlist->recv_from[k]);
		// sites with parity = EVEN come first
		if (parity == EVEN) {
			// update even
			send_offset = 0;
			recv_offset = 0;
			send_count = send->even;
			recv_count = recv->even;
			send_max = send->even;
			recv_max = recv->even;
		} else {
			// update odd
			send_offset = send->even;
			recv_offset = recv->even;
			send_count = send->odd;
			recv_count = recv->odd;
			send_max = send->sites;
			recv_max = recv->sites;
		}
		send_count *= dofs;
		recv_count *= dofs;


		send_buf = malloc(send_count * sizeof(*send_buf));
		recv_buf = malloc(recv_count * sizeof(*recv_buf));

		// copy field values into send_buf
		long j = 0;
		for (long i = send_offset; i<send_max; i++) {
			index = send->sitelist[i];

			for (int d=0; d<dofs; d++) {
				send_buf[j + d] = field[index][d];
			}
			j += dofs;
		}

		start = clock();
		MPI_Sendrecv(send_buf, send_count, MPI_DOUBLE, send->node, tag, recv_buf, recv_count, MPI_DOUBLE, recv->node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		end = clock();
		
		time += ((double) (end - start)) / CLOCKS_PER_SEC;
		
		// copy recv_buf into field array
		j = 0;
		for (long i = recv_offset; i<recv_max; i++) {
			index = recv->sitelist[i];
			
			for (int d=0; d<dofs; d++) {
				field[index][d] = recv_buf[j + d];
			}
			j += dofs;
		}

		free(send_buf);
		free(recv_buf);

	} // end loop over neighbors

	return time;
}



/* Same as sendrecv_field, but for gauge fields.
* Specifically this routine updates gauge halos in direction dir
* for sites with given parity.
*/
double sendrecv_gaugefield(comlist_struct* comlist, char parity, double*** field, int dofs, int dir) {
	
	double start, end, time = 0.0;
	
	long send_offset, recv_offset, send_max, recv_max,index;
	long send_count, recv_count; // how many MPI_DOUBLEs do we send/receive
	int tag=0;

	double* send_buf;
	double* recv_buf;
	sendrecv_struct* send;
	sendrecv_struct* recv;

	for (int k=0; k<comlist->neighbors; k++) {

		send = &(comlist->send_to[k]);
		recv = &(comlist->recv_from[k]);
		// sites with parity = EVEN come first
		if (parity == EVEN) {
			// update even
			send_offset = 0;
			recv_offset = 0;
			send_count = send->even;
			recv_count = recv->even;
			send_max = send->even;
			recv_max = recv->even;
		} else {
			// update odd
			send_offset = send->even;
			recv_offset = recv->even;
			send_count = send->odd;
			recv_count = recv->odd;
			send_max = send->sites;
			recv_max = recv->sites;
		}
		send_count *= dofs;
		recv_count *= dofs;

		send_buf = malloc(send_count * sizeof(*send_buf));
		recv_buf = malloc(recv_count * sizeof(*recv_buf));

		// copy field values into send_buf. 
		long j = 0;
		for (long i = send_offset; i<send_max; i++) {
			index = send->sitelist[i];
			for (int d=0; d<dofs; d++) {		
				// this bit is slightly different from sendrecv_field()
				send_buf[j + d] = field[index][dir][d];
			}
			j += dofs;
		}

		start = clock();
		MPI_Sendrecv(send_buf, send_count, MPI_DOUBLE, send->node, tag, recv_buf, recv_count, MPI_DOUBLE, recv->node, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		end = clock();
		
		time += ((double) (end - start)) / CLOCKS_PER_SEC;
		
		// copy recv_buf into field array
		j = 0;
		for (long i = recv_offset; i<recv_max; i++) {
			index = recv->sitelist[i];
			for (int d=0; d<dofs; d++) {
				// this bit is slightly different from sendrecv_field(){
				field[index][dir][d] = recv_buf[j + d];
			}
			j += dofs;
		}

		free(send_buf);
		free(recv_buf);

	} // end loop over neighbors

	return time;
}

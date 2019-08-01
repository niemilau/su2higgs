#ifndef COMMS_H
#define COMMS_H

// did we compile with MPI?
#ifdef MPI
	#include <mpi.h>
#endif

// TODO make this completely separate from su2.h and its .c files.

typedef struct {
  int node;               /* node index to send to/received from */
  long even, odd, sites;   /* number of sites to be sent or received*/
  long *sitelist;   	/* list of sites to be sent or received */
	double* buf; 			/* buffer for sending and receiving */
} sendrecv_struct;


typedef struct {
  sendrecv_struct * send_to; // where to send
  sendrecv_struct * recv_from; // where to collect from
  int neighbors;
	long max_sendrecv;
} comlist_struct;

// comms.c 


#endif

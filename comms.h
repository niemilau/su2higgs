#ifndef COMMS_H
#define COMMS_H

#define SEND 0
#define RECV 1


typedef struct {
  int node;               /* node index to send to/received from */
  long even, odd, sites;   /* number of sites to be sent or received */
  long *sitelist;   	/* list of sites to be sent or received */
	double* buf; 			/* buffer for sending and receiving */
} sendrecv_struct;


typedef struct {
  sendrecv_struct * send_to; // where to send
  sendrecv_struct * recv_from; // where to collect from
  int sends; // how many nodes to send to
	int recvs; // how many nodes to receive from
	/* for ordinary comlists, recvs and sends are the same:
	* each neighbor both sends and receives.
	* For comlists used for blocking, they can differ. */

} comlist_struct;

// comms.c


#endif

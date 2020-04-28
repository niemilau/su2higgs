
/** @file checkpoint.c
*
* Routines for storing lattice configuration in a file,
* and printing useful information about the simulation.
*
* TODO
*
*/

#include "su2.h"

/* Print acceptance rates of relevant algorithms
* in a compact form. Called at each checkpoint.
*/
void print_acceptance(params p, counters c) {
	int print_gauge = 0;
	#ifdef TRIPLET
		print_gauge = 1;
	#endif
	printf("\n-- Acceptance rates --\n	");
	if (p.algorithm_su2link == METROPOLIS || print_gauge) {
		printf("SU(2) link %.2lf%%, ",
			100.0*c.accepted_su2link/c.total_su2link);
	}

	#ifdef U1
		if (p.algorithm_u1link == METROPOLIS) {
			printf("U(1) link %.2lf%%, ",
				100.0*c.accepted_u1link/c.total_u1link);
		}
	#endif

	#ifdef HIGGS
		if (p.algorithm_su2doublet == METROPOLIS) {
			printf("Higgs %.2lf%%, ",
				100.0*c.accepted_doublet/c.total_doublet);

		} else if (p.algorithm_su2doublet == OVERRELAX) {
			printf("Higgs overrelax %.2lf%%, Higgs Metropolis %.2lf%%, ",
				100.0*c.acc_overrelax_doublet/c.total_overrelax_doublet,
				100.0*c.accepted_doublet/c.total_doublet);
		}
	#endif

	#ifdef TRIPLET
		if (p.algorithm_su2triplet == METROPOLIS) {
			printf("Triplet %.2lf%%, ",
				100.0*c.accepted_triplet/c.total_triplet);
		} else if (p.algorithm_su2triplet == OVERRELAX) {
			printf("Triplet overrelax %.2lf%%, Triplet Metropolis %.2lf%%, ",
				100.0*c.acc_overrelax_triplet/c.total_overrelax_triplet,
				100.0*c.accepted_triplet/c.total_triplet);
		}
	#endif
	if (p.multicanonical) {
		printf("multicanonical %.2lf%%",
				100.0*c.accepted_muca/c.total_muca);
	}
	printf("\n");

}

/* Write all fields to a file.
* Also stores lattice dimensions and current iteration number.
* Theory parameters such as beta_G and masses are NOT stored!
* Neither are model-specific acceptance rates.
*/
void save_lattice(params p, fields f, counters c) {

	FILE *file;
	if (p.rank == 0) {
		file = fopen(p.latticefile, "wb");
		// first line: p.size p.dim L1 L2 ... Ln
		fwrite(&p.size, sizeof(p.size), 1, file);
		fwrite(&p.dim, sizeof(p.dim), 1, file);
		fwrite(p.L, sizeof(p.L[0]), p.dim, file);

		// second line: iteration total_time comms_time
		fwrite(&c.iter, sizeof(c.iter), 1, file);
		fwrite(&c.total_time, sizeof(c.total_time), 1, file);
		fwrite(&c.comms_time, sizeof(c.comms_time), 1, file);
	}

	// fields. file is only open in root node, so others cannot use it here.
	write_field(p, file, &f.su2link[0][0][0], p.dim * SU2LINK);
	#ifdef U1
	write_field(p, file, &f.u1link[0][0], p.dim);
	#endif
	#ifdef HIGGS
	write_field(p, file, &f.su2doublet[0][0], SU2DB);
	#endif
	#ifdef TRIPLET
	write_field(p, file, &f.su2triplet[0][0], SU2TRIP);
	#endif

	if (p.rank == 0) {
		fclose(file);
		printf("Written fields to %s.\n", p.latticefile);
	}

}


/* Read all fields from the file created by save_lattice().
* Note that this assumes that layouting is exactly the same
* as when the file was written, including the relative "locations"
* of MPI nodes. Hence, perform a crosscheck here.
* ALL nodes read the first few lines of the latticefile so that
* counters and iteration number can be kept in sync, while only the root node
* reads fields and distributes them to others.
*/
void load_lattice(params const* p, fields* f, counters* c, comlist_struct* comlist) {

	FILE *file;

	file = fopen(p->latticefile, "rb");
	// first line: p.size p.dim L1 L2 ... Ln
	// note that these are int type (need to be careful with this)
	int dim, size, read = 0;
	// compiler gives warning if return value is not used, so count the reads here
	read += fread(&size, sizeof(size), 1, file);
	read += fread(&dim, sizeof(dim), 1, file);

	int L[dim];
	read += fread(L, sizeof(L[0]), dim, file);

	if (!read) {
		// did not read anything...
		printf0(*p, "Error reading latticefile!\n");
		die(500);
	}

	int ok = 1;
	// check that dimensions of the lattice file match those in our config
	if (p->dim != dim || p->size != size)
		ok = 0;

	if (ok) {
		for (int d=0; d<p->dim; d++) {
			if (L[d] != p->L[d])
				ok = 0;
		}
	}

	if (!ok) {
		printf0(*p, "Dimensions in latticefile do not match! Got:\n");
		printf0(*p, " 	MPI size %d, dimension %d, volume ", size, dim);
		for (int d=0; d<dim; d++) {
			printf0(*p, "%d x ", L[d]);
		}
		printf0(*p, "\b\b \b\n\nWas supposed to be: \n");
		printf0(*p, "	MPI size %d, dimension %d, volume ", p->size, p->dim);
		for (int d=0; d<p->dim; d++) {
			printf0(*p, "%d x ", p->L[d]);
		}
		printf0(*p, "\b\b \b\n");
		die(501);
	}

	read = 0;
	// second line: iteration total_time comms_time
	read += fread(&c->iter, sizeof(c->iter), 1, file);
	read += fread(&c->total_time, sizeof(c->total_time), 1, file);
	read += fread(&c->comms_time, sizeof(c->comms_time), 1, file);

	// if not root node, can close the file here
	if (p->rank != 0)
		fclose(file);

	// Read fields. Ordering HAS to be same as in save_lattice()
	read_field(*p, file, &f->su2link[0][0][0], p->dim * SU2LINK);
	#ifdef U1
	read_field(*p, file, &f->u1link[0][0], p->dim);
	#endif
	#ifdef HIGGS
	read_field(*p, file, &f->su2doublet[0][0], SU2DB);
	#endif
	#ifdef TRIPLET
	read_field(*p, file, &f->su2triplet[0][0], SU2TRIP);
	#endif

	// finally, sync all halo fields; these were not loaded from the file
	sync_halos(f, p, comlist);

	if (p->rank == 0) {
		fclose(file);
	}
}

#ifdef MPI

/* Write a field to file. All nodes send their entire field array
* to the root node, which performs the writing in the order of MPI ranks.
* Note that we also write all halo fields for each node.
* Routine assumes that file is open only in the root node.
*
* field = pointer to first element of the field array (e.g. &f.su2link[0][0][0])
* size = how many components the field has at a single site
*/
void write_field(params p, FILE *file, double *field, int size) {

	MPI_Barrier(MPI_COMM_WORLD);
	// how many sites in my node, including halos?
	long max = p.sites_total * size;

	if (p.rank == 0) {
		// start by writing the field in root node
		fwrite(field, sizeof(*field), max, file);
	} else {
		// other nodes: send field to root
		MPI_Send(field, max, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}


	if (p.rank == 0) {
		// in root node, receive from all nodes one at a time
		double* buf = malloc(max * sizeof(*buf));

		for (int rank=1; rank<p.size; rank++) {
			// Probe message size here, in case the other node has different number of sites
			MPI_Status status;
			int count;
			MPI_Probe(rank, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			long newmax = count; // got to hope this does not overflow...

			if (newmax != max) {
				buf = realloc(buf, newmax * sizeof(*buf));
				if (buf == NULL) {
					printf("WARNING! Failed to realloc buffer in write_field()\n");
				}
			}

			MPI_Recv(buf, newmax, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			fwrite(buf, sizeof(*buf), newmax, file);
			max = newmax;
		}

		free(buf);
	}
	MPI_Barrier(MPI_COMM_WORLD);

}

/* Reads a field from latticefile. As in write_field(), reading is done
* by the root node, which then sends the field to the other nodes.
* Note that we send and receive TWO messages per node, so need to use tags.
*/
void read_field(params p, FILE *file, double *field, int size) {

	MPI_Barrier(MPI_COMM_WORLD);
	int maxtag = 1, fieldtag = 2;
	// how many sites in my node, including halos?
	long max = p.sites_total * size;

	if (p.rank == 0) {
		// first read the field belonging to root node
		long read = fread(field, sizeof(*field), max, file);
		if (read != max) {
			printf("Error reading field for root node!\n");
			die(505);
		}
	} else {
		// other nodes: send max to root, in case number of sites is different
		MPI_Send(&max, 1, MPI_LONG, 0, maxtag, MPI_COMM_WORLD);
		// start receiving the field
		MPI_Recv(field, max, MPI_DOUBLE, 0, fieldtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if (p.rank == 0) {

		double* buf = malloc(max * sizeof(*buf));

		for (int rank=1; rank<p.size; rank++) {
			long newmax;
			// in root node, first receive their max
			MPI_Recv(&newmax, 1, MPI_LONG, rank, maxtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			if (newmax != max) {
				buf = realloc(buf, newmax * sizeof(*buf));
				if (buf == NULL) {
					printf("WARNING! Failed to realloc buffer in read_field()\n");
				}
			}
			max = newmax;

			// read in their field
			long read = fread(buf, sizeof(*buf), max, file);
			if (read != max) {
				printf("Error reading field for node %d!\n", rank);
				die(506);
			}

			// read OK, send to the other node
			MPI_Send(buf, max, MPI_DOUBLE, rank, fieldtag, MPI_COMM_WORLD);
		}

		free(buf);
	}
	MPI_Barrier(MPI_COMM_WORLD);

}

#else // no MPI; simplified write and read routines

/* Write a field to latticefile.
* field = pointer to first element of the field array (e.g. &f.su2link[0][0][0])
* size = how many components the field has at a single site
*/
void write_field(params p, FILE *file, double *field, int size) {

	long max = p.sites_total * size;
	// write the whole field at once
	fwrite(field, sizeof(*field), max, file);

}

/* Reads in a field, assuming the format in write_field().
*/
void read_field(params p, FILE *file, double *field, int size) {

	long max = p.sites_total * size;
	long read = fread(field, sizeof(*field), max, file);

	if (read != max) {
		printf0(p, "Error reading field!\n");
		die(505);
	}

}

#endif

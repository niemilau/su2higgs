
/** @file checkpoint.c
*
* Routines for storing lattice configuration in a file,
* and printing useful information about the simulation.
* Halos will not be stored, so need to sync those separately.
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

	#if (NHIGGS > 0)

		for (int db=0; db<NHIGGS; db++) {
			if (p.algorithm_su2doublet == METROPOLIS) {
				printf("Higgs#%d %.2lf%%, ", db+1,
					100.0*c.accepted_doublet[db]/c.total_doublet[db]);

			} else if (p.algorithm_su2doublet == OVERRELAX) {
				printf("Higgs#%d overrelax %.2lf%%, Higgs#%d Metropolis %.2lf%%, ",
					db+1, 100.0*c.acc_overrelax_doublet[db]/c.total_overrelax_doublet[db],
					db+1, 100.0*c.accepted_doublet[db]/c.total_doublet[db]);
			}
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

	#ifdef SINGLET
	if (p.algorithm_singlet == METROPOLIS) {
		printf("Singlet %.2lf%%, ",
			100.0*c.accepted_singlet/c.total_singlet);
	} else if (p.algorithm_singlet == OVERRELAX) {
		printf("Singlet overrelax %.2lf%%, Singlet Metropolis %.2lf%%, ",
			100.0*c.acc_overrelax_singlet/c.total_overrelax_singlet,
			100.0*c.accepted_singlet/c.total_singlet);
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
* Neither are model-specific acceptance rates. */
void save_lattice(lattice const* l, fields f, counters c, char* fname) {

	FILE *file;
	if (l->rank == 0) {
		file = fopen(fname, "wb");
		// first line: p.size p.dim L1 L2 ... Ln
		fwrite(&l->size, sizeof(l->size), 1, file);
		fwrite(&l->dim, sizeof(l->dim), 1, file);
		fwrite(l->L, sizeof(l->L[0]), l->dim, file);

		// second line: iteration total_time comms_time
		fwrite(&c.iter, sizeof(c.iter), 1, file);
		fwrite(&Global_total_time, sizeof(Global_total_time), 1, file);
		fwrite(&Global_comms_time, sizeof(Global_comms_time), 1, file);
	}

	// fields. file is only open in root node, so others cannot use it here.
	write_field(l, file, &f.su2link[0][0][0], l->dim * SU2LINK);
	#ifdef U1
		write_field(l, file, &f.u1link[0][0], l->dim);
	#endif
	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) write_field(l, file, &f.su2doublet[db][0][0], SU2DB);
	#endif

	#ifdef TRIPLET
		write_field(l, file, &f.su2triplet[0][0], SU2TRIP);
	#endif

	#ifdef SINGLET
		write_field(l, file, &f.singlet[0][0], 1);
	#endif

	if (l->rank == 0) {
		fclose(file);
		printf("Wrote fields to %s.\n", fname);
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
void load_lattice(lattice* l, fields* f, counters* c, char* fname) {

	FILE *file;

	file = fopen(fname, "rb");

	// first line: p.size p.dim L1 L2 ... Ln
	int dim, size, read = 0;
	// compiler gives warning if return value is not used, so count the reads here
	read += fread(&size, sizeof(size), 1, file);
	read += fread(&dim, sizeof(dim), 1, file);

	int L[dim];
	read += fread(L, sizeof(L[0]), dim, file);

	if (!read) {
		// did not read anything...
		printf0(*l, "Error reading latticefile!\n");
		die(500);
	}

	int ok = 1;
	// check that dimensions of the lattice file match those in our config
	if (l->dim != dim || l->size != size)
		ok = 0;

	if (ok) {
		for (int d=0; d<l->dim; d++) {
			if (L[d] != l->L[d])
				ok = 0;
		}
	}

	if (!ok) {
		printf0(*l, "Dimensions in latticefile do not match! Got:\n");
		printf0(*l, " 	MPI size %d, dimension %d, volume ", size, dim);
		for (int d=0; d<dim; d++) {
			printf0(*l, "%d x ", L[d]);
		}
		printf0(*l, "\b\b \b\n\nWas supposed to be: \n");
		printf0(*l, "	MPI size %d, dimension %d, volume ", l->size, l->dim);
		for (int d=0; d<l->dim; d++) {
			printf0(*l, "%d x ", l->L[d]);
		}
		printf0(*l, "\b\b \b\n");
		die(501);
	}

	read = 0;
	// second line: iteration total_time comms_time
	read += fread(&c->iter, sizeof(c->iter), 1, file);
	read += fread(&Global_total_time, sizeof(Global_total_time), 1, file);
	read += fread(&Global_comms_time, sizeof(Global_comms_time), 1, file);

	// if not root node, can close the file here
	if (l->rank != 0)
		fclose(file);

	// Read fields. Ordering HAS to be same as in save_lattice()
	read_field(l, file, &f->su2link[0][0][0], l->dim * SU2LINK);
	#ifdef U1
		read_field(l, file, &f->u1link[0][0], l->dim);
	#endif
	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) read_field(l, file, &f->su2doublet[db][0][0], SU2DB);
	#endif

	#ifdef TRIPLET
		read_field(l, file, &f->su2triplet[0][0], SU2TRIP);
	#endif

	#ifdef SINGLET
		read_field(l, file, &f->singlet[0][0], 1);
	#endif

	// finally, sync all halo fields; these were not loaded from the file
	sync_halos(l, f);

	if (l->rank == 0) {
		fclose(file);
	}
}

#ifdef MPI

/* Write a field to file. All nodes send their entire field array
* to the root node, which performs the writing in the order of MPI ranks.
* Note that halos are not stored.
* Routine assumes that file is open only in the root node.
*
* field = pointer to first element of the field array (e.g. &f.su2link[0][0][0])
* size = how many components the field has at a single site
*/
void write_field(lattice const* l, FILE *file, double *field, int size) {

	MPI_Barrier(MPI_COMM_WORLD);
	// how many sites in my node, excluding halos?
	long max = l->sites * size;

	if (l->rank == 0) {
		// start by writing the field in root node
		fwrite(field, sizeof(*field), max, file);
	} else {
		// other nodes: send field to root
		MPI_Send(field, max, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}


	if (l->rank == 0) {
		// in root node, receive from all nodes one at a time
		double* buf = malloc(max * sizeof(*buf));

		for (int rank=1; rank<l->size; rank++) {
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
void read_field(lattice const* l, FILE *file, double *field, int size) {

	MPI_Barrier(MPI_COMM_WORLD);
	int maxtag = 1, fieldtag = 2;
	// how many sites in my node
	long max = l->sites * size;

	if (l->rank == 0) {
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

	if (l->rank == 0) {

		double* buf = malloc(max * sizeof(*buf));

		for (int rank=1; rank<l->size; rank++) {
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
void write_field(lattice const* l, FILE *file, double *field, int size) {

	long max = l->sites * size;
	// write the whole field at once
	fwrite(field, sizeof(*field), max, file);

}

/* Reads in a field, assuming the format in write_field().
*/
void read_field(lattice const* l, FILE *file, double *field, int size) {

	long max = l->sites * size;
	long read = fread(field, sizeof(*field), max, file);

	if (read != max) {
		printf0(*l, "Error reading field!\n");
		die(505);
	}

}

#endif

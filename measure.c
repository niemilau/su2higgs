/** @file measure.c
*
* Routines for measuring volume averages and the total action,
* and writing them to a file.
*
*
* TODO
*
*/

#include "su2.h"

/* Print data labels for measurements (into a separate label file)
*/
void print_labels() {
	FILE* f = fopen("labels", "w+");

	int k = 1;
	fprintf(f, "%d weight\n", k); k++; // multicanonical weight
	fprintf(f, "%d muca param\n", k); k++; // multicanonical order parameter value
	fprintf(f, "%d action\n", k); k++;
	fprintf(f, "%d SU(2) Wilson\n", k); k++;
	#ifdef HIGGS
		fprintf(f, "%d hopping_phi\n", k); k++;
		fprintf(f, "%d phi^2\n", k); k++;
		fprintf(f, "%d phi^4\n", k); k++;
	#endif
	#ifdef TRIPLET
		fprintf(f, "%d hopping_Sigma\n", k); k++;
		fprintf(f, "%d Sigma^2\n", k); k++;
		fprintf(f, "%d Sigma^4\n", k); k++;
	#endif
	#if defined HIGGS && defined TRIPLET
		fprintf(f, "%d phi^2 Sigma^2\n", k); k++;
	#endif
	#ifdef U1
		fprintf(f, "%d U(1) Wilson\n", k); k++;
	#endif
	#ifdef TRIPLET
		fprintf(f, "%d total magnetic charge density\n", k); k++;
		fprintf(f, "%d number of magnetic monopoles\n", k); k++;
	#endif

	fclose(f);
}


/* Measure observables and write them to file.
* Argument file is typically p->resultsfile, but is anyway specified here in case
* we want to use the same function for writing different files (e.g. for realtime trajectories)
*/
void measure(FILE* file, fields const* f, params const* p, weight* w) {

	double start, end, time;

	// observables that we want to measure
	double action = 0.0;
	double wilson = 0.0;
	double u1wilson = 0.0;

	// Higgs doublet:
	double hopping_phi = 0.0;
	double phi2 = 0.0;
	double phi4 = 0.0;

	// Triplet
	double phi2Sigma2 = 0.0;
	double hopping_Sigma = 0.0;
	double Sigma2 = 0.0;
	double Sigma4 = 0.0;
	double mag_charge = 0.0;
	double mag_charge_abs = 0.0;


	double mod = 0.0;

	// some overlap here. action_local() already calculates local wilson action, hopping terms etc.
	for (long i=0; i<p->sites; i++) {
		action += action_local(f, p, i);
		wilson += local_su2wilson(f, p, i);
		#ifdef U1
			u1wilson += local_u1wilson(f, p, i);
		#endif

		#ifdef HIGGS
			for (int dir=0; dir<p->dim; dir++) {
				hopping_phi += hopping_doublet_forward(f, p, i, dir);
			}
			mod = doubletsq(f->su2doublet[i]);
			phi2 += mod;
			phi4 += mod*mod;
			#ifdef TRIPLET
				phi2Sigma2 += mod * tripletsq(f->su2triplet[i]);
			#endif
		#endif

		#ifdef TRIPLET
			double tripletmod = tripletsq(f->su2triplet[i]);
			Sigma2 += tripletmod;
			Sigma4 += tripletmod * tripletmod;
			for (int dir=0; dir<p->dim; dir++) {
				hopping_Sigma += hopping_triplet_forward(f, p, i, dir);
			}

			// calculate charge density of magnetic monopoles
			double charge = magcharge_cube(p, f, i);
			mag_charge += charge;
			mag_charge_abs += fabs(charge);
		#endif
	}

	// combine results from all nodes.
	start = clock();

	// multicanonical weight and order parameter value
	double weight = 0.0;
	double muca_param = 0.0;
	if (p->multicanonical) {
		/* if w->do_acceptance == 0, order param is likely not up to date
		* because it is not used/updated in field update sweeps. Recalculate it here.
		*/
		if (!w->do_acceptance) {
			calc_orderparam(p, f, w, EVEN); // updates EVEN contribution only
			muca_param = calc_orderparam(p, f, w, ODD); // updates ODD and returns the full value
		} else {
			muca_param = w->param_value[EVEN] + w->param_value[ODD];
		}
		weight = get_weight(w, muca_param);
	}
	// our muca action is S' = S + W, but Kari uses S = S - W.
	// store weight with a minus sign here to ensure compability with Kari's tools
	weight = -1.0 * weight;

	// collect from other nodes
	action = reduce_sum(action);
	wilson = reduce_sum(wilson);
	#ifdef U1
	u1 = reduce_sum(u1wilson);
	#endif
	#ifdef HIGGS
	hopping_phi = reduce_sum(hopping_phi);
	phi2 = reduce_sum(phi2);
	phi4 = reduce_sum(phi4);
	#endif
	#ifdef TRIPLET
	hopping_Sigma = reduce_sum(hopping_Sigma);
	Sigma2 = reduce_sum(Sigma2);
	Sigma4 = reduce_sum(Sigma4);
	mag_charge = reduce_sum(mag_charge);
	mag_charge_abs = reduce_sum(mag_charge_abs);
	// magnetic charge should be quantized in units of 4pi/g
	mag_charge_abs /= (2.0*M_PI*sqrt(p->betasu2)); // this is now the total number of monopoles
		#ifdef HIGGS
		phi2Sigma2 = reduce_sum(phi2Sigma2);
		#endif
	#endif


	#ifdef GRADFLOW
		// for debugging
		Global_current_action = action;
	#endif

	end = clock();
	time = ((double) (end - start)) / CLOCKS_PER_SEC;
	Global_comms_time += time;

	// write to the file from root node. This is very fast performance wise
	// the ordering here should be the same as in print_labels()
	if (!p->rank) {
		fprintf(file, "%g %g ", weight, muca_param);
		fprintf(file, "%g %g ",
			action, wilson/((double)p->vol) );
		#ifdef HIGGS
			fprintf(file, "%g %g %g ",
				hopping_phi/((double)p->vol), phi2/((double)p->vol), phi4/((double)p->vol)
			);
		#endif
		#ifdef TRIPLET
			fprintf(file, "%g %g %g ",
				hopping_Sigma/((double)p->vol), Sigma2/((double)p->vol), Sigma4/((double)p->vol)
			);
			#ifdef HIGGS
			fprintf(file, "%g ",
				phi2Sigma2/((double)p->vol)
			);
			#endif
		#endif
		#ifdef U1
			fprintf(file, "%g ", u1/((double)p->vol) );
		#endif
		#ifdef TRIPLET
			// store total magnetic charge density (should be ~0)
			// and number of monopoles + antimonopoles (should be an integer)
			fprintf(file, "%g %g ", mag_charge, mag_charge_abs );
		#endif
		fprintf(file, "\n");
		fflush(file);
	}

}


/* Calculate local action for the system at site i.
*	The construction is so that a loop over i gives the total action.
*/
double action_local(fields const* f, params const* p, long i) {

	double tot = 0.0;
	tot += local_su2wilson(f, p, i);

	#ifdef U1
		tot += local_u1wilson(f, p, i);
	#endif

	#ifdef HIGGS
		tot += covariant_doublet(f, p, i) + higgspotential(f, p, i);
	#endif
	#ifdef TRIPLET
		double mod = tripletsq(f->su2triplet[i]);
		// potential + covariant derivative. Higgs portal term is included in higgspotential().
		tot += p->msq_triplet * mod + p->b4 * mod * mod;
		tot += covariant_triplet(f, p, i);
	#endif

	return tot;
}


/* Make a label file for local measurements.
*/
void print_labels_local(params* p, char* fname) {
	FILE* f = fopen(fname, "w+");

	int k = 1;
	// first columns are the site coordinates
	for (int dir=0; dir<p->dim; dir++) {
		fprintf(f, "%d x%d\n", k, dir); k++;
	}
	#ifdef TRIPLET
		fprintf(f, "%d Sigma^2\n", k); k++;
		fprintf(f, "%d magnetic charge (integer)\n", k); k++;
	#endif

	fclose(f);
}

/* Measures and writes quantities locally at each site. Extensive!
* Output is not in any particular order in terms of the coordinates.
* Assumes that sites in each node are ordered in the same fashion
* and that all nodes have the same number of (real) sites!
* Could also just send the coordinates, but this needs more comms and/or memory...
*/
void measure_local(char* fname, params const* p, fields const* f) {

	FILE* file;

	int n_meas = 1; // includes offset!
	#ifdef TRIPLET
		// Tr Sigma^2 = 0.5*Sigma^a Sigma^a
		double Sigma2[p->sites]; n_meas++;
		double magcharge[p->sites]; n_meas++;
	#endif

	for (long i=0; i<p->sites; i++) {

		#ifdef TRIPLET

			Sigma2[i] = tripletsq(f->su2triplet[i]);
			magcharge[i] = magcharge_cube(p, f, i) / (2.0*M_PI*sqrt(p->betasu2)); // integer!
		#endif
	} // end i

	if (!p->rank) {
		file = fopen(fname, "wb");
	}

	int coord_offset[p->dim];
	for (int dir=0; dir<p->dim; dir++) {
		coord_offset[dir] = p->offset[dir];
	}

	#ifdef MPI
	/* Send everything to root node for writing.
	* There will be one send per measurement array, incl. offset.
	* So need tags to avoid errors, but I use blocking sends. */
	int tag = 0;
	if (p->rank != 0) {
		MPI_Send(&coord_offset, p->dim, MPI_INT, 0, tag, MPI_COMM_WORLD); tag++;
		#ifdef TRIPLET
			MPI_Send(&Sigma2, p->sites, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD); tag++;
			MPI_Send(&magcharge, p->sites, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD); tag++;
		#endif
	}

	#endif



	if (!p->rank) {
		// loop over MPI ranks
		for (int r=0; r<p->size; r++) {

			#ifdef MPI
			// get the data and offsets from rank == r node (if r=0, just write own meas)
			if (r != 0) {
				tag = 0;
				MPI_Recv(&coord_offset, p->dim, MPI_INT, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
				#ifdef TRIPLET
					MPI_Recv(&Sigma2, p->sites, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
					MPI_Recv(&magcharge, p->sites, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
				#endif
			}
			#endif

			for (long i=0; i<p->sites; i++) {


				// Binary: smaller files
				int coords[p->dim];
				for (int dir=0; dir<p->dim; dir++) {
					coords[dir] = p->coords[i][dir] + coord_offset[dir];
				}

				fwrite(coords, sizeof(*coords), p->dim, file);

				// merge measurements at site=i to an array for writing
				int k = 0;
				double meas[n_meas - 1];

				#ifdef TRIPLET
					meas[k] = Sigma2[i]; k++;
					meas[k] = magcharge[i]; k++;
				#endif
				fwrite(meas, sizeof(*meas), n_meas - 1, file);

				fwrite("\n", sizeof(char), 1, file);
				// end binary

				// Formatted: text files only
				/*
				// write coords
				for (int dir=0; dir<p->dim; dir++) {
					fprintf(file, "%ld ", p->coords[i][dir] + coord_offset[dir]);
				}
				// then the measurements
				#ifdef TRIPLET
					fprintf(file, "%g %g ", Sigma2[i], magcharge[i]);
				#endif

				fprintf(file, "\n");
				*/

			} // end i

		} // end r

	}

	if (!p->rank) {
		fclose(file);
	}

	/*
	#ifdef MPI // wait for nonblocking sends to be received
		for (int k=0; k<n_meas; k++) {
			MPI_Wait(&req[k], MPI_STATUS_IGNORE);
		}
		barrier();
	#endif
	*/
	barrier();

}

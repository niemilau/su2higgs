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
	fprintf(f, "%d SU(2) Wilson (divided by beta)\n", k); k++;
	#ifdef U1
		fprintf(f, "%d U(1) Wilson (divided by beta)\n", k); k++;
	#endif
	#if (NHIGGS > 0)
		fprintf(f, "%d 2*(phi^2 - phi^+(x+i) U_i(x) phi(x)) (avg over directions)\n", k); k++;
		fprintf(f, "%d phi^2\n", k); k++;
		fprintf(f, "%d phi^4\n", k); k++;
	#endif
	#if (NHIGGS == 2)
		fprintf(f, "%d 2*(phi^2 - phi^+(x+i) U_i(x) phi(x)) for phi_2 (avg over directions)\n", k); k++;
		fprintf(f, "%d phi2^2\n", k); k++;
		fprintf(f, "%d phi2^4\n", k); k++;
		fprintf(f, "%d R = Re phi1^+ phi2\n", k); k++;
		fprintf(f, "%d I = Im phi1^+ phi2\n", k); k++;
	#endif
	#ifdef TRIPLET
		fprintf(f, "%d hopping_Sigma (avg over directions)\n", k); k++;
		fprintf(f, "%d Sigma^2\n", k); k++;
		fprintf(f, "%d Sigma^4\n", k); k++;
	#endif
	#if (NHIGGS > 0) && defined TRIPLET
		fprintf(f, "%d phi^2 Sigma^2\n", k); k++;
	#endif
	#ifdef TRIPLET
		fprintf(f, "%d total magnetic charge density\n", k); k++;
		fprintf(f, "%d number of magnetic monopoles\n", k); k++;
	#endif
	#ifdef SINGLET
		fprintf(f, "%d S^2 - S(x)S(x+i) (avg over directions)\n", k); k++;
		fprintf(f, "%d S\n", k); k++;
		fprintf(f, "%d S^2\n", k); k++;
		fprintf(f, "%d S^3\n", k); k++;
		fprintf(f, "%d S^4\n", k); k++;
		#if (NHIGGS==1)
			fprintf(f, "%d S phi^2\n", k); k++;
			fprintf(f, "%d S^2 phi^2\n", k); k++;
		#endif
	#endif


	fclose(f);
}


/* Measure observables and write them to file.
* Argument file is typically l->resultsfile, but is anyway specified here in case
* we want to use the same function for writing different files (e.g. for realtime trajectories)
*/
void measure(FILE* file, lattice const* l, fields const* f, params const* p, weight* w) {

	double start, end, time;

	// observables that we want to measure
	double action = 0.0;
	double wilson = 0.0;
	double u1wilson = 0.0;

	#if (NHIGGS > 0)
		// Higgs doublets:
		double covariant_phi[NHIGGS] = {0.0};
		double phi2[NHIGGS] = {0.0};
		double phi4[NHIGGS] = {0.0};
	#endif

	#if (NHIGGS == 2)
		complex phi12;
		phi12.re = 0.0; phi12.im = 0.0;
	#endif

	// Triplet
	double phi2Sigma2 = 0.0;
	double hopping_Sigma = 0.0;
	double Sigma2 = 0.0;
	double Sigma4 = 0.0;
	double mag_charge = 0.0;
	double mag_charge_abs = 0.0;

	// singlet
	double singlet_covariant = 0.0;
	double singlet = 0.0;
	double singlet2 = 0.0;
	double singlet3 = 0.0;
	double singlet4 = 0.0;
	double Sphisq = 0.0;
	double S2phisq = 0.0;

	// some overlap here. action_local() already calculates local wilson action, hopping terms etc.
	for (long i=0; i<l->sites; i++) {
		action += action_local(l, f, p, i);
		wilson += local_su2wilson(l, f, p, i);
		#ifdef U1
			u1wilson += local_u1wilson(l, f, p, i);
		#endif

		double mod = 0.0;

		#if (NHIGGS > 0)

			for (int db=0; db<NHIGGS; db++) {

				mod = doubletsq(f->su2doublet[db][i]);
				// Covariant derivative, includes all directions
				covariant_phi[db] += covariant_doublet(l, f, i, db);
				phi2[db] += mod;
				phi4[db] += mod*mod;
			}

			#ifdef TRIPLET // only implemented for one Higgs!!
				mod = doubletsq(f->su2doublet[0][i]);
				phi2Sigma2 += mod * tripletsq(f->su2triplet[i]);
			#endif

		#endif

		// some specific operators for 2 Higgs potential
		#if (NHIGGS == 2)
			complex f12 = get_phi12(f->su2doublet[0][i], f->su2doublet[1][i]);
			phi12.re += f12.re;
			phi12.im += f12.im;
		#endif

		#ifdef TRIPLET
			double tripletmod = tripletsq(f->su2triplet[i]);
			Sigma2 += tripletmod;
			Sigma4 += tripletmod * tripletmod;
			for (int dir=0; dir<l->dim; dir++) {
				hopping_Sigma += hopping_triplet_forward(l, f, p, i, dir);
			}

			// calculate charge density of magnetic monopoles
			double charge = magcharge_cube(l, f, p, i);
			mag_charge += charge;
			mag_charge_abs += fabs(charge);
		#endif

		#ifdef SINGLET
			double S = f->singlet[i][0];
			for (int dir=0; dir<l->dim; dir++) {
				double S_next = f->singlet[ l->next[i][dir] ][0]; 
				singlet_covariant += S*S - S*S_next;
			}
			singlet += S;
			singlet2 += S*S;
			singlet3 += S*S*S;
			singlet4 += S*S*S*S;
			#if (NHIGGS == 1)
				mod = doubletsq(f->su2doublet[0][i]);
				Sphisq += mod * S;
				S2phisq += mod * S * S;
			#endif
		#endif
	}

	// combine results from all nodes.
	start = clock();

	// multicanonical weight and order parameter value
	double weight = 0.0;
	double muca_param = 0.0;
	if (p->multicanonical) {
		muca_param = w->param_value[EVEN] + w->param_value[ODD];
		weight = get_weight(w, muca_param);
	}
	// our muca action is S' = S + W, but Kari uses S = S - W.
	// store weight with a minus sign here to ensure compability with Kari's tools
	weight = -1.0 * weight;

	// collect from other nodes
	action = reduce_sum(action, l->comm);
	wilson = reduce_sum(wilson, l->comm);
	#ifdef U1
	u1wilson = reduce_sum(u1wilson, l->comm);
	#endif

	#if (NHIGGS > 0)
	for (int db=0; db<NHIGGS; db++) {
		covariant_phi[db] = reduce_sum(covariant_phi[db], l->comm);
		phi2[db] = reduce_sum(phi2[db], l->comm);
		phi4[db] = reduce_sum(phi4[db], l->comm);
	}
	#endif

	#if (NHIGGS == 2)
		phi12.re = reduce_sum(phi12.re, l->comm);
		phi12.im = reduce_sum(phi12.im, l->comm);
	#endif

	#ifdef TRIPLET
	hopping_Sigma = reduce_sum(hopping_Sigma, l->comm);
	Sigma2 = reduce_sum(Sigma2, l->comm);
	Sigma4 = reduce_sum(Sigma4, l->comm);
	mag_charge = reduce_sum(mag_charge, l->comm);
	mag_charge_abs = reduce_sum(mag_charge_abs, l->comm);
	// magnetic charge should be quantized in units of 4pi/g
	mag_charge_abs /= (2.0*M_PI*sqrt(p->betasu2)); // this is now the total number of monopoles
		#if (NHIGGS > 0)
		phi2Sigma2 = reduce_sum(phi2Sigma2, l->comm);
		#endif
	#endif

	#ifdef SINGLET
		singlet_covariant = reduce_sum(singlet_covariant, l->comm);
		singlet = reduce_sum(singlet, l->comm);
		singlet2 = reduce_sum(singlet2, l->comm);
		singlet3 = reduce_sum(singlet3, l->comm);
		singlet4 = reduce_sum(singlet4, l->comm);
		#if (NHIGGS==1)
			Sphisq = reduce_sum(Sphisq, l->comm);
			S2phisq = reduce_sum(S2phisq, l->comm);
		#endif
	#endif


	#ifdef GRADFLOW
		// for debugging
		Global_current_action = action;
	#endif

	end = clock();
	time = ((double) (end - start)) / CLOCKS_PER_SEC;
	Global_comms_time += time;

	/* write volume averages to file. For hopping terms and gauge actions,
	* take also the average over directions */

	double vol = ((double) l->vol);

	// Write to the file from root node. Use same ordering here as in print_labels()
	// For Wilson gauge actions I divide by their respective beta, which was included in local_su2wilson() etc
	if (!l->rank) {
		fprintf(file, "%g %g ", weight, muca_param);
		fprintf(file, "%g %g ",
			action, wilson / (vol * p->betasu2)
		);
		#ifdef U1
			fprintf(file, "%g ", u1wilson / (vol * p->betau1));
		#endif

		#if (NHIGGS > 0 )
		for (int db=0; db<NHIGGS; db++) {
			fprintf(file, "%g %g %g ",
				covariant_phi[db]/(vol * l->dim), phi2[db]/vol, phi4[db]/vol
			);
		}
		#endif

		#if (NHIGGS == 2)
		fprintf(file, "%g %g ",
			phi12.re/vol, phi12.im/vol
		);
		#endif

		#ifdef TRIPLET
			fprintf(file, "%g %g %g ",
				hopping_Sigma/(vol * l->dim), Sigma2/vol, Sigma4/vol
			);
			#if (NHIGGS > 0)
			fprintf(file, "%g ",
				phi2Sigma2/vol
			);
			#endif
		#endif
		#ifdef TRIPLET
			// store total magnetic charge density (should be ~0)
			// and number of monopoles + antimonopoles (should be an integer)
			fprintf(file, "%g %g ", mag_charge, mag_charge_abs );
		#endif

		#ifdef SINGLET
			fprintf(file, "%g ",
				singlet_covariant / (vol * l->dim)
			);
			fprintf(file, "%g %g %g %g ",
				singlet/vol, singlet2/vol, singlet3/vol, singlet4/vol
			);
			#if (NHIGGS == 1)
				fprintf(file, "%g %g ",
					Sphisq/vol, S2phisq/vol
				);
			#endif
		#endif

		fprintf(file, "\n"); // end write
		fflush(file);
	}

}


/* Calculate local action for the system at site i.
*	The construction is so that a loop over i gives the total action. */
double action_local(lattice const* l, fields const* f, params const* p, long i) {

	double tot = 0.0;
	tot += local_su2wilson(l, f, p, i);

	#ifdef U1
		tot += local_u1wilson(l, f, p, i);
	#endif

	tot += higgspotential(f, p, i); // does nothing if no scalars are present
	#if (NHIGGS > 0)
		for (int db=0; db<NHIGGS; db++) tot += covariant_doublet(l, f, i, db);
	#endif

	#ifdef TRIPLET
		tot += covariant_triplet(l, f, p, i);
	#endif

	#ifdef SINGLET
		// add just the kinetic term, rest is in higgspotential()
		double S = f->singlet[i][0];
		tot += l->dim * S*S;
		for (int dir=0; dir<l->dim; dir++) {
			long next = l->next[i][dir];
			tot -= S * f->singlet[next][0];
		}
	#endif

	return tot;
}


/* Make a label file for local measurements. */
void print_labels_local(lattice const* l, char* fname) {
	FILE* f = fopen(fname, "w+");

	int k = 1;
	// first columns are the site coordinates
	for (int dir=0; dir<l->dim; dir++) {
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
* Could also just send the coordinates, but this needs more comms and/or memory... */
void measure_local(char* fname, lattice const* l, fields const* f, params const* p) {

	FILE* file;

	int n_meas = 1; // includes offset!
	#ifdef TRIPLET
		// Tr Sigma^2 = 0.5*Sigma^a Sigma^a
		double Sigma2[l->sites]; n_meas++;
		double magcharge[l->sites]; n_meas++;
	#endif

	for (long i=0; i<l->sites; i++) {

		#ifdef TRIPLET

			Sigma2[i] = tripletsq(f->su2triplet[i]);
			magcharge[i] = magcharge_cube(l, f, p, i) / (2.0*M_PI*sqrt(p->betasu2)); // integer!
		#endif
	} // end i

	if (!l->rank) {
		file = fopen(fname, "wb");
	}

	int coord_offset[l->dim];
	for (int dir=0; dir<l->dim; dir++) {
		coord_offset[dir] = l->offset[dir];
	}

	#ifdef MPI
	/* Send everything to root node for writing.
	* There will be one send per measurement array, incl. offset.
	* So need tags to avoid errors, but I use blocking sends. */
	int tag = 0;
	if (l->rank != 0) {
		MPI_Send(&coord_offset, l->dim, MPI_INT, 0, tag, MPI_COMM_WORLD); tag++;
		#ifdef TRIPLET
			MPI_Send(&Sigma2, l->sites, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD); tag++;
			MPI_Send(&magcharge, l->sites, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD); tag++;
		#endif
	}

	#endif



	if (!l->rank) {
		// loop over MPI ranks
		for (int r=0; r<l->size; r++) {

			#ifdef MPI
			// get the data and offsets from rank == r node (if r=0, just write own meas)
			if (r != 0) {
				tag = 0;
				MPI_Recv(&coord_offset, l->dim, MPI_INT, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
				#ifdef TRIPLET
					MPI_Recv(&Sigma2, l->sites, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
					MPI_Recv(&magcharge, l->sites, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); tag++;
				#endif
			}
			#endif

			for (long i=0; i<l->sites; i++) {


				// Binary: smaller files
				int coords[l->dim];
				for (int dir=0; dir<l->dim; dir++) {
					coords[dir] = l->coords[i][dir] + coord_offset[dir];
				}

				fwrite(coords, sizeof(*coords), l->dim, file);

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
				for (int dir=0; dir<l->dim; dir++) {
					fprintf(file, "%ld ", l->coords[i][dir] + coord_offset[dir]);
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

	if (!l->rank) {
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

}

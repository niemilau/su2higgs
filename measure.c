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

	fclose(f);
}


/* Measure observables and write them to file.
*
*/
void measure(fields f, params p, counters* c, weight* w) {

	double start, end, time;

	// observables that we want to measure, plus variables for their MPI reduced versions

	double action = 0.0, action_tot;
	double wilson = 0.0, wilson_tot;
	double u1wilson = 0.0, u1_tot;

	// Higgs doublet:
	double hopping_phi = 0.0, hopping_phi_tot;
	double phi2 = 0.0, phi2_tot;
	double phi4 = 0.0, phi4_tot;

	// Triplet
	double phi2Sigma2 = 0.0, phi2Sigma2_tot;
	double hopping_Sigma = 0.0, hopping_Sigma_tot;
	double Sigma2 = 0.0, Sigma2_tot;
	double Sigma4 = 0.0, Sigma4_tot;



	double mod = 0.0;

	// some overlap here. action_local() already calculates local wilson action, hopping terms etc.
	for (long i=0; i<p.sites; i++) {
		action += action_local(f, p, i);
		wilson += local_su2wilson(f, p, i);
		#ifdef U1
			u1wilson += local_u1wilson(f, p, i);
		#endif

		#ifdef HIGGS
			for (int dir=0; dir<p.dim; dir++) {
				hopping_phi += hopping_doublet_forward(f, p, i, dir);
			}
			mod = doubletsq(f.su2doublet[i]);
			phi2 += mod;
			phi4 += mod*mod;
			#ifdef TRIPLET
				phi2Sigma2 += mod * tripletsq(f.su2triplet[i]);
			#endif
		#endif

		#ifdef TRIPLET
			double tripletmod = tripletsq(f.su2triplet[i]);
			Sigma2 += tripletmod;
			Sigma4 += tripletmod * tripletmod;
			for (int dir=0; dir<p.dim; dir++) {
				hopping_Sigma += hopping_triplet_forward(f, p, i, dir);
			}

			// TEMPORARY: monopoles
			magcharge_cube(&p, &f, i);
		#endif
	}

	// combine results from all nodes.
	start = clock();

	double weight = 0.0; // if not multicanonical, just use zero weight
	if (p.multicanonical) {
		double muca_param = w->param_value[EVEN] + w->param_value[ODD];
		weight = get_weight(*w, muca_param);
	}
	// our muca action is S' = S + W, but Kari uses S = S - W.
	// store weight with a minus sign here to ensure compability with Kari's tools
	weight = -1.0 * weight;

	action_tot = reduce_sum(action);
	wilson_tot = reduce_sum(wilson);
	#ifdef U1
	u1_tot = reduce_sum(u1wilson);
	#endif
	#ifdef HIGGS
	hopping_phi_tot = reduce_sum(hopping_phi);
	phi2_tot = reduce_sum(phi2);
	phi4_tot = reduce_sum(phi4);
	#endif
	#ifdef TRIPLET
	hopping_Sigma_tot = reduce_sum(hopping_Sigma);
	Sigma2_tot = reduce_sum(Sigma2);
	Sigma4_tot = reduce_sum(Sigma4);
		#ifdef HIGGS
		phi2Sigma2_tot = reduce_sum(phi2Sigma2);
		#endif
	#endif


	end = clock();
	time = ((double) (end - start)) / CLOCKS_PER_SEC;
	c->comms_time += time;

	// write to resultsfile. This is very fast performance wise
	// the ordering here should be the same as in print_labels()
	if (!p.rank) {
		fprintf(p.resultsfile, "%g ", weight);
		fprintf(p.resultsfile, "%g %g ",
			action_tot, wilson_tot/((double)p.vol) );
		#ifdef HIGGS
			fprintf(p.resultsfile, "%g %g %g ",
				hopping_phi_tot/((double)p.vol), phi2_tot/((double)p.vol), phi4_tot/((double)p.vol)
			);
		#endif
		#ifdef TRIPLET
			fprintf(p.resultsfile, "%g %g %g ",
				hopping_Sigma_tot/((double)p.vol), Sigma2_tot/((double)p.vol), Sigma4_tot/((double)p.vol)
			);
			#ifdef HIGGS
			fprintf(p.resultsfile, "%g ",
				phi2Sigma2_tot/((double)p.vol)
			);
			#endif
		#endif
		#ifdef U1
			fprintf(p.resultsfile, "%g ", u1_tot/((double)p.vol) );
		#endif
		fprintf(p.resultsfile, "\n");
		fflush(p.resultsfile);
	}

}


/* Calculate local action for the system at site i.
*	The construction is so that a loop over i gives the total action.
*/
double action_local(fields f, params p, long i) {

	double tot = 0.0;
	tot += local_su2wilson(f, p, i);

	#ifdef U1
		tot += local_u1wilson(f, p, i);
	#endif

	#ifdef HIGGS
		tot += covariant_doublet(f, p, i) + higgspotential(f, p, i);
	#endif
	#ifdef TRIPLET
		double mod = tripletsq(f.su2triplet[i]);
		// potential + covariant derivative. Higgs portal term is included in higgspotential().
		tot += p.msq_triplet * mod + p.b4 * mod * mod;
		tot += covariant_triplet(f, p, i);
	#endif

	return tot;
}

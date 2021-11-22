/** @file multicanonical.c
*
* Routines for multicanonical simulations.
*
* Weight is stored in struct weight, w.pos[i] being the "starting" value of
* i. bin and w.W[i] being the value of the weight function at w.pos[i].
* Inside the bin, a linear approximation is used to obtain a continuous weight function.
* The sign of W is chosen so that our weight function INCREASES when we want to suppress
* a particular configuration (ensemble ~ exp(-S - W).
*
* Summary of the algorithm:
*
* 	1. Do a standard checkerboard sweep, updating the multicanonical field with
*			 local updates WITHOUT caring about the weight at all. This will move the system
*			 towards a local minimum.
*		2. After the sweep, recalculate the order parameter and do a global accept/reject
*			 based on the weight change. If rejected, ALL updates in the sweep are undone.
*			 This will create a bias towards non-minimum configurations.
*
*	Additional steps for calculating the weight function, if w.mode != READONLY (== 0)
*
*		2. Each time a global muca acc/rej is performed, find the current bin and accumulate
*			 w.hits in that bin AND its nearby bins (even if update is rejected!). See muca_accumulate_hits().
*		3. After certain number of multicanonical acc/rej steps (according to Kari autocorrelations
*			 don't matter that much here), increase weight in each bin by value proportional to w.hits[bin],
*			 but relative to the weight in first bin which is kept constant (see update_weight()).
*
*		4. Once the simulation has visited both the first and last bin, decrease w.delta to make the
*			 iteration converge. This is necessary because at first the system prefers the canonical
*			 minima, so those obtain large weight very quickly. Afterwards, the mixed-phase configurations are
*		   preferred and the system spends most of the time there, so smaller weight updates are
*			 then required to prevent the weight from becoming flat again.
*
*
***************** On different order parameters ********************
* Most of the routines here for updating the weight and accepting a multicanonical update are
* insensitive to what the actual order parameter is. To implement a new multicanonical order parameter,
* the following steps need to be taken:
*
*		1. #define a constant label for it in a header (su2.h) and add it to get_weight_parameters().
*		2. Tell calc_orderparam(), alloc_muca_arrays(), free_muca_arrays(),
*			 store_muca_fields() and reset_muca_fields() what to do for this case.
*		3. In update_lattice_muca(), call multicanonical_acceptance() after updating all the
*			 fields contributing to the order param and reject the updates if necessary (this is the tricky part).
*/

/* Everywhere in this code, w.bins refers to the number of bins in the strict weighting
* range [w.min, w.max]. However all weight arrays come with w.bins+2 elements:
* the two 'extra' bins have indices 0 and w.bins+1 with w.pos[0] = -infty,
* w.pos[w.bins+1] = w.max (starts at w.max, extends to +infinity). The weight in
* these bins is always w.W[0] = w.W[1], w.W[w.bins+1] = w.W[w.bins], i.e. muca update
* is always accepted if moving outside the weighting range. With this setup
* the weight in the final proper bin w.bins cannot be linearized and constant weight is used there.
* If weight is to be updated, the program may shift w.max slightly so that w.wrk_max does
* not lie in the extra bin starting at w.max; otherwise we would end up updating the extra bin..
*
* Finally, the extra bin w.bins+1 is also written to the weight file.
* Example weight file with 4 bins in range [1.0, 2.0], initial delta=0.3,
* last_max = 0 and weight updating in range [1.0, 1.5]:
*	  4 0.3 0 1.0 1.5
*		1.00 0
*		1.25 0
*		1.50 0
*		1.75 0
*		2.00 0
*/

/* Update 28.7.2020: implemented a slower but 'safer' method for weight calculation
* as described in hep-lat/9804019.
*/

#include "su2.h"


/* Save the current weight into file.
* There will be w.bins + 1 lines, the first one being:
* 	w.bins w.delta w.last_max w.wrk_min w.wrk_max
* the w.bins+1 lines after that read:
* 	w.pos[i] w.W[i]
*
* Note that w.hits is NOT stored! */
void save_weight(lattice const* l, weight const* w) {

	// for read only run, do nothing
	if (w->mode == READONLY)
		return;

	if(!l->rank) {
		FILE *wfile = fopen(w->weightfile, "w");

		fprintf(wfile, "%d %.12lf %d %lf %lf\n", w->bins, w->delta, w->last_max, w->wrk_min, w->wrk_max);

		for(long i=1; i<w->bins+2; i++) {
			fprintf(wfile, "%lf %lf\n", w->pos[i], w->W[i]);
		}

		fclose(wfile);
	}

 if (w->mode == SLOW && !l->rank) {
	 // store arrays needed for 'slow' updating in a separate file
	 FILE *file = fopen("weight_params", "w");
	 fprintf(file, "%d \n", w->bins);

	 // extra bins not written
	 for(long i=1; i<w->bins+1; i++) {
		 fprintf(file, "%lf %ld %ld\n", w->pos[i], w->gsum[i], w->nsum[i]);
	 }
	 fclose(file);
 }

}

/* Load multicanonical weight from weightfile.
* If file does not exist, initializes a new flat weight.
* Assumes that bins, min and max have been read from config file,
* but overrides these if an existing weight file is found.
* DO NOT call this more than once per run. */
void load_weight(lattice const* l, weight *w) {

	printf0("\nLoading weightfile: %s\n", w->weightfile);

	// alloc arrays, accounting for the two 'virtual' bins
  w->pos = malloc( (w->bins+2)* sizeof(*(w->pos)));
  w->W = malloc((w->bins+2) * sizeof(*(w->W)));
	w->hits = malloc((w->bins+2) * sizeof(*(w->hits)));
	w->slope = malloc((w->bins+2) * sizeof(*(w->slope)));
	w->b = malloc((w->bins+2) * sizeof(*(w->b)));

	// how often to update the weight?
	if (w->mode != SLOW) w->update_interval = 100;
	else w->update_interval = 5000;

	w->do_acceptance = 1;

	if (w->mode == SLOW) {
		// slow but safer weight calculation, affects only runs with mode = SLOW
		w->hgram = malloc((w->bins+2) * sizeof(*w->hgram));
		w->gsum = malloc((w->bins+2) * sizeof(*w->gsum));
		w->nsum = malloc((w->bins+2) * sizeof(*w->nsum));
	}

  long i;

	int read_ok = access(w->weightfile, R_OK); // 0 if OK, -1 otherwise
  if( read_ok != 0) {

		// no weight found; initialize a flat weight according to config file
    printf0("Unable to access weightfile!!\n");
		if (w->mode == READONLY) {
			printf0("No multicanonical weight given for a read-only run! Exiting...\n");
			die(20);
		}

		if (w->bins <= 0) {
			printf0("\n!!!!! Multicanonical error: bins <= 0\n");
			die(-1141);
		}
		// create equally sized bins. last (extra) bin in weight file STARTS at w.max:
		double dbin = (w->max - w->min)/((double) (w->bins));

		// initialize all bins except for bin=0
    for(i=0; i<w->bins+1; i++) {
      w->pos[i+1] = w->min + ((double) i) * dbin;
      w->W[i+1] = 0.0;
    }
		w->last_max = 0; // assume starting from min (but see init_last_max())

		/* Assume that weight update range is the same as weighting range, but to avoid
		* complications shift w.max just a little bit, so that w.wrk_max is never inside
		* the last extra bin. */
		w->wrk_max = w->max;
		w->wrk_min = w->min;
		w->max += 0.001 * dbin;
		w->pos[w->bins + 1] = w->max;

		printf0("Initialized new weight \n");

  } else {
		// found existing weight, use it instead of the one specified in config
		FILE *wfile = fopen(w->weightfile, "r");

		// read = how many elements read per in a line
		int read;
		int bins_read;

		// first line: bins, delta, last_max and range for weight updating
		read = fscanf(wfile, "%d %lf %d %lf %lf", &bins_read, &w->delta,
			&w->last_max, &w->wrk_min, &w->wrk_max);

		if (read != 5) {
			printf0("Error reading first line of weightfile! \n");
			die(22);
		}

		// override initial options with the values read from the actual weightfile
		w->bins = bins_read;
		if (w->bins <= 0) {
			printf0("\n!!!!! Multicanonical error: bins <= 0\n");
			die(-1141);
		}

		// realloc accounting for the 2 virtual bins
		w->pos = realloc(w->pos, (w->bins+2) * sizeof(*(w->pos)));
		w->W = realloc(w->W, (w->bins+2) * sizeof(*(w->W)));
		w->hits = realloc(w->hits, (w->bins+2) * sizeof(*(w->hits)));
		w->slope = realloc(w->slope, (w->bins+2) * sizeof(*(w->slope)));
		w->b = realloc(w->b, (w->bins+2) * sizeof(*(w->b)));

		if (w->mode == SLOW) {
			w->hgram = realloc(w->hgram, (w->bins+2) * sizeof(*(w->hgram)));
			w->gsum = realloc(w->gsum, (w->bins+2) * sizeof(*(w->gsum)));
			w->nsum = realloc(w->nsum, (w->bins+2) * sizeof(*(w->nsum)));
		}

		// read bins
    for(i=1; i<w->bins+2; i++) {
      read = fscanf(wfile, "%lf %lf", &(w->pos[i]), &(w->W[i]));
      if(read != 2) {
				printf0("Error reading weightfile! Got %d values at line %ld\n", read, i+1);
				die(22);
      }
    }
		fclose(wfile);

		// sync min and max, although these shouldn't be needed after this routine finishes
		w->min = w->pos[1];
		w->max = w->pos[w->bins+1];

		if (!w->mode != READONLY && (w->wrk_min < w->min || w->wrk_max > w->max)) {
			printf0("\n!!! Multicanonical error: weight update range [%lf, %lf] larger than binning range [%lf, %lf] !!!\n\n",
			 		w->wrk_min, w->wrk_max, w->min, w->max);
			die(233);
		}
		/* again, avoid complications by shifting the last bin so that w.wrk_max is
		* never inside the extra bin */
		if (!w->mode != READONLY && fabs(w->wrk_max - w->max) < 1e-10) {
			double dbin = w->max - w->pos[w->bins];
			w->max += 0.001 * dbin;
			w->pos[w->bins+1] = w->max;
		}

  } // end weight loading

	// sync the virtual bins: constant weight in ranges [-infty, w.min], [w.max, infty]
	w->pos[0] = w->min - 1000.0 * (w->max - w->min); // "-infinity"
	w->W[0] = w->W[1];
	w->W[w->bins+1] = w->W[w->bins]; // position set to w.max

	// prepare w->b and w->slope
	linearize_weight(w);

	printf0("Using weight function with %d bins in range %lf, %lf\n", w->bins, w->min, w->max);
	printf0("Global multicanonical check %d times per even/odd update sweep\n", w->checks_per_sweep);

	int mode = w->mode;
	if (mode != READONLY) printf0("Will modify weight in range %lf, %lf\n", w->wrk_min, w->wrk_max);
	if (mode == FAST) printf0("Fast weight update: starting with delta %.12lf, last_max %d\n", w->delta, w->last_max);
	if (mode == SLOW) printf0("Using slow but safe weight updating: parameter file 'weight_params'\n");

	// restart accumulation of muca hits even if existing weight is loaded
	w->muca_count = 0;
	for (i=0; i<w->bins+2; i++) {
		w->hits[i] = 0;
		if (w->mode == SLOW) w->hgram[i] = exp(w->W[i]);
	}

	if (w->mode == SLOW) {
		/* read in work arrays from separate file */
		load_weight_params(l->rank, "weight_params", w);
	}

	save_weight(l, w);

}

/* Read param file required for "slow" weight updating.
* Here rank is just used for printing errors */
void load_weight_params(int rank, char* fname, weight* w) {

	if(access(fname, R_OK) != 0) {
		if (!rank) printf("!!! Unable to access weight parameter file. Initializing new\n");
		for (int i=0; i<w->bins+2; i++) {
			// cannot be 0!!
			w->gsum[i] = 5;
			w->nsum[i] = 1;
		}
	} else {

		FILE *file = fopen(fname, "r");

		int bins_read;
		int read = fscanf(file, "%d", &bins_read);
		if (bins_read != w->bins) {
			// mismatch!
			if (!rank) printf( "Error reading parameter file!! Expected %d bins, got %d\n\n", w->bins, bins_read);
			die(991);
		}

		for(int i=1; i<w->bins+1; i++) {
			double temp;

			read = fscanf(file, "%lf %ld %ld", &temp, &(w->gsum[i]), &(w->nsum[i]));
			if(read != 3) {
				if (!rank) printf("Error reading 'weight_params'!! Got %d values at line %d\n\n", read, i+1);
				die(992);
			}
		}

		fclose(file);

	}

}

// Make arrays 'slope' and 'b' such that W(R) = b + R*slope in each bin
void linearize_weight(weight* w) {
	// extra bins: constant weight
	w->slope[0] = 0; w->slope[w->bins+1] = 0;
	w->b[0] = w->W[0]; w->b[w->bins+1] = w->W[w->bins+1];

	for (int i=1; i<w->bins+1; i++) {
		w->slope[i] = (w->W[i+1] - w->W[i]) / (w->pos[i+1] - w->pos[i]);
		w->b[i] = w->W[i] - w->slope[i] * w->pos[i];
	}
}

/* Get linearized weight corresponding to order
* parameter value val. */
double get_weight(weight const* w, double val) {

	// no need to linearize if outside the binning range
	if (val >= w->max) {
		return w->b[w->bins+1];
	} else if (val < w->min) {
		return w->b[0];
	}

	int bin = whichbin(w, val);
	// should not return 0 or w.bins+1 because these correspond to checks above
	if (bin >= w->bins+1) {
		printf("Error in get_weight (multicanonical.c) !!!\n");
		return w->b[w->bins+1];
	}

	// W(R) = b + slope*R
	return w->b[bin] + w->slope[bin] * val;

	/*
	// linearize the weight function in this bin
	double val_prev = w->pos[bin];
	double nextW, val_next;

	nextW = w->W[bin+1];
	val_next = w->pos[bin+1];

	return w->W[bin] + (val - val_prev) * (nextW - w->W[bin]) / (val_next - val_prev);
	*/
}


/* Global accept/reject step for multicanonical updating.
* oldval is the old order parameter value before field was updated locally
* Return 1 if update was accepted, 0 otherwise.*/
int multicanonical_acceptance(lattice const* l, weight* w, double oldval, double newval) {

	// if we call this function while w->do_acceptance is 0 then something went wrong
	if (!w->do_acceptance) {
		printf0("Should not get here!! in multicanonical.c\n");
		die(-1000);
	}

	// acc/rej only in root node
	int accept;
	if (l->rank == 0) {

		double W_new, W_old;

		W_new = get_weight(w, newval);
		W_old = get_weight(w, oldval);

		double diff = W_new - W_old;

		if(exp(-(diff)) > dran()) {
      accept = 1;
    } else {
      accept = 0;
    }
	}

	// broadcast outcome to all nodes
	bcast_int(&accept, l->comm);

	// update hits and call update_weight() if necessary.
	// do this even if the update was rejected
	if (w->mode != READONLY) {

		if (!accept) {
			newval = oldval;
		}

		muca_accumulate_hits(w, newval);
		w->muca_count++;

		/* update weight if necessary */
		if (w->muca_count % w->update_interval == 0) {

			if (w->mode == FAST) {
				int tunnel = update_weight(w);
				save_weight(l, w);
				if (tunnel) {
					printf0("\nReducing weight update factor! Now %.12lf \n", w->delta);
				}

			} else if (w->mode == SLOW) {
				update_weight_slow(w);
				save_weight(l, w);
			}

			w->muca_count = 0;
		}
	} // end !readonly

	return accept;
}

/* Update w.hits array. This is done in a "Gaussian" fashion:
* if the order parameter has value 'val' corresponding to bin 'i',
* then we update w.hits[i] by a lot but also nearby bins for a smaller amount.
* The logic here is adapted from Kari's susy code (multican_generic.c) */
void muca_accumulate_hits(weight* w, double val) {

	int bin = whichbin(w, val);

	/* Update bin based on condition R_{i-1/2} <= R < R_{i+1/2}, so that the effective
	* bin width is (pos[i+1] - pos[i-1]) / 2 ; R is the order parameter.
	* Exception: first bin (i=1) and the last extra bin. I take these bins to have size d_min = R_2 - R_1
	* and d_max = R_imax - R_{imax-1}, where the weight range ends at R_imax.
	* If R >= R_imax + 0.5*d_max or R < R_1 - 0.5*d_min, do nothing,
	* but if R_{imax-1/2} <= R < R_imax update the last extra bin (sync later)
	* but also its neighbors by a smaller amount */

	if (bin == 0) {
		double dmin = (w->pos[2] - w->pos[1]);
		if (val < w->pos[1] - 0.5*dmin) return;
 		else bin += 1;

	} else if (bin < w->bins+1) {
		double diff1 = val - w->pos[bin];
		double diff2 = w->pos[bin+1] - val;
		if (diff1 > diff2) {
			bin = bin + 1;
		}
	} else {
		// extra bin
		double dmax = w->pos[w->bins+1] - w->pos[w->bins];
		if (val > w->pos[w->bins+1] + 0.5*dmax) return;
	}

	if (w->mode == FAST) {
		// fast weight update
		w->hits[bin] += 5;
		// nearby bins, but don't accumulate hits in extra bins
		if (bin > 1) w->hits[bin-1] += 3;
		if (bin > 2) w->hits[bin-2] += 1;
		if (bin < w->bins) w->hits[bin+1] += 3;
		if (bin < w->bins-1) w->hits[bin+2] += 1;

	} else if (w->mode == SLOW) {
		// slow update, calculate histogram estimate
		w->hits[bin]++;
		double weight = get_weight(w, val);
		w->hgram[bin] += exp(weight); // no minus sign with my weight convention
	}

	/* Note that any hits in bins outside the work range
	* are actually not used by update_weight(), including those in the extra bins */
}

/* Update weight in each bin based on hits in the bin,
* and reset the "accumulation" array w->hits.
* After updating, checks whether
* w->delta should be decreased for the next loop */
int update_weight(weight* w) {

	if (w->mode == READONLY) {
		return 0;
	}

	int firstbin = -1;
	int lastbin = -1;
	double delta = 0.0;

	/* update all bins relative to the first bin in the update range.
	* This means that we don't need to change anything in bins before firstbin */
	for (int i=0; i<w->bins+2; i++) {
		if (w->pos[i] >= w->wrk_min && w->wrk_max >= w->pos[i])  {

			if (firstbin < 0) {
				firstbin = i;
			}

			delta = (w->hits[i] - w->hits[firstbin]) * w->delta / ((double) w->bins);
			w->W[i] += delta;
			lastbin = i;

		} else if (w->wrk_max < w->pos[i]) {
			// bins after wrk_max: increase by the same delta that was used in the last work bin
			w->W[i] += delta;
		}
	}
	// sync extra bins, although this should be automatic in the above loop
	w->W[0] = w->W[1];
	w->W[w->bins+1] = w->W[w->bins];

	linearize_weight(w);

	/* Reduce w->delta if we tunneled from one end of the work range to the other */
	int tunnel = 0;
	if (w->last_max && w->hits[1] > 0) {
		w->last_max = 0;
		tunnel = 1;
	} else if (!w->last_max && w->hits[w->bins] > 0) {
		w->last_max = 1;
		tunnel = 1;
	}

	if (tunnel > 0) {
		w->delta /= 1.5; // seems to work well
	}

	// reset hit accumulation
	for (int i=0; i<w->bins+2; i++) {
		w->hits[i] = 0;
	}
	return tunnel;
}

/* Weight update according to the recipe in hep-lat/9804019
* Remember that in my implementation, bin = 0 extends effectively from -infty to w.min */
void update_weight_slow(weight* w) {

	double w_old[w->bins+2];
	memcpy(w_old, w->W, (w->bins+2) * sizeof(w_old[0]));

	/* Calculate canonical probabilities per unit 'length' (necessary if bins have diff sizes).
	* Here the first and last bins get the same treatment as in accumulate_hits() */
	double dmin = w->pos[2] - w->pos[1];
	//double dmax = w->pos[w->bins+1] - w->pos[w->bins];

	w->hgram[1] /= dmin;
	for (int i=2; i<w->bins+1; i++) w->hgram[i] /= (w->pos[i+1] - w->pos[i-1]);

	// calculate the weight relative to the first bin in the update range
	int firstbin = -1;
	int lastbin = -1;

	// minimum number of 'hits' in a bin before taking it into account
	int hits_min = 8;

	for (int i=1; i<w->bins+1; i++) {
		if (w->pos[i] >= w->wrk_min && w->wrk_max >= w->pos[i]) {

			// do not update the first bin in work range, or bins before that
			if (firstbin < 0) {
				firstbin = i;
				continue;
			}

			w->nsum[i] += w->hits[i];

			int g = 0;
			if (w->hits[i] >= hits_min && w->hits[i-1] >= hits_min) {
				g = w->hits[i] + w->hits[i-1];
			} else {
				// less than minimum number of hits, this run does not contribute
				w->W[i] = w_old[i] - w_old[i-1] + w->W[i-1];
				continue;
			}

			long gsum_new = w->gsum[i] + g; // w->gsum[i] = sum of all g factors from earlier runs
			double ln = w->hgram[i-1] / w->hgram[i];
			ln = -1.0*g*log(ln); // minus sign with my weight convention

			w->W[i] = w->W[i-1] + ((w_old[i]-w_old[i-1]) * w->gsum[i] + ln) / ((double) gsum_new);

			w->gsum[i] = gsum_new;

		} else if (w->wrk_max < w->pos[i]) {
			// outside work range, preserve weight relative to the last bin in work range
			if (lastbin < 0) lastbin = i-1;
			double diff = w->W[lastbin] - w_old[lastbin];
			w->W[i] += diff;
		}

	} // end i

	// reset hits
	for (int i=0; i<w->bins+2; i++) {
		w->hits[i] = 0;
		w->hgram[i] = 0.0;
	}
	// sync extra bins
	w->W[0] = w->W[1];
	w->W[w->bins+1] = w->W[w->bins];

	// overcorrect: make less visited bins more likely (but preserve relative weight outside update range)
	double c = 2.0;
	/* will first modify w->W, then use those to calculate w->b, w->slope,
	* then finally undo the overcorrection to w->W but not to w->b, w->slope
	* => simulation uses the corrected weight but afterwards updates the original w->W */
	memcpy(w_old, w->W, (w->bins+2) * sizeof(w_old[0]));

	double corr = 0.0;

	/* Again, calculate probability per unit 'length'
	* Take dmin = width of first bin in the work range.
	* firstbin does not get overcorrected => overcorrect relative to it.
	* This preserves relative weights also for bins < firstbin   */
	if (firstbin == 1) {
		dmin = w->pos[firstbin+1] - w->pos[firstbin];
	} else {
		dmin = w->pos[firstbin+1] - w->pos[firstbin-1];
	}

	for (int i=firstbin+1; i<w->bins+1; i++) {

		if (i <= lastbin) {

			double diff1 = w->pos[i+1] - w->pos[i-1];
			corr = c * log( (dmin / diff1) *  w->nsum[i] / w->nsum[firstbin]);

			w->W[i] += corr;
		} else {
			// outside update range, add corr calculated in the last update bin
			w->W[i] += corr;
		}
	}

	w->W[0] = w->W[1];
	w->W[w->bins+1] = w->W[w->bins];
	linearize_weight(w);

	// w->slope and w->b done, restore the uncorrected weight
	memcpy(w->W, w_old, (w->bins+2) * sizeof(w_old[0]));
}


/* Return bin index corresponding to given value of order parameter.
* If outside the binning range, returns the extra bin 0 or w.bins+1 */
int whichbin(weight const* w, double val) {

	if (val < w->min) {
		return 0;
	} else if (val >= w->max) {
		return w->bins+1;

	} else {
		// do a quick binary search
		int bin;
		int bmin = 1; int bmax = w->bins;
		double current, next;

		while (bmin <= bmax) {
			bin = bmin + (bmax - bmin) / 2;
			current = w->pos[bin];
			next = w->pos[bin+1];

			if ((val < next && val > current) || (fabs(val - current) < 1e-12) ) {
				return bin;
			}

			if (val > current) {
				bmin = bin + 1;
			} else {
				bmax = bin - 1;
			}
		}

		// end binary search, if we got here then something went wrong!
		printf("\n!!! WARNING: Error in multicanonical.c: failed to find bin for value = %lf!!!\n\n", val);
		return w->bins;
	}
}

/* Calculate muca order parameter and distribute to all nodes.
* Only the contribution from sites with parity = par is recalculated
* while the other parity contribution is read from w.param_value */
double calc_orderparam(lattice const* l, fields const* f, params const* p, weight* w, char par) {
	double tot = 0.0;
	long offset, max;
	if (par == EVEN) {
		offset = 0; max = l->evensites;
	} else {
		offset = l->evensites; max = l->sites;
	}

	switch(w->orderparam) {

#ifdef TRIPLET
		case SIGMASQ :
			for (long i=offset; i<max; i++) tot += tripletsq(f->su2triplet[i]);
			break;
#endif

#if defined (TRIPLET) && (NHIGGS > 0)
    case PHI2MINUSSIGMA2 :
      for (long i=offset; i<max; i++) tot += doubletsq(f->su2doublet[0][i]) - tripletsq(f->su2triplet[i]);
      break;
#endif

#if (NHIGGS > 0)
		case PHISQ :
			for (long i=offset; i<max; i++) tot += doubletsq(f->su2doublet[0][i]);
			break;
#endif

#if (NHIGGS > 1)
		case PHI2SQ :
			for (long i=offset; i<max; i++) tot += doubletsq(f->su2doublet[1][i]);
			break;
#endif

	} // end switch

	tot = allreduce(tot, l->comm) / l->vol;

	w->param_value[par] = tot;
	// add other parity contribution
	return tot + w->param_value[ otherparity(par) ];
}


/* Allocate field backup arrays (no halos to save memory) */
void alloc_muca_backups(lattice const* l, weight* w) {
	switch(w->orderparam) {

#ifdef TRIPLET
		case SIGMASQ :
			w->fbu.su2triplet = make_field(l->sites, SU2TRIP);
			break;

#if (NHIGGS > 0)
    case PHI2MINUSSIGMA2 :
			w->fbu.su2triplet = make_field(l->sites, SU2TRIP);
			w->fbu.su2doublet[0] = make_field(l->sites, SU2DB);
			break;
#endif
#endif

#if (NHIGGS > 0)
		case PHISQ :
			w->fbu.su2doublet[0] = make_field(l->sites, SU2DB);
			break;
#endif
#if (NHIGGS > 1)
		case PHI2SQ :
			w->fbu.su2doublet[1] = make_field(l->sites, SU2DB);
			break;
#endif

	} // end switch
}

// Free memory allocated for multicanonical
void free_muca_arrays(fields* f, weight *w) {
	free(w->pos);
  free(w->W);
	free(w->hits);
	free(w->slope);
	free(w->b);
	if (w->mode == SLOW) {
		free(w->hgram);
		free(w->gsum);
		free(w->nsum);
	}

	// free backup arrays
	switch(w->orderparam) {
#ifdef TRIPLET
		case SIGMASQ :
			free_field(w->fbu.su2triplet);
			break;
#if (NHIGGS > 0)
    case PHI2MINUSSIGMA2 :
			free_field(w->fbu.su2triplet);
			free_field(w->fbu.su2doublet[0]);
			break;
#endif
#endif

#if (NHIGGS > 0)
		case PHISQ :
			free_field(w->fbu.su2doublet[0]);
			break;
#endif
#if (NHIGGS > 1)
		case PHI2SQ :
			free_field(w->fbu.su2doublet[1]);
			break;
#endif
	}
}

/* Initialize last_max to 1 if we start closer to max than min (e.g. after thermalization) */
void init_last_max(weight* w) {
	double val = w->param_value[EVEN] + w->param_value[ODD];
	double diff_min = fabs(val - w->wrk_min);
	double diff_max = fabs(val - w->wrk_max);
	w->last_max = (diff_min >= diff_max);
}

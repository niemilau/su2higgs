/** @file multicanonical.c
*
* Routines for multicanonical simulations.
*
* We assume that a range of order parameter values (min, max) is given
* in the config file, and restrict weighting to this range.
* This range is divided into given number of bins of equal width.
*
* Weight is stored in struct weight, w.pos[i] being the "starting" value of
* i. bin and w.W[i] being the value of the weight function at w.pos[i].
* Inside the bin, we perform a linear approximation to obtain a continuous weight function.
*
*	At least two good options to calculate the weight:
* 	1. Recursive computation described in hep-lat/9804019. Hard to implement, but
*			 will eventually converge to the exact weight function.
*		2. Simpler method similar to the Landau-Wang algorithm, but information from
*			 previous short runs is "forgotten" in each iteration. Kari uses this in his
*			 modern codes. Easy to implement, but may not give the exact weight if the
*			 "delta" parameter becomes too small too early. Usually still sufficient.
*
* Following Kari, our multicanonical action is S = \sum_x L(x) - W(R),
* i.e. weight DECREASES when a configuration is suppressed.
*
*	We follow the parallel method described in hep-lat/9804019. In short, we update sites locally using
* canonical updates, and perform a global accept/reject based on the weight change only after a full
* even or odd sweep (half volume). The weight function is calculated recursively using the recipe of
* hep-lat/9804019, including the overcorrection improvements. NB! Kari measures hits in bin based on
* R_{i - 1/2} <= R <= R_{i + 1/2}, while I consider the range R_{i} <= R < R_{i+1}
*
*
* Kari tips: Behavior of the weight function outside the given range can in principle be chosen arbitrary.
* Usually it's enough to set it to a constant value there, so that the total weight function is still continuous.
* For bubble nucleation it may be better to set an infinite weight there instead.
* We choose "infinite" weight: w.outsideW_min = w.outsideW_max = 1000, but this can be adjusted.
*	Weight in last bin is chosen constant. Note that only relative weights matter!!
* For the weight range, it's best to choose it so that the peaks are mostly included
*
*
*/

#include "su2.h"


/* Load multicanonical weight from weightfile;
* binning range is read earlier in parameters.c.
* If file does not exist, initializes a new flat weight.
* The weight can be saved with save_weight().
* DO NOT call this more than once per run
* (unless memory leaks are your thing).
*
* Original implementation by David Weir.
*/
void load_weight(params p, weight *w) {

	printf0(p, "\nLoading weightfile: %s\n", w->weightfile);

  w->pos = malloc(w->bins * sizeof(*(w->pos)));
  w->W = malloc(w->bins * sizeof(*(w->W)));
	w->hits = malloc(w->bins * sizeof(*(w->hits)));
	w->h = malloc(w->bins * sizeof(*(w->h)));
	w->gsum = malloc(w->bins * sizeof(*(w->gsum)));

	w->dbin = (w->max - w->min)/((double)w->bins);
	w->update_interval = 1000; // Kari used 2000 in MSSM

  long i;
	for (i=0; i<w->bins; i++) {
		w->pos[i] = w->min + ((double) i) * w->dbin;
	}

  if(access(w->weightfile, R_OK) != 0) {

    printf0(p,"Unable to access weightfile!!\n");

		if (w->readonly) {
			printf0(p, "No multicanonical weight given for a read-only run! Exiting...\n");
			die(20);
		}

    for(i=0; i<w->bins; i++) {
      w->W[i] = 0.0;
			w->hits[i] = 0;
			w->h[i] = 0.0;
			w->gsum[i] = 0.0;
    }
		w->last_max = 1; // assume starting from max
		w->m = 0;
		if (w->absolute_bounds) {
			w->outsideW_max = -1000.0;
			w->outsideW_min = -1000.0;
		} else {
			w->outsideW_max = 0.0;
			w->outsideW_min = 0.0;
		}
		printf0(p, "Initialized new weight \n");

  } else {

		// found existing weight, check that bins and the range matches to those in config file
		FILE *wfile = fopen(w->weightfile, "r");

		// read character by character and look for newline
		char c = getc(wfile);
		long lines = 0;
		while (c != EOF) {
        if (c == '\n') {
          lines++;
        }
        // next
        c = getc(wfile);
    }

		if (lines != w->bins) {
			printf0(p, "Error: Bin count weightfile does not match! Found %ld, was supposed to be %ld\n", lines, w->bins);
			die(21);
		}

    rewind(wfile);
    int read;

    for(i=0; i<w->bins; i++) {

			// read = how many elements read per in a line
      read = fscanf(wfile, "%lf %lf", &(w->pos[i]), &(w->W[i]));

      // if read != 2, there was an error reading the line, so we exit
      if(read != 2) {
				printf0(p, "Error reading weightfile! At line %ld\n", i+1);
				die(22);
      }
    }

    fclose(wfile);

		// check that the range in weight matches to what was given in config
		double readmax = w->pos[w->bins-1] + w->dbin;
		if ( fabs(w->pos[0] - w->min) > 0.00001 || fabs(readmax - w->max) > 0.00001) {
			printf0(p, "Error: Wrong range weightfile!\n Found min=%lf, max=%lf was supposed to be min=%lf, max=%lf\n", w->pos[0], readmax, w->min, w->max);
			die(23);
		}

		printf("WARNING: WEIGHT LOADING IS NOT IMPLEMENTED!!!!!!!!!\n");
		// also: FREE all multi arrays somewhere...
  }

}

/* Save the current weight into file.
* There will be w.bins lines, each containing:
* w.pos[i] w.W[i]
*
*	w.hits and w.h are NOT saved. They will reset if
* the simulation is restarted later.
*
* Original implementation by David Weir.
*/
void save_weight(params p, weight w) {

	// for readonly run, do nothing
	if (!w.readonly) {

		//printf0(p, "Saving weightfile: %s\n", w.weightfile);
		if(!p.rank) {
			FILE *wfile = fopen(w.weightfile, "w");

			for(long i=0; i<w.bins; i++) {
				fprintf(wfile, "%lf %lf\n", w.pos[i], w.W[i]);
			}

			fclose(wfile);

		}
	}
}


/* Get linearized weight corresponding to order
* parameter value val. Our convention: bin with i = 0
* contains weight for value = w.min and the last bin
* corresponds to the interval ending in w.max.
* We use constant weight in the last bin.
*/
double get_weight(weight w, double val) {

	// If out of bounds, return w.outsideW_min or w.outsideW_max
	if (val > w.max) {
		return w.outsideW_max;
	} else if (val < w.min) {
		return w.outsideW_min;
	}

	// which bin is val in?
	long bin = whichbin(w, val);

	// linearize the weight function in this bin
	double val_prev = w.min + bin * w.dbin;
	double nextW, val_next;
	if (bin >= w.bins - 1) {
		// last bin, use dummy weight at the end of range
		nextW = w.W[w.bins - 1];
		val_next = w.max;
	} else {
		nextW = w.W[bin+1];
		val_next = val_prev + w.dbin;
	}

	//printf("val = %lf, valprev = %lf, valnext = %lf, w.W[bin] = %lf, nextW = %lf, returns = %lf \n,",
	//val, val_prev, val_next, w.W[bin], nextW, w.W[bin] + (val - val_prev) * (nextW - w.W[bin]) / (val_next - val_prev));

	return w.W[bin] + (val - val_prev) * (nextW - w.W[bin]) / (val_next - val_prev);
}


/* Measure current value of the order parameter and its canonical histogram.
* This routine is meant to be ran together with the usual measure() routine.
* Could ran more often, but no point if autocorrelations are large.
*/
void measure_muca(params p, fields f, weight* w) {

	if (w->readonly)
		return;

	double val = calc_orderparam(p, f);

	// only measure if we are in binning range
	if (!(val > w->max || val < w->min)) {

		long bin = whichbin(*w, val);
		w->hits[bin]++;
		double current_W = get_weight(*w, val);
		w->h[bin] += exp(-current_W);
	}

	w->m++;
	if (w->m >= w->update_interval) {
		update_weight(p, w);
	}

}

/* Update weight at all bins using the recursive method.
*/
void update_weight(params p, weight* w) {

	if (w->readonly) {
		return;
	}

	/* recursion for bin i is:
	*		w^(k+1)_i - w^(k-1)_i = (sum_{n<=k} g^n_{i} ln(h^n_{i-1} / h^n_i)) / sum_{n<=k} g^n_{i},
	* where sum_{n<=k} g^n_{i} is stored in w.gsum[i].
	* We need to be careful with the boundary bins!
 	*/

	int nmin = 2;
	long g;

	// old w_1 - w_0
	double dw = w->W[1] - w->W[0];
	// boundary test: some value...
	w->W[0] += w->hits[0];

	// now calculate weight in other bins relative to the first one
	for (long i=1; i<w->bins; i++) {
		// old value for the sum in the numerator
		double num = dw * w->gsum[i];
		// next iteration needs OLD dw = w_{i+1} - w_i (unless last bin)
		if (i < w->bins - 1) {
			dw = w->W[i+1] - w->W[i];
		}

		// new g factor? g_i = n_{i-1} + n_i if n_i, n_{i-1} > n_min, g=0 otherwise
		if (w->hits[i-1] <= nmin || w->hits[i] <= nmin) {
			g = 0;
		} else {
			g = w->hits[i-1] + w->hits[i];
		}
		w->gsum[i] += g;

		// if gsum is 0, we do not have enough measurement for this bin
		if (w->gsum[i] == 0 || w->h[i] == 0 || w->h[i-1] == 0) {
			// do not update??
			//w->W[i] = w->W[i-1];
		} else {
			// new numerator
			num += g * log(w->h[i-1] / w->h[i]);
			w->W[i] = w->W[i-1] + num / w->gsum[i];
		}

	}

	// reset muca loop and save weight
	w->m = 0;
	for (long i=0; i<w->bins; i++) {
		w->hits[i] = 0;
		w->h[i] = 0.0;
	}
	save_weight(p, *w);

/* TODO tunneling check
	// did we tunnel? (the numbers here are arbitrary...)
	// First and last how many bins count as being in the other minimum?
	long tunnel_threshold = (long) (w->bins / 100.0 * 2.0);
	int tunnel = 0;
	if (w->last_max && bin < tunnel_threshold ) {
		// tunneled from order param = max to order param = min
		w->last_max = 0;
		tunnel++;
	} else if (!w->last_max && bin > w->bins - tunnel_threshold) {
		// tunneled from min to max
		w->last_max = 1;
		tunnel++;
	}

	if (tunnel > 0) {
		w->increment /= 1.5;
		printf0(p, "Reducing weight update factor! Now %lf \n", w->increment);
	}

*/

}



/* Global accept/reject step for multicanonical updating.
* oldval is the old order parameter value before field was updated locally
* Return 1 if update was accepted, 0 otherwise.
* TODO specify arbitrary order param. need some function that calculates the order param value when I give it some label as an argument...
* TODO optimize calculation of order param diff
*/
int multicanonical_acceptance(params p, weight* w, double oldval, double newval) {
	// acc/rej only in root node
	int accept;
	if (p.rank == 0) {

		double W_new, W_old;

		W_new = get_weight(*w, newval);
		W_old = get_weight(*w, oldval);

		// exp(-S + W) ->
		double diff = W_new - W_old;
		if (diff >= 0) {
			accept = 1;
		} else if(exp((diff)) > drand48()) {
      accept = 1;
    } else {
      accept = 0;
    }
	}

	// broadcast outcome to all nodes
	bcast_int(&accept);

	/*
	// if accepted, increase weight for the new order param value to make it less attractive
	if (!w->readonly && accept) {
		update_weight(p, w, newval);
	}
	*/

	return accept;
}


/* Return bin index corresponding to given value of order parameter.
* If out of range, we return the closest bin index. The calling
* functions should ensure that this does not happen, however,
* as it may lead to accidents....
*/
long whichbin(weight w, double val) {
	if (val < w.min) {
		return 0;
	} else if (val >= w.max) {
		return w.bins - 1;
	} else {
		return (long) ((val - w.min) / w.dbin);
	}
}


// Calculate muca order parameter and distribute to all nodes
double calc_orderparam(params p, fields f) {
	double tot = 0.0;
	for (long i=0; i<p.sites; i++) {
		tot += doubletsq(f.su2doublet[i]);
	}
	tot /= p.vol;
	return allreduce(tot);
}

// Backup a field in case a sweep needs to be undone later

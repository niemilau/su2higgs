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
* Following the sign convention of Kari and David,
* our weight function INCREASES when we want to suppress a particular configuration (ensemble ~ exp(-S - W).
* When a bin is visited by the order parameter, the weight in that bin is increased by w.increment.
*
*	At least two good options to calculate the weight:
* 	1. Recursive computation described in hep-lat/9804019. Hard to implement, but
*			 will eventually converge to the exact weight function. Here weight in the first bin can
*			 be chosen arbitrary (e.g. 0), as weight elsewhere is computed relative to the previous bin.
*		2. Simpler method similar to the Landau-Wang algorithm, but information from
*			 previous short runs is "forgotten" in each iteration. Kari uses this in his
*			 modern codes. Easy to implement, but may not give the exact weight if the
*			 "delta" parameter becomes too small too early. Usually still sufficient.

* We use method 2 here. To summarize the algorithm:
*
* 	1. Do a standard checkerboard sweep, updating the multicanonical field with
*			 local updates WITHOUT caring about the weight at all. This will move the system
*			 towards a local minimum.
*		2. After the sweep, recalculate the order parameter and do a global accept/reject
*			 based on the weight change. If rejected, ALL updates in the sweep are undone.
*			 This will create a bias towards non-minimum configurations. Increase w.hits[bin] by one.
*
*			 Note that if phisq is the order param
*			 and we first sweep and update phi, then e.g. gauge links and only then
*			 do multicanonical accept/reject, even the link updates need to be undone
*			 in case of reject. In practice such updating does not make much sense, though.

*		3. After certain number of multicanonical measurements (can measure after every sweep,
*			 according to Kari; autocorrelations don't matter that much here), increase weight in
*			 each bin by value proportional to w.hits[bin]. Then reset w.hits; this makes the simulation
*			 "forget" the earlier runs.
*		4. Once all bins have been visited, decrease w.increment to make the iteration converge.
*			 This is necessary because at first the system prefers the canonical minima, so those
*			 obtain large weight very quickly. Afterwards, the mixed-phase configurations are
*		   preferred and the system spends most of the time there, so smaller weight updates are
*			 then required to prevent the weight from becoming flat again.
*
*			(TODO once the system tunnels once, decrease w.increment)
*
* Kari tips: Behavior of the weight function outside the given range can in principle be chosen arbitrary.
* Usually it's enough to set it to a constant value there, so that the total weight function is still continuous.
* For bubble nucleation it may be better to set an infinite weight there instead.
* For the weight range, it's best to choose it so that the peaks are mostly included
*
* TODO different order parameters. For this we need to modify calc_orderparam()
* but not others (apart from possible minor changes). We also need to place the
* call to muca_acceptance() appropriately (this is hard to automatize!)
*
*/

#include "su2.h"


/* Save the current weight into file.
* There will be w.bins + 1 lines, the first one being:
* 	w.bins w.increment w.last_max
* the lines after that read:
* 	w.pos[i] w.W[i]
*
* Note that w.hits is NOT stored!
*/
void save_weight(params p, weight w) {

	// for readonly run, do nothing
	if (w.readonly)
		return;

	if(!p.rank) {
		FILE *wfile = fopen(w.weightfile, "w");

		fprintf(wfile, "%ld %lf %d\n", w.bins, w.increment, w.last_max);
		for(long i=0; i<w.bins; i++) {
			fprintf(wfile, "%lf %lf\n", w.pos[i], w.W[i]);
		}

		fclose(wfile);

	}
}

/* Load multicanonical weight from weightfile.
* Assumes that binning range has already been read from config.
* If file does not exist, initializes a new flat weight.
* DO NOT call this more than once per run.
*
* Original implementation by David Weir.
*/
void load_weight(params p, weight *w) {

	printf0(p, "\nLoading weightfile: %s\n", w->weightfile);

  w->pos = malloc(w->bins * sizeof(*(w->pos)));
  w->W = malloc(w->bins * sizeof(*(w->W)));
	w->hits = malloc(w->bins * sizeof(*(w->hits)));

	w->dbin = (w->max - w->min)/((double)w->bins);

	w->update_interval = 500; // Kari used 500-2000 in MSSM

  long i;

  if(access(w->weightfile, R_OK) != 0) {
		// no weight found
    printf0(p,"Unable to access weightfile!!\n");
		if (w->readonly) {
			printf0(p, "No multicanonical weight given for a read-only run! Exiting...\n");
			die(20);
		}

    for(i=0; i<w->bins; i++) {
      w->pos[i] = w->min + ((double) i) * w->dbin;
      w->W[i] = 0.0;
    }
		w->last_max = 0; // assume starting from min

		printf0(p, "Initialized new weight \n");
		
  } else {
		// found existing weight, check that bins and the range matches those in config file
		FILE *wfile = fopen(w->weightfile, "r");

		// check that binning matches to values in config
		// read = how many elements read per in a line
		int read;
		long bins_read;


		// first line:
		read = fscanf(wfile, "%ld %lf %d", &bins_read, &w->increment, &w->last_max);

		if (read != 3 || bins_read != w->bins) {
			printf0(p, "Error reading first line of weightfile! Check bins \n");
			die(22);
		}

    for(i=0; i<w->bins; i++) {
      read = fscanf(wfile, "%lf %lf", &(w->pos[i]), &(w->W[i]));
      if(read != 2) {
				printf0(p, "Error reading weightfile! Got %d values at line %ld\n", read, i+2);
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

		printf0(p, "Starting with weight increment %lf, last_max %d\n", w->increment, w->last_max);
		// also: free muca arrays somewhere

  }

	// restart accumulation of muca hits even if existing weight is loaded
	w->m = 0;
	for (i=0; i<w->bins; i++) {
		w->hits[i] = 0;
	}

	// weight outside the binning range?
	if (w->absolute_bounds) {
		w->outsideW_max = 500;
		w->outsideW_min = 500;
	} else {
		w->outsideW_max = w->W[w->bins-1];
		w->outsideW_min = w->W[0];
	}

}


/* Get linearized weight corresponding to order
* parameter value val. Our convention: bin with i = 0
* contains weight for value = w.min and the last bin
* corresponds to the interval ending in w.max.
* For simplicity, we use constant weight in the last bin.
*/
double get_weight(weight w, double val) {

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
		// last bin, use constant weight
		nextW = w.W[w.bins - 1];
		val_next = w.max;
	} else {
		nextW = w.W[bin+1];
		val_next = val_prev + w.dbin;
	}

	return w.W[bin] + (val - val_prev) * (nextW - w.W[bin]) / (val_next - val_prev);
}



/* Update weight in each bin based on hits in the bin,
* and reset the "accumulation" array w->hits. The increment factor
* is normalized by the total number of bins. This gives a much
* smoother weight function, avoiding "spikes" that can lead to low
* multicanonical acceptance. After updating, checks whether
* w->increment should be decreased for the next loop.
*/
void update_weight(params p, weight* w) {

	if (w->readonly) {
		return;
	}

	for (long i=0; i<w->bins; i++) {
		w->W[i] += w->hits[i] * w->increment / w->bins;
	}

	if (!w->absolute_bounds) {
		w->outsideW_max = w->W[w->bins - 1];
		w->outsideW_min = w->W[0];
	}

	check_tunnel(p, w);

	for (long i=0; i<w->bins; i++) {
		w->hits[i] = 0;
	}

}

/* Routine for decreasing weight increment factor.
*	The condition we assume here is that once the system has visited
* both first and last bins, we consider the system to have tunneled
* through the barrier and decrease w->increment.
* Note that this count is not reset when weight is updated.
*/
void check_tunnel(params p, weight *w) {

	// First and last how many bins count as the endpoint?
	//long tunnel_threshold = (long) (w->bins / 100.0 * 2.0);
	int tunnel = 0;
	if (w->last_max && w->hits[0] > 0) {
		// tunneled from order param = max to order param = min
		w->last_max = 0;
		tunnel++;
	} else if (!w->last_max && w->hits[w->bins-1] > 0) {
		// tunneled from min to max
		w->last_max = 1;
		tunnel++;
	}

	if (tunnel > 0) {
		w->increment /= 1.5;
		printf0(p, "\nReducing weight update factor! Now %lf \n", w->increment);
	}

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

		// TODO optimize: get bin index here?

		W_new = get_weight(*w, newval);
		W_old = get_weight(*w, oldval);

		double diff = W_new - W_old;

		if(exp(-(diff)) > drand48()) {
      accept = 1;
    } else {
      accept = 0;
    }
	}

	// broadcast outcome to all nodes
	bcast_int(&accept);

	// update w->hits and check if w->increment should be modified
	// do this even if the update was rejected
	if (!w->readonly && newval <= w->max && newval >= w->min) {
		long bin = whichbin(*w, newval);

		w->hits[bin]++;
		w->m++;
		if (w->m >= w->update_interval) {
			// update weight function, save it and start over
			update_weight(p, w);
			save_weight(p, *w);
			w->m = 0;
		}
	}

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

// Free memory allocated for multicanonical
void free_muca_arrays(weight *w) {
	free(w->pos);
  free(w->W);
	free(w->hits);
}

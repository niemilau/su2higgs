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
* We use constant weight also in the last bin!!
*
* Following the sign convention of Kari and David,
* our weight function INCREASES when we want to suppress a particular configuration (ensemble ~ exp(-S - W).
* When a bin is visited by the order parameter, the weight in that bin is increased by w.increment.
*
* It is very difficult to fully parallelize the multicanonical method, because local updates depend on the global
* weight function. Could perhaps try treating each node as a sublattice that calculates their own weight, and
* then combine them every so often? The sublattices are not independent, however...
* The implementation here essentially serializes updates for the field that acts as the multicanonical order parameter
*
* Plan: Use recursive computation of the weight function described
* in Kari's notes https://www.mv.helsinki.fi/home/rummukai/lectures/montecarlo_oulu/lectures/mc_notes6.pdf
* (see also https://www.mv.helsinki.fi/home/rummukai/simu/multican.pdf).
* Detailed procedure, inspired by Kari:
*
*		1. Write the weight as a piecewise linear function in each bin: weight[i] = k[i] * x + W[i],
*			 where x is the order parameter. Outside the binning range we choose a constant weight,
*			 equal to b[i] in the nearest bin. (For bubble nucleation it's better to choose infinite
*			 weight instead).
*		2. Do a standard checkerboard sweep (fixed parity) on the multicanonical field,
*			 updating sites locally using standard metropolis/overrelax. Here the only difference
*			 to canonical updating is that the weight DOES need to be taken into account. We do this by
*			 a redefinition of the mass term: if order parameter is phisq, then msq -> msq + k[i] in bin i. (WHICH SIGN?!)
*			 We DO NOT update the weight function itself here, because this would not parallelize well.
*		3. After the sweep, we calculate the total difference in the order parameter, and perform
*			 a global accept/reject (let root node do this) based on the weight change. If accepted,
*			 we increment the weight in the new bin to make this bin less attractive in the future.
*			 If rejected, ALL local updates need to be undone.
*		4. Once all bins have been visited, decrease the weight increment factor so that the iteration converges.
*
*	Basic routines:
*		multicanonical_acceptance() does the acc/rej and calls update_weight() if necessary
*		calculate_weight() linearizes weight in a given bin and redefines action parameters as required
*
* Note: For the weight range, it's best to choose it so that the peaks are mostly included
*
* What Kari does: sweep over half of the sites and update scalars locally, then perform global weight check
* Note that if also other fields, eg gauge links, are updated BEFORE the muca check,
* their changes need to be undone as well if the muca check is rejected
*
*/

#include "su2.h"


/* Load multicanonical weight from weightfile.
* If file does not exist, initializes a new flat weight.
* The weight can be saved with save_weight().
* DO NOT call this more than once per run.
*
* Original implementation by David Weir.
*/
void load_weight(params p, weight *w) {

	printf0(p, "\nLoading weightfile: %s\n", w->weightfile);

  w->pos = malloc(w->bins * sizeof(*(w->pos)));
  w->W = malloc(w->bins * sizeof(*(w->W)));
	w->slope = malloc(w->bins * sizeof(*(w->slope)));
	w->visited = malloc(w->bins * sizeof(*(w->visited)));

	w->dbin = (w->max - w->min)/((double)w->bins);


  long i;

  if(access(w->weightfile, R_OK) != 0) {

    printf0(p,"Unable to access weightfile!!\n");

		if (w->readonly) {
			printf0(p, "No multicanonical weight given for a read-only run! Exiting...\n");
			die(20);
		}

    for(i=0; i<w->bins; i++) {
      w->pos[i] = w->min + ((double) i) * w->dbin;
      w->W[i] = 0.0;
			w->visited[i] = 0;
			w->slope[i] = 0.0;
    }
		w->visited_total = 0;

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
      read = fscanf(wfile, "%lf %lf %d", &(w->pos[i]), &(w->W[i]), &(w->visited[i]));

      // if read != 3, there was an error reading the line, so we exit
      if(read != 3) {
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

		// weight loaded, now calculate slopes and redefine action parameters, if necessary
		for(i=0; i<w->bins; i++) {
			w->slope[i] = calc_weight_slope(*w, i);
		}
  }
}

/* Save the current weight into file.
* There will be w.bins lines, each containing:
* w.pos[i] w.W[i] w.visited[i]
*
* Original implementation by David Weir.
*/
void save_weight(params p, weight w) {

	// for readonly run, do nothing
	if (!w.readonly) {

		printf0(p, "Saving weightfile: %s\n", w.weightfile);
		if(!p.rank) {
			FILE *wfile = fopen(w.weightfile, "w");

			for(long i=0; i<w.bins; i++) {
				fprintf(wfile, "%lf %lf %d\n", w.pos[i], w.W[i], w.visited[i]);
			}

			fclose(wfile);

		}
	}
}


/* Get linearized weight corresponding to order
* parameter value val. Our convention: bin with i = 0
* contains weight for value = w.min and the last bin
* corresponds to the interval ending in w.max.
* To linearize the last bin, we thus need some dummy value
* for the weight at w.max
*/
double get_weight(weight w, double pos) {

	// outside the weighting range we have constant weight so that
	// the total weight function is continuous. Hence, return weight of the
	// nearest bin if input value is out of range.
	if (pos > w.max) {
		return w.W[w.bins - 1];
	} else if (pos < w.min) {
		return w.W[0];
	}

	// which bin is val in?
	long bin = whichbin(w, pos);

	if (bin >= w.bins - 1) {
		// last bin, constant weight
		return w.W[w.bins - 1];
	} else {
		// linearize weight in this bin
		return w.slope[bin] * (pos - w.pos[bin]) + w.W[bin];
	}

}



/* Update weight at a given bin by the increment factor.
* This is essentially a "visit bin" -routine, meaning that once
* a bin is updated, we also set w.visited[i] = 1 and increase visited_total.
* If all sites have now been visited, decrease the increment factor and reset visited array.
*/
void update_weight(params* p, weight* w, double pos) {

	// do not update if order param is out of bounds
	if (w->readonly || pos > w->max || pos < w->min) {
		return;
	}

	long bin = whichbin(*w, pos);

	w->W[bin] += w->increment;
	// calculate new slope in new and redefine action parameter(s).
	// note that changing w.W[i] also changes the slope in bin i-1
	double oldslope = w->slope[bin];
	w->slope[bin] = calc_weight_slope(*w, bin);
	if (bin != 0) {
		w->slope[bin-1] = calc_weight_slope(*w, bin-1);
	}
		printf("msq now = %lf\n", p->msq_phi);
	muca_param_shift(p, *w, bin, oldslope);

	if (!(w->visited[bin])) {
		// first time visiting this bin
		w->visited[bin] = 1;
		w->visited_total++;

		if (w->visited_total >= w->bins) {
			// visited all, decrease increment factor and reset
			w->increment /= 1.5;
			printf0(*p, "Reducing weight update factor! Now %lf \n", w->increment);
			for (long i=0; i<w->bins; i++) {
				w->visited[i] = 0;
			}
			w->visited_total = 0;
		}
	}
}



/* Global accept/reject step for multicanonical updating.
* oldval is the old order parameter value before field was updated locally
* Return 1 if update was accepted, 0 otherwise.
* TODO specify arbitrary order param. need some function that calculates the order param value when I give it some label as an argument...
* TODO optimize calculation of order param diff
*/
int multicanonical_acceptance(params* p, weight* w, double oldval, double newval) {
	// acc/rej only in root node
	int accept;
	if (p->rank == 0) {

		double W_new, W_old;

		W_new = get_weight(*w, newval);
		W_old = get_weight(*w, oldval);

		printf0(*p,"oldval = %lf, newval = %lf, W_old = %lf, W_new = %lf \n", oldval, newval, W_old, W_new);

		double diff = W_new - W_old;
		if (diff <= 0) {
			accept = 1;
		} else if(exp(-(diff)) > drand48()) {
      accept = 1;
    } else {
      accept = 0;
    }
	}

	// broadcast outcome to all nodes
	bcast_int(&accept);

	// if accepted, increase weight for the new order param value to make it less attractive
	if (!w->readonly && accept) {
		update_weight(p, w, newval);
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


/* Calculate linearized weight slope in a given bin, but DO NOT
* update it in w->slope. Updating should only be done in update_weight()
* or in initialization, which is where this routine is useful.
*/
double calc_weight_slope(weight w, long bin) {
	if (bin < 0 || bin >= w.bins) {
		// should not get here!
		return 0;
	}

	double k;
	if (bin == w.bins - 1) {
		// last bin, use constant weight
		k = 0;
	} else {
		k =	(w.W[bin+1] - w.W[bin] ) / (w.pos[bin+1] - w.pos[bin]);
	}

	return k;
}


/* Redefine a parameter in the action to account for the multicanonical weight.
* The point here is that since our distribution is exp(-S - W), W can affect
* local updates of the sites even if it stays constant in the update. To account
* for this, we linearize W as weight = w.slope[i] * x + w.W[i] in bin i (x = order param),
* and absorb the slope into a parameter redefinition, ex: mass term correction.
* This is useful because then the local update routines do not need to be modified at all,
* since we keep the weight constant there.
* The implementation of this routine depends heavily on the choice of the order parameter!
* PROBLEM: what if the local update changes the order param so much that the bin changes??
* PROBLEM: p->msq may accumulate errors??
*/
void muca_param_shift(params* p, weight w, long bin, double oldslope) {

	// if orderparam == phisq ( WHICH SIGN?!?!)
	p->msq_phi = p->msq_phi + 0*oldslope - 0*w.slope[bin];
}

// Backup a field in case a sweep needs to be undone later

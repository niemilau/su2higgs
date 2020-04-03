/** @file multicanonical.c
*
* Routines for multicanonical simulations.
*
* We assume that a range of order parameter values (min, max) is given
* in the config file, and restrict weighting to this range.
* This range is divided into given number of bins.
* Bin width does NOT need to be the same in all bins!
*
* Weight is stored in struct weight, w.pos[i] being the "starting" value of
* i. bin and w.W[i] being the value of the weight function at w.pos[i].
* Inside the bin, we perform a linear approximation to obtain a continuous weight function.
*
* Following the sign convention of Kari and David,
* our weight function INCREASES when we want to suppress a particular configuration (ensemble ~ exp(-S - W).
* When a bin is visited by the order parameter, the weight in that bin is increased by w.increment.
*
* Summary of the algorithm:
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
*
*	Additional steps for calculating the weight function, if w.readonly != 1:

*		3. After certain number of multicanonical measurements (can measure after every sweep,
*			 according to Kari; autocorrelations don't matter that much here), increase weight in
*			 each bin by value proportional to w.hits[bin]. Then reset w.hits; this makes the simulation
*			 "forget" the earlier runs.
*		4. Once the system has visited both the first and last bin, decrease w.increment to make the
*			 iteration converge. This is necessary because at first the system prefers the canonical
*			 minima, so those obtain large weight very quickly. Afterwards, the mixed-phase configurations are
*		   preferred and the system spends most of the time there, so smaller weight updates are
*			 then required to prevent the weight from becoming flat again.
*
*
* Practical note: For the weight range, it's best to choose it so that the peaks are mostly included
*
***************** On different order parameters ********************
* Most of the routines here for updating the weight and accepting a multicanonical update are
* insensitive to what the actual order parameter is. To implement a new multicanonical order parameter,
* the following steps need to be taken:
*
*		1. #define a constant label for it in a header (su2.h) and add it to get_weight_parameters().
*		2. Tell calc_orderparam(), alloc_backup_arrays(), free_muca_arrays(),
*			 store_muca_fields() and reset_muca_fields() what to do for this case.
*		3. In update_lattice_muca(), call multicanonical_acceptance() after updating all the
*			 fields contributing to the order param and reject the updates if necessary (this is the tricky part).
*/

/* --- Weight logic ---
*	1. w.min and w.max determine the order parameter range where the weight is to be
		 updated (if readonly=0). the weight update factor is reduced each time the system tunnels from w.min to w.max
* 2. w.min_abs and w.max_abs determine the absolute range of weighting.
*		 Last bin should end in w.max_abs.
*		 Outside this range the weight is set to the same value as in the first or last bin.
*	2.5. these ranges are setup in load_weight(). If no weight is provided, the program will generate a flat weight
*  		 with equal bin size and assumes w.min = w.min_abs, w.max = w.max_abs, which are read from the config file.
*	3. inside a bin, a linear interpolation is used for the weight. Last bin uses constant weight.
*/

#include "su2.h"


/* Save the current weight into file.
* There will be w.bins + 1 lines, the first one being:
* 	w.bins w.increment w.last_max w.min w.min_abs w.max w.max_abs
* the lines after that read:
* 	w.pos[i] w.W[i]
*
* Note that w.hits is NOT stored!
*/
void save_weight(params const* p, weight const* w) {

	// for readonly run, do nothing
	if (w->readonly)
		return;

	if(!p->rank) {
		FILE *wfile = fopen(w->weightfile, "w");

		fprintf(wfile, "%ld %lf %d %lf %lf %lf %lf\n", w->bins, w->increment, w->last_max, w->min, w->max, w->min_abs, w->max_abs);
		for(long i=0; i<w->bins; i++) {
			fprintf(wfile, "%lf %lf\n", w->pos[i], w->W[i]);
		}

		fclose(wfile);

	}
}

/* Load multicanonical weight from weightfile.
* If file does not exist, initializes a new flat weight.
* Assumes that bins, min and max have been read from config file,
* but overrides these if an existing weight file is found.
* DO NOT call this more than once per run.
*/
void load_weight(params const* p, weight *w) {

	printf0(*p, "\nLoading weightfile: %s\n", w->weightfile);

  w->pos = malloc(w->bins * sizeof(*(w->pos)));
  w->W = malloc(w->bins * sizeof(*(w->W)));
	w->hits = malloc(w->bins * sizeof(*(w->hits)));

	w->update_interval = 1000; // Kari used 500-2000 in MSSM

	w->do_acceptance = 1;

  long i;

  if(access(w->weightfile, R_OK) != 0) {
		// no weight found; initialize a flat weight according to config file
    printf0(*p,"Unable to access weightfile!!\n");
		if (w->readonly) {
			printf0(*p, "No multicanonical weight given for a read-only run! Exiting...\n");
			die(20);
		}

		// assume that the weight update range is the same as weighting range
		w->max_abs = w->max;
		w->min_abs = w->min;

		// equally sized bins
		double dbin = (w->max - w->min)/((double) w->bins);

    for(i=0; i<w->bins; i++) {
      w->pos[i] = w->min + ((double) i) * dbin;
      w->W[i] = 0.0;
    }
		w->last_max = 0; // assume starting from min

		printf0(*p, "Initialized new weight \n");
		save_weight(p, w);

  } else {
		// found existing weight, use it instead of the one specified in config
		FILE *wfile = fopen(w->weightfile, "r");

		// read = how many elements read per in a line
		int read;
		long bins_read;

		// first line: increment, last_max and range for weight updating
		read = fscanf(wfile, "%ld %lf %d %lf %lf %lf %lf ", &bins_read, &w->increment, &w->last_max, &w->min, &w->max, &w->min_abs, &w->max_abs);

		if (read != 7) {
			printf0(*p, "Error reading first line of weightfile! \n");
			die(22);
		}

		// override initial options with the values read from the actual weightfile
		w->bins = bins_read; // includes the extra bin
		w->pos = realloc(w->pos, w->bins * sizeof(*(w->pos)));
		w->W = realloc(w->W, w->bins * sizeof(*(w->W)));
		w->hits = realloc(w->hits, w->bins * sizeof(*(w->hits)));


    for(i=0; i<w->bins; i++) {
      read = fscanf(wfile, "%lf %lf", &(w->pos[i]), &(w->W[i]));
      if(read != 2) {
				printf0(*p, "Error reading weightfile! Got %d values at line %ld\n", read, i+2);
				die(22);
      }
    }

		fclose(wfile);

		if (!w->readonly) {
			if (w->min < w->min_abs || w->max > w->max_abs) {
				printf0(*p, "Weight error! Weight update range is larger than weighting range!\n");
				printf0(*p, "Got w.min %lf, w.min_abs %lf ; w.max %lf, w.max_abs %lf \n", w->min, w->min_abs, w->max, w->max_abs);
				die(24);
			}
		}

  } // weight loading/initialization OK

	printf0(*p, "Using weight function with %ld bins in range %lf, %lf\n", w->bins, w->min_abs, w->max_abs);
	if (!w->readonly) {
		printf0(*p, "Will modify weight in range %lf, %lf\n", w->min, w->max);
		printf0(*p, "starting with increment %lf, last_max %d\n", w->increment, w->last_max);
	}

	// restart accumulation of muca hits even if existing weight is loaded
	w->m = 0;
	for (i=0; i<w->bins; i++) {
		w->hits[i] = 0;
	}

	// find bin indices for the bins containing min and max values of the update range
	w->min_bin = whichbin(w, w->min);
	w->max_bin = whichbin(w, w->max);

}


/* Get linearized weight corresponding to order
* parameter value val.
*/
double get_weight(weight const* w, double val) {

	// if outside binning range, use same weight as in first or last bin
	if (val > w->max_abs) {
		return w->W[w->bins-1];
	} else if (val < w->min_abs) {
		return w->W[0];
	}

	// which bin is val in?
	long bin = whichbin(w, val);

	// linearize the weight function in this bin
	double val_prev = w->pos[bin];
	double nextW, val_next;
	if (bin >= w->bins - 1) {
		// last bin, use constant weight
		nextW = w->W[w->bins - 1];
		val_next = w->max_abs;
	} else {
		nextW = w->W[bin+1];
		val_next = w->pos[bin+1];
	}

	return w->W[bin] + (val - val_prev) * (nextW - w->W[bin]) / (val_next - val_prev);
}



/* Update weight in each bin based on hits in the bin,
* and reset the "accumulation" array w->hits. The increment factor
* is normalized by the total number of bins. This gives a much
* smoother weight function, avoiding "spikes" that can lead to low
* multicanonical acceptance. After updating, checks whether
* w->increment should be decreased for the next loop.
*/
void update_weight(params const* p, weight* w) {

	if (w->readonly) {
		return;
	}

	for (long i=0; i<w->bins; i++) {
		w->W[i] += w->hits[i] * w->increment / w->bins;
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
void check_tunnel(params const* p, weight* w) {

	int tunnel = 0;
	if (w->last_max && w->hits[w->min_bin] > 0) {
		// tunneled from order param = max to order param = min
		w->last_max = 0;
		tunnel++;
	} else if (!w->last_max && w->hits[w->max_bin] > 0) {
		// tunneled from min to max
		w->last_max = 1;
		tunnel++;
	}

	if (tunnel > 0) {
		w->increment *= w->reduction_factor;
		printf0(*p, "\nReducing weight update factor! Now %lf \n", w->increment);
	}

}


/* Global accept/reject step for multicanonical updating.
* oldval is the old order parameter value before field was updated locally
* Return 1 if update was accepted, 0 otherwise.
*/
int multicanonical_acceptance(params const* p, weight* w, double oldval, double newval) {

	// if we call this function while w->do_acceptance is 0 then something went wrong
	if (!w->do_acceptance) {
		printf0(*p, "Should not get here!! in multicanonical.c\n");
		die(-1000);
	}

	// if the new value is outside the weight range and the restrict flags are on,
	// automatically reject the update. No comms needed for this
	if (w->restrict_min && oldval < w->min_abs) {
		return 0;
	} else if (w->restrict_max && oldval > w->max_abs) {
		return 0;
	}

	// acc/rej only in root node
	int accept;
	if (p->rank == 0) {

		double W_new, W_old;

		// TODO optimize: get bin index here?

		W_new = get_weight(w, newval);
		W_old = get_weight(w, oldval);

		double diff = W_new - W_old;

		if(exp(-(diff)) > drand48()) {
      accept = 1;
    } else {
      accept = 0;
    }
	}

	// broadcast outcome to all nodes
	bcast_int(&accept);

	// update hits and call update_weight() if necessary.
	// do this even if the update was rejected
	if (!accept)
		newval = oldval;

	// if outside the actual weight update range, move on.
	// this also ensures that hits is never increased in bins outside the update range
	// (so hits in the extra bin should always stay at 0)
	if (!w->readonly && newval <= w->max && newval >= w->min) {
		long bin = whichbin(w, newval);

		w->hits[bin]++;
		w->m++;
		if (w->m >= w->update_interval) {
			// update weight function, save it and start over
			update_weight(p, w);
			save_weight(p, w);
			w->m = 0;
		}
	}

	return accept;
}


/* Return bin index corresponding to given value of order parameter.
* If out of binning range, we return the closest bin index. The calling
* functions should ensure that this does not happen, however.
*/
long whichbin(weight const* w, double val) {
	if (val < w->min_abs) {
		return 0;
	} else if (val >= w->max_abs) {
		// return the last bin
		return w->bins-1;
	} else {
		// do a quick binary search
		long bin;
		long bmin = 0; long bmax = w->bins-1;
		double current, next;

		while (bmin <= bmax) {
			bin = bmin + (bmax - bmin) / 2;
			current = w->pos[bin];

			if (bin == w->bins-1) {
				// last proper bin
				next = w->max_abs;
			} else {
				next = w->pos[bin+1];
			}

			if ((val < next && val > current) || (fabs(val - current) < 1e-10) ) {
				return bin;
			}

			if (val > current) {
				bmin = bin + 1;
			} else {
				bmax = bin - 1;
			}
		}

		// end binary search, if we got here then something went wrong!
		printf("Error in multicanonical.c: failed to find bin!\n");
		return w->bins - 1;
	}
}


/* Calculate muca order parameter and distribute to all nodes.
* Only the contribution from sites with parity = par is recalculated
* while the other parity contribution is read from w.param_value
*/
double calc_orderparam(params const* p, fields const* f, weight* w, char par) {
	double tot = 0.0;
	long offset, max;
	if (par == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	switch(w->orderparam) {
		case SIGMASQ :
			for (long i=offset; i<max; i++) {
				tot += tripletsq(f->su2triplet[i]);
			}
			break;
		case PHI2SIGMA2 :
			for (long i=offset; i<max; i++) {
				tot += doubletsq(f->su2doublet[i]) * tripletsq(f->su2triplet[i]);
			}
			break;
    case PHI2MINUSSIGMA2 :
      for (long i=offset; i<max; i++) {
				tot += doubletsq(f->su2doublet[i]) - tripletsq(f->su2triplet[i]);
      }
      break;
		case PHISQ :
			for (long i=offset; i<max; i++) {
				tot += doubletsq(f->su2doublet[i]);
			}
			break;
	}

	tot = allreduce(tot) / p->vol;

	w->param_value[par] = tot;
	// add other parity contribution
	return tot + w->param_value[ otherparity(par) ];
}


/* Backup all fields contributing to the order parameter,
* in case an update sweep is rejected later by the global acc/rej.
* Ordering in the backup array is exactly the same as in the original field array.
* TODO: only store sites of a given parity?
*/
void store_muca_fields(params const* p, fields* f, weight* w) {

	switch(w->orderparam) {
		case SIGMASQ :
			for (long i=0; i<p->sites; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->backup_triplet[i][dof] = f->su2triplet[i][dof];
				}
			}
			break;
    case PHI2MINUSSIGMA2 :
		case PHI2SIGMA2 :
			for (long i=0; i<p->sites; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->backup_triplet[i][dof] = f->su2triplet[i][dof];
				}
			}
			// continue to store Higgs
		case PHISQ :
			for (long i=0; i<p->sites; i++) {
				for (int dof=0; dof<SU2DB; dof++) {
					f->backup_doublet[i][dof] = f->su2doublet[i][dof];
				}
			}
			break;
	}
}


/* Undo field updates if the multicanonical step is rejected.
*/
void reset_muca_fields(params const* p, fields* f, weight* w, char par) {

	long offset, max;
	if (par == EVEN) {
		offset = 0; max = p->evensites;
	} else {
		offset = p->evensites; max = p->sites;
	}

	switch(w->orderparam) {
		case SIGMASQ :
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->su2triplet[i][dof] = f->backup_triplet[i][dof];
				}
			}
			break;
    case PHI2MINUSSIGMA2 :
		case PHI2SIGMA2 :
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->su2triplet[i][dof] = f->backup_triplet[i][dof];
				}
			}
		// continue to undo Higgs changes
		case PHISQ :
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2DB; dof++) {
					f->su2doublet[i][dof] = f->backup_doublet[i][dof];
				}
			}
			break;
	}
}


// Allocate field backup arrays
void alloc_backup_arrays(params const* p, fields* f, weight const* w) {
	switch(w->orderparam) {
		case SIGMASQ :
			f->backup_triplet = make_field(p->sites, SU2TRIP);
			break;
		case PHI2SIGMA2 :
    case PHI2MINUSSIGMA2 :
			f->backup_triplet = make_field(p->sites, SU2TRIP);
		case PHISQ :
			f->backup_doublet = make_field(p->sites, SU2DB);
			break;
	}
}

// Free memory allocated for multicanonical
void free_muca_arrays(fields* f, weight *w) {
	free(w->pos);
  free(w->W);
	free(w->hits);

	// free backup arrays
	switch(w->orderparam) {
		case SIGMASQ :
			free_field(f->backup_triplet);
			break;
		case PHI2SIGMA2 :
    case PHI2MINUSSIGMA2 :
			free_field(f->backup_triplet);
		case PHISQ :
			free_field(f->backup_doublet);
			break;
	}
}

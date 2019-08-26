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

	// update w->hits and call update_weight() if necessary.
	// do this even if the update was rejected
	if (!accept) 
		newval = oldval;
	
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


/* Calculate muca order parameter and distribute to all nodes.
* Only the contribution from sites with parity = par is recalculated
* while the other parity contribution is read from w.param_value
*/
double calc_orderparam(params p, fields f, weight* w, char par) {
	double tot = 0.0;
	long offset, max;
	if (par == EVEN) {
		offset = 0; max = p.evensites;
	} else {
		offset = p.evensites; max = p.sites;
	}

	switch(w->orderparam) {
		case SIGMASQ :
			for (long i=offset; i<max; i++) {
				tot += tripletsq(f.su2triplet[i]);
			}
			break;
		case PHI2SIGMA2 :
			for (long i=offset; i<max; i++) {
				tot += doubletsq(f.su2doublet[i]) * tripletsq(f.su2triplet[i]);
			}
			break;
		case PHISQ :
			for (long i=offset; i<max; i++) {
				tot += doubletsq(f.su2doublet[i]);
			}
			break;
	}
	
	tot = allreduce(tot) / p.vol;
	
	w->param_value[par] = tot;
	// add other parity contribution 
	return tot + w->param_value[ otherparity(par) ];
}


/* Backup all fields contributing to the order parameter, 
* in case an update sweep is rejected later by the global acc/rej.
* Ordering in the backup array is exactly the same as in the original field array.
* Depending on the order parameter, we may have to also backup halo fields.
* Some optimization is achieved by only including halos here if it is scrictly necessary.
*/
void store_muca_fields(params p, fields* f, weight* w) {
	
	switch(w->orderparam) {
		case SIGMASQ : 
			// no halos
			for (long i=0; i<p.sites; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->backup_triplet[i][dof] = f->su2triplet[i][dof];
				}
			}
			break;
		case PHI2SIGMA2 :
			// need real sites AND halos
			for (long i=0; i<p.sites_total; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->backup_triplet[i][dof] = f->su2triplet[i][dof];
				}
			}
			// continue to store Higgs
		case PHISQ : 
			// no halos 
			for (long i=0; i<p.sites; i++) {
				for (int dof=0; dof<SU2DB; dof++) {
					f->backup_doublet[i][dof] = f->su2doublet[i][dof];
				}
			}
			break;
	}
}


/* Undo field updates if the multicanonical step is rejected.
* Here we may need to undo halo updates as well, depending on 
* the order parameter (see store_muca_fields()). For optimization, 
* we do NOT include halos for the Higgs, so it should be updated last in a sweep. 
*/
void reset_muca_fields(params p, fields* f, weight* w, char par) {

	long offset, max, halo_offset, halo_max;
	if (par == EVEN) {
		offset = 0; max = p.evensites; halo_offset = p.sites; halo_max = p.sites + p.evenhalos;
	} else {
		offset = p.evensites; max = p.sites; halo_offset = p.sites + p.evenhalos; halo_max = p.sites_total;
	}
	
	switch(w->orderparam) {
		case SIGMASQ : 
			// no need to undo halos
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->su2triplet[i][dof] = f->backup_triplet[i][dof];
				}
			}
			break;
		case PHI2SIGMA2 :
			// need real sites AND halos
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->su2triplet[i][dof] = f->backup_triplet[i][dof];
				}
			}
			// then halos 
			for (long i=halo_offset; i<halo_max; i++) {
				for (int dof=0; dof<SU2TRIP; dof++) {
					f->su2triplet[i][dof] = f->backup_triplet[i][dof];
				}
			}
		// continue to undo Higgs changes
		case PHISQ : 
			// no need to undo halos
			for (long i=offset; i<max; i++) {
				for (int dof=0; dof<SU2DB; dof++) {
					f->su2doublet[i][dof] = f->backup_doublet[i][dof];
				}
			}
			break;
	} 
}


// Allocate field backup arrays
void alloc_backup_arrays(params p, fields* f, weight w) {
	switch(w.orderparam) {
		case SIGMASQ : 
			f->backup_triplet = make_field(p.sites_total, SU2TRIP);
			break;
		case PHI2SIGMA2 : 
			f->backup_triplet = make_field(p.sites_total, SU2TRIP);
		case PHISQ : 
			f->backup_doublet = make_field(p.sites, SU2DB);
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
			free_field(f->backup_triplet);
		case PHISQ : 
			free_field(f->backup_doublet);
			break;
	}
}

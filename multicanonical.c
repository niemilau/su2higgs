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
* our weight function INCREASES when we want to suppress a particular configuration. 
* When a bin is visited by the order parameter, the weight in that bin is increased by w.increment.
* (to make simulation converge, increment needs to become SMALLER when certain conditions are met, how to implement?)
*
* Plan for updating sites: 
* 	1. First do a standard metropolis or overrelax update by taking action + weight to be the full action. 
*		2. 
* 
//TODO 
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

/* Update weight at a given bin by the increment factor.
* This is essentially a "visit bin" -routine, meaning that once 
* a bin is updated, we also set w.visited[i] = 1 and increase visited_total.
* If all sites have now been visited, decrease the increment factor and reset visited array.
*	(TODO)
*/
void update_weight(weight* w, long bin) {
	
	if (w->readonly) {
		return;
	}
	
	w->W[bin] += w->increment;
	
	if (!(w->visited[bin])) {
		// first time visiting this bin
		w->visited[bin] = 1;
		w->visited_total++;
		
		if (w->visited_total >= w->bins) {
			// visited all, decrease increment factor and reset 
			w->increment /= sqrt(2);
			for (long i=0; i<w->bins; i++) {
				w->visited[i] = 0;
			}
		}
	} 
}
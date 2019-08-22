/** @file parameters.c Parses parameters from a file.
*		Original implementation by David Weir.
 */

 // TODO when parallelizing, modify the resultsfile part so that only one node does writing? See David

#include "su2.h"


/** Helper routine, checks a parameter is set and if not prints error.
 */
inline void check_set(int set, char *name) {
  if(!set) {
    fprintf(stderr,
	    "Not set parameter \"%s\"! Exiting...\n", name);
    die(1);
  }

}


/** Load parameters into struct from file.
 *
 * Parses parameters from file given by string `filename`.
 *
 * Each parameter is a key-value pair separated by whitespace. A line
 * beginning with # is treated as a comment and ignored.
 *
 * This routine also checks that no essential parameters are missing.
 */
void get_parameters(char *filename, params *p) {

	int set_dim = 0;
  int set_L = 0;

  int set_iterations = 0;
  int set_checkpoint = 0;
	int set_interval = 0;
  int set_checks = 0;

  int set_su2alg = 0;
  int set_su2DBalg = 0;
	int set_su2triplet_alg = 0;

	int set_multicanonical = 0;

  int set_seed = 0;

  int set_resultsfile = 0;

	int set_betasu2 = 0;
	// doublet
	int set_lambda_phi = 0;
	int set_msq_phi = 0;
	int set_phi0 = 0;
	// triplet
	int set_msq_triplet = 0;
	int set_a2 = 0;
	int set_b4 = 0;
	int set_sigma0 = 0;


  int set_update_doublet = 0;
	int set_update_triplet = 0;

  char key[100];
  char value[100];
  char total[200];

  int ret;

  int i;



  FILE *config;

  if(strlen(filename) == 1 && filename[0] == '-') {
    fprintf(stderr,"reading config from stdin\n");
    config = stdin;
  } else {
    if(access(filename,R_OK) != 0) {
      fprintf(stderr,"unable to access config file %s\n", filename);
      exit(101);
    }

    config = fopen(filename, "r");
  }

  while(!feof(config)) {

    if(fgets(total,198,config) == NULL) {
      // end of file?
      break;
    }

    ret = sscanf(total,"%99s%99s",key,value);

    if(ret == EOF) {
      continue;
    }

    if(key[0] == '#') {
      continue;
    }

    if(ret != 2)
      fprintf(stderr,"Unable to read input line: %s",total);

		else if(!strcasecmp(key,"dim")) {
      p->dim = strtol(value,NULL,10);
      set_dim = 1;
    }
		else if(!strcasecmp(key,"multicanonical")) {
      p->multicanonical = strtol(value,NULL,10);
      set_multicanonical = 1;
    }
		else if(!strcasecmp(key,"iterations")) {
      p->iterations = strtol(value,NULL,10);
      set_iterations = 1;
    }
		else if(!strcasecmp(key,"interval")) {
      p->interval = strtol(value,NULL,10);
      set_interval = 1;
    }
    else if(!strcasecmp(key,"checkpoint")) {
      p->checkpoint = strtol(value,NULL,10);
      set_checkpoint = 1;
    }
    else if(!strcasecmp(key,"run_checks")) {
      p->run_checks = strtol(value,NULL,10);
      set_checks = 1;
    }
    else if(!strcasecmp(key,"update_doublet")) {
      p->update_su2doublet = strtol(value,NULL,10);
      set_update_doublet = 1;
    }
		else if(!strcasecmp(key,"update_triplet")) {
      p->update_su2triplet = strtol(value,NULL,10);
      set_update_triplet = 1;
    }
		else if(!strcasecmp(key,"betasu2")) {
      p->betasu2 = strtod(value,NULL);
      set_betasu2 = 1;
    }
		// Higgs parameters
		else if(!strcasecmp(key,"msq")) {
      p->msq_phi = strtod(value,NULL);
      set_msq_phi = 1;
    }
		else if(!strcasecmp(key,"lambda")) {
      p->lambda_phi = strtod(value,NULL);
      set_lambda_phi = 1;
    }
		else if(!strcasecmp(key,"phi0")) {
      p->phi0 = strtod(value,NULL);
      set_phi0 = 1;
    }
		// triplet parameters
		else if(!strcasecmp(key,"msq_triplet")) {
      p->msq_triplet = strtod(value,NULL);
      set_msq_triplet = 1;
    }
		else if(!strcasecmp(key,"b4")) {
      p->b4 = strtod(value,NULL);
      set_b4 = 1;
    }
		else if(!strcasecmp(key,"a2")) {
      p->a2 = strtod(value,NULL);
      set_a2 = 1;
    }
		else if(!strcasecmp(key,"sigma0")) {
      p->sigma0 = strtod(value,NULL);
      set_sigma0 = 1;
    }


		else if(!strcasecmp(key,"resultsfile")) {
      if(!strcasecmp(value,"stdout")) {
				p->resultsfile = stdout;
      } else if(!strcasecmp(value,"stderr")) {
				p->resultsfile = stderr;
      } else {
				// open results file for the master node only
				if(!p->rank) {
					p->resultsfile = fopen(value, "a");
				} else {
					p->resultsfile = NULL;
				}
      }
      set_resultsfile = 1;
    }
    // Read update algorithms
    else if(!strcasecmp(key,"algorithm_su2link")) {
      if(!strcasecmp(value,"metropolis")) {
	       p->algorithm_su2link = METROPOLIS;
      } else if(!strcasecmp(value,"heatbath")) {
	       p->algorithm_su2link = HEATBATH;
      } else {
        printf("Unknown algorithm for SU(2) links, using metropolis\n");
        p->algorithm_su2link = METROPOLIS;
      }
      set_su2alg = 1;
    }
		// doublet
    else if(!strcasecmp(key,"algorithm_su2doublet")) {
      if(!strcasecmp(value,"metropolis")) {
	       p->algorithm_su2doublet = METROPOLIS;
      } else if(!strcasecmp(value,"overrelax")) {
	       p->algorithm_su2doublet = OVERRELAX;
      } else {
        printf("Unknown algorithm for SU(2) doublets, using metropolis\n");
        p->algorithm_su2doublet = METROPOLIS;
      }
      set_su2DBalg = 1;
    }
		// triplet
		else if(!strcasecmp(key,"algorithm_su2triplet")) {
      if(!strcasecmp(value,"metropolis")) {
	       p->algorithm_su2triplet = METROPOLIS;
      } else if(!strcasecmp(value,"overrelax")) {
	       p->algorithm_su2triplet = OVERRELAX;
      } else {
        printf("Unknown algorithm for SU(2) triplets, using metropolis\n");
        p->algorithm_su2triplet = METROPOLIS;
      }
      set_su2triplet_alg = 1;
    }
    // end reading update algorithms
  }

	check_set(set_dim, "dim");

	// Allocate memory for lattice side lengths.
	// This couldn't be done earlier because we didn't know the dimension
	// p.parity is allocated in alloc.c, alloc_neighbors().
	p->L = malloc(p->dim * sizeof(p->L));

	// read the file again for lattice sizes L_i. This could use some polishing
	char istr[255];
	int set = 0;
	for (int i=1; i<=p->dim; i++) {
		char LL[256] = "L";
		rewind(config);
		sprintf(istr, "%d", i);
		strcat(LL, istr);

		while(!feof(config)) {

			if(fgets(total,198,config) == NULL) {
				// end of file?
				break;
			}

			ret = sscanf(total,"%99s%99s",key,value);

			if(ret == EOF) {
				continue;
			}

			if(key[0] == '#') {
				continue;
			}

			if(ret != 2)
				fprintf(stderr,"Unable to read input line: %s",total);

			if(!strcasecmp(key, LL)) {
				p->L[i-1] = strtol(value,NULL,10);
				set++;
				continue;
			}
		}
	}

	if (set != p->dim) {
		fprintf(stderr,
		"Failed to set lattice lengths! Only set %d sides.\n", set);
	} else {
		set_L = 1;
	}

	check_set(set_L, "L");


  // calculate total volume now that we know the side lengths
  p->vol = 1;
  for (int i=0; i<p->dim; i++) {
      p->vol *= p->L[i];
			if (p->L[i] % 2 != 0) {
				printf0(*p, "WARNING!! Lattice side length L%d is an odd number! Checkerboard updating is not well defined...\n", i);
			}
  }

	check_set(set_multicanonical, "multicanonical");

  check_set(set_iterations, "iterations");
	check_set(set_interval, "measure interval");
  check_set(set_checkpoint, "checkpoint");
  check_set(set_checks, "run_checks");

  check_set(set_su2alg, "algorithm_su2link");
	check_set(set_betasu2, "betasu2");

	#ifdef HIGGS // now these parameters are attempted to read but the checks are skipped
  check_set(set_su2DBalg, "algorithm_su2doublet");
	check_set(set_update_doublet, "update_doublet");
	check_set(set_phi0, "phi0");
	check_set(set_lambda_phi, "lambda");
	check_set(set_msq_phi, "msq");
	#endif
	#ifdef TRIPLET
	check_set(set_su2triplet_alg, "algorithm_su2triplet");
	check_set(set_update_triplet, "update_triplet");
	check_set(set_sigma0, "sigma0");
	check_set(set_b4, "b4");
	check_set(set_msq_triplet, "msq_triplet");
	#endif
  #if defined TRIPLET && defined DOUBLET
  check_set(set_a2, "a2");
  #endif

  check_set(set_resultsfile, "resultsfile");


  fclose(config);

}


/* Just like get_parameters(), but reads params related to
* multicanonical weighting.
*/
void get_weight_parameters(char *filename, params *p, weight* w) {

	if (!p->multicanonical) {
		// no multicanonical, so just set dummy values
		w->bins = 0;
		w->min = 0;
		w->max = 0;
		w->readonly = 1;
		w->increment = 0;
		strcpy(w->weightfile,"weight");
	} else {
		// multicanonical run, read from config

		int set_bins = 0;
		int set_min = 0, set_max = 0;
		int set_increment = 0;
		int set_readonly = 0;
		int set_weightfile = 0;

    char key[100];
    char value[100];
    char total[200];
    int ret;
    int i;

		FILE *config;

		if(strlen(filename) == 1 && filename[0] == '-') {
			fprintf(stderr,"reading config from stdin\n");
			config = stdin;
		} else {
			if(access(filename,R_OK) != 0) {
				fprintf(stderr,"unable to access config file %s\n", filename);
				exit(101);
			}

			config = fopen(filename, "r");
		}

		while(!feof(config)) {

			if(fgets(total,198,config) == NULL) {
				// end of file?
				break;
			}

			ret = sscanf(total,"%99s%99s",key,value);

			if(ret == EOF) {
				continue;
			}

			if(key[0] == '#') {
				continue;
			}

			if(ret != 2)
				fprintf(stderr,"Unable to read input line: %s",total);

			else if(!strcasecmp(key,"bins")) {
				w->bins = strtol(value,NULL,10);
				set_bins = 1;
			} else if(!strcasecmp(key,"min")) {
				w->min = strtod(value,NULL);
				set_min = 1;
			} else if(!strcasecmp(key,"max")) {
				w->max = strtod(value,NULL);
				set_max = 1;
			} else if(!strcasecmp(key,"increment")) {
				w->increment = strtod(value,NULL);
				set_increment = 1;
			} else if(!strcasecmp(key,"readonly")) {
				w->readonly = strtol(value,NULL,10);
				set_readonly = 1;
			} else if(!strcasecmp(key,"weightfile")) {
				strcpy(w->weightfile,value);
				set_weightfile = 1;
			}

		}

		check_set(set_bins, "bins");
		check_set(set_min, "min");
		check_set(set_max, "max");
		check_set(set_readonly, "readonly");
		check_set(set_weightfile, "weightfile");

		fclose(config);

	}

}

/*
/** Write a summary of the loaded parameters to stderr.
 */

void print_parameters(params p) {
	printf("Volume ");
	for (int i = 0; i<p.dim; i++) {
		printf("%d x ", p.L[i]);
	}
	printf( "\b\b\b, total %lu\n",p.vol);

  printf("iterations %lu, measurement interval %lu, checkpoint interval %lu\n",
		p.iterations, p.interval, p.checkpoint);

	printf("-------------------------- Lattice parameters --------------------------\n");
	printf("SU(2) beta %.1lf\n", p.betasu2);
	#ifdef HIGGS
	printf("msq (Higgs) %lf, lambda (Higgs) %lf, ", p.msq_phi, p.lambda_phi);
	printf("initial phi0 %.2lf\n",p.phi0);
	#endif
	#ifdef TRIPLET
	printf("msq (triplet) %lf, b4 %lf, a2 %lf, ", p.msq_triplet, p.b4, p.a2);
	printf("initial sigma0 %.2lf\n",p.sigma0);
	#endif
	printf("\n");
  //fprintf(stderr, "msq %lf, lambda %lf, g %lf, msigmasq %lf, a2 %lf, b4 %lf\n",
	  //p.msq, p.lambda, p.g, p.msigmasq, p.a2, p.b4);

}

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
 * Original implementation by David Weir.
 */
void get_parameters(char *filename, params *p) {

	int set_dim = 0;
  int set_L = 0;

  int set_reset = 0;
  int set_iterations = 0;
  int set_checkpoint = 0;
	int set_interval = 0;
  int set_n_thermalize = 0;
  int set_checks = 0;

  int set_su2alg = 0;
  int set_u1alg = 0;
  int set_su2DBalg = 0;
	int set_su2triplet_alg = 0;

	int set_multicanonical = 0;

  int set_seed = 0;

  int set_resultsfile = 0;
  int set_latticefile = 0;

	int set_betasu2 = 0;
  int set_betau1 = 0;
	// doublet
	int set_lambda_phi = 0;
	int set_msq_phi = 0;
	int set_phi0 = 0;
	// triplet
	int set_msq_triplet = 0;
	int set_a2 = 0;
	int set_b4 = 0;
	int set_sigma0 = 0;

  int set_update_links = 0;
  int set_scalar_sweeps = 0;
  int set_update_doublet = 0;
	int set_update_triplet = 0;

  #ifdef NUCLEATION
  int set_traj_min = 0;
  int set_traj_max = 0;
  int set_n_traj = 0;
  #endif

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
    else if(!strcasecmp(key,"reset")) {
      p->reset = strtol(value,NULL,10);
      set_reset = 1;
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
    else if(!strcasecmp(key,"n_thermalize")) {
      p->n_thermalize = strtol(value,NULL,10);
      set_n_thermalize = 1;
    }
    else if(!strcasecmp(key,"run_checks")) {
      p->run_checks = strtol(value,NULL,10);
      set_checks = 1;
    }
		else if(!strcasecmp(key,"betasu2")) {
      p->betasu2 = strtod(value,NULL);
      set_betasu2 = 1;
    }
    else if(!strcasecmp(key,"betau1")) {
      p->betau1 = strtod(value,NULL);
      set_betau1 = 1;
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
    else if(!strcasecmp(key,"latticefile")) {
      strcpy(p->latticefile,value);
      set_latticefile = 1;
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
    } else if(!strcasecmp(key,"update_links")) {
      p->update_links = strtol(value,NULL,10);
      set_update_links = 1;
    }
    else if(!strcasecmp(key,"scalar_sweeps")) {
      p->scalar_sweeps = strtol(value,NULL,10);
      set_scalar_sweeps = 1;
    }

    #ifdef U1
      else if(!strcasecmp(key,"algorithm_u1link")) {
        if(!strcasecmp(value,"metropolis")) {
  	       p->algorithm_u1link = METROPOLIS;
        } else if(!strcasecmp(value,"heatbath")) {
          printf("U(1) heatbath not implemented!!!!!!\n");
  	      p->algorithm_u1link = METROPOLIS;
        } else {
          printf("Unknown algorithm for U(1) gauge fields, using metropolis\n");
          p->algorithm_u1link = METROPOLIS;
        }
        set_u1alg = 1;
      } else if(!strcasecmp(key,"update_links")) {
        p->update_links = strtol(value,NULL,10);
        set_update_links = 1;
      }
    #endif

    #ifdef HIGGS
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
      } else if(!strcasecmp(key,"update_doublet")) {
        p->update_su2doublet = strtol(value,NULL,10);
        set_update_doublet = 1;
      }
    #endif

    #ifdef TRIPLET
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
      } else if(!strcasecmp(key,"update_triplet")) {
        p->update_su2triplet = strtol(value,NULL,10);
        set_update_triplet = 1;
      }
    #endif
    // end reading update algorithms

    #ifdef NUCLEATION
      else if(!strcasecmp(key,"traj_min")) {
        traj_min = strtod(value,NULL);
        set_traj_min = 1;
      } else if(!strcasecmp(key,"traj_max")) {
        traj_max = strtod(value,NULL);
        set_traj_max = 1;
      } else if(!strcasecmp(key,"n_traj")) {
        n_traj = strtol(value,NULL);
        set_n_traj = 1;
      }
    #endif

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

  check_set(set_reset, "reset");
  check_set(set_iterations, "iterations");
	check_set(set_interval, "measure interval");
  check_set(set_n_thermalize, "n_thermalize");
  check_set(set_checkpoint, "checkpoint");
  check_set(set_checks, "run_checks");

  check_set(set_su2alg, "algorithm_su2link");
	check_set(set_betasu2, "betasu2");
  check_set(set_update_links, "update_links");
  check_set(set_scalar_sweeps, "scalar_sweeps");
  #ifdef U1
  check_set(set_betau1, "betau1");
  check_set(set_u1alg, "algorithm_u1link");
  #endif
  // these parameters are attempted to read but the checks are skipped
  // if preprocessor identifier not defined
	#ifdef HIGGS
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
  #if defined (TRIPLET) && defined (HIGGS)
  check_set(set_a2, "a2");
  #endif

  // bubble nucleation:
  #ifdef NUCLEATION
    check_set(set_traj_min, "traj_min");
    check_set(set_traj_max, "traj_max");
    check_set(set_n_traj, "n_traj");
  #endif

  check_set(set_resultsfile, "resultsfile");
  check_set(set_latticefile, "latticefile");

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
    w->reduction_factor = 1;
		strcpy(w->weightfile,"weight");
	} else {
		// multicanonical run, read from config

		int set_bins = 0;
		int set_min = 0, set_max = 0;
		int set_increment = 0;
		int set_readonly = 0;
		int set_weightfile = 0;
		int set_absolute_bounds = 0;
    int set_orderparam = 0;
    int set_reduction_factor = 0;

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
			} else if(!strcasecmp(key,"reduction_factor")) {
				w->reduction_factor = strtod(value,NULL);
				set_reduction_factor = 1;
			} else if(!strcasecmp(key,"absolute_bounds")) {
				w->absolute_bounds = strtol(value,NULL, 10);
				set_absolute_bounds = 1;
			} else if(!strcasecmp(key,"readonly")) {
				w->readonly = strtol(value,NULL,10);
				set_readonly = 1;
			} else if(!strcasecmp(key,"weightfile")) {
				strcpy(w->weightfile,value);
				set_weightfile = 1;
			}
      // read multicanonical order parameter
      else if(!strcasecmp(key,"orderparam")) {
				if (!strcasecmp(value,"wilson")) {
					// TODO
					w->orderparam = 0;
				}
				#ifdef HIGGS
        else if(!strcasecmp(value,"phisq")) {
  	       w->orderparam = PHISQ;
           printf0(*p, "Multicanonical order parameter: phi^2\n");
        }
				#endif
				#ifdef TRIPLET
				else if (!strcasecmp(value,"Sigmasq")) {
  	       w->orderparam = SIGMASQ;
           printf0(*p, "Multicanonical order parameter: Tr Sigma^2 \n");
        }
				else if (!strcasecmp(value,"phi2Sigma2")) {
  	       w->orderparam = PHI2SIGMA2;
           printf0(*p, "Multicanonical order parameter: phi^2 Tr Sigma^2 \n");
        }
          #ifdef HIGGS
          else if (!strcasecmp(value,"phi2minusSigma2")) {
    	       w->orderparam = PHI2MINUSSIGMA2;
             printf0(*p, "Multicanonical order parameter: phi^2 - Tr Sigma^2 \n");
          }
          #endif
				#endif
				else {
          printf0(*p, "Unknown multicanonical order parameter!!\n");
          w->orderparam = 0;
        }
      } // end order param

		}
    if (w->orderparam != 0)
      set_orderparam = 1;

		check_set(set_bins, "bins");
		check_set(set_min, "min");
		check_set(set_max, "max");
		check_set(set_readonly, "readonly");
		check_set(set_weightfile, "weightfile");
		check_set(set_absolute_bounds, "absolute_bounds");
    check_set(set_orderparam, "orderparam");
    check_set(set_increment, "increment");
    check_set(set_reduction_factor, "reduction_factor");

		fclose(config);

	}

}

/*
/** Write a summary of the loaded parameters to stderr.
 */

void print_parameters(params p) {
	printf("Volume ");
	for (int dir = 0; dir<p.dim; dir++) {
    if (dir>0) printf(" x ");
		printf("%d", p.L[dir]);
	}
	printf( ", total %lu\n",p.vol);

  printf("iterations %lu, measurement interval %lu, checkpoint interval %lu\n",
		p.iterations, p.interval, p.checkpoint);
  printf("Will perform %d gauge sweeps, %d scalar sweeps per iteration\n",
    p.update_links, p.scalar_sweeps);

	printf("-------------------------- Lattice parameters --------------------------\n");
	printf("SU(2) beta %.1lf\n", p.betasu2);
  #ifdef U1
  printf("U(1) beta %.1lf\n", p.betau1);
  #endif
	#ifdef HIGGS
	printf("msq (Higgs) %lf, lambda (Higgs) %lf, ", p.msq_phi, p.lambda_phi);
	printf("initial phi0 %.2lf, update_su2doublet %d\n",p.phi0, p.update_su2doublet);
	#endif
	#ifdef TRIPLET
	printf("msq (triplet) %lf, b4 %lf, a2 %lf, ", p.msq_triplet, p.b4, p.a2);
	printf("initial sigma0 %.2lf, update_su2triplet %d\n",p.sigma0, p.update_su2triplet);
	#endif
	printf("\n");
  //fprintf(stderr, "msq %lf, lambda %lf, g %lf, msigmasq %lf, a2 %lf, b4 %lf\n",
	  //p.msq, p.lambda, p.g, p.msigmasq, p.a2, p.b4);

}

/* Read config file again and update certain values. This is called at every checkpoint to
* see if the iterations limit has been changed by the user.
*/
void read_updated_parameters(char *filename, params *p) {

	FILE *config;
	long new;
	char key[100];
  char value[100];
	char total[200];
  int ret;

	if(access(filename,R_OK) == 0) {
		config = fopen(filename, "r");

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

			if(!strcasecmp(key,"iterations")) {
				new = strtol(value,NULL,10);
				if (p->iterations != new) {
					p->iterations = new;
					printf0(*p, "Updated max iterations to %ld\n", new);
				}
			}
			else if(!strcasecmp(key,"interval")) {
				new = strtol(value,NULL,10);
				if (p->interval != new) {
					p->interval = new;
					printf0(*p, "Updated measurement interval to %ld\n", new);
				}
			}
			else if(!strcasecmp(key,"checkpoint")) {
				new = strtol(value,NULL,10);
				if (p->checkpoint != new) {
					p->checkpoint = new;
					printf0(*p, "Updated checkpoint interval to %ld\n", new);
				}
			}

		}
		fclose(config);
	}

}

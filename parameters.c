/** @file parameters.c Parses parameters from a file.
*		Original implementation by David Weir.
 */

 // TODO when parallelizing, modify the resultsfile part so that only one node does writing? See David

#include "su2.h"


/** Helper routine, checks a parameter is set and if not prints error.
 */
void check_set(int set, char *name) {
  if(!set) {
    fprintf(stderr,
	    "Not set parameter \"%s\"! Exiting...\n", name);
    die(1);
  }
}

/* quick function for reading couplings etc into a memory address*/
int read_double(char* key, char* name, char* value, double* address) {
  if(!strcasecmp(key, name)) {
    *address = strtod(value,NULL);
    return 1;
  } else return 0;
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
void get_parameters(char *filename, lattice* l, params *p) {


	int set_dim = 0;
  int set_L = 0;

  int set_reset = 0;
  int set_iterations = 0;
  int set_checkpoint = 0;
	int set_interval = 0;
  int set_n_thermalize = 0;
  int set_checks = 0;
  int set_local_meas = 0;

  int set_interval_z = 0;

  #ifdef CORRELATORS
    int set_do_correlators = 0;
    int set_correlator_interval = 0;
  #endif

  #ifdef BLOCKING
    int set_blocks = 0;
  #endif

  #ifdef GRADFLOW
    int set_do_flow = 0;
    int set_flow_dt = 0;
    int set_flow_t_max = 0;
    int set_flow_interval = 0;
    int set_flow_meas_interval = 0;
  #endif

  #ifdef HB_TRAJECTORY
    int set_do_trajectory = 0;
  #endif

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
  int set_r_u1 = 0;
	// doublet
	int set_lambda_phi = 0;
	int set_msq_phi = 0;
	int set_phi0 = 0;
	// triplet
	int set_msq_triplet = 0;
	int set_a2 = 0;
	int set_b4 = 0;
	int set_sigma0 = 0;

  #if (NHIGGS == 2)
    // too many params to write all out, just count if we set them all or not
    int higgs2_params = 14 - 2; // lam1 and m1sq are accounted above
    int set_higgs2_params = 0;
  #endif

  #ifdef SINGLET
    int set_update_singlet = 0;
    int set_singletalg = 0;
    int set_singlet0 = 0;
    int set_b1s=0, set_msq_s=0, set_b3s=0, set_b4s=0;
    #if (NHIGGS == 1)
      int set_a1s=0, set_a2s=0;
    #endif
  #endif

  int set_update_links = 0;
  int set_scalar_sweeps = 0;
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
      l->dim = strtol(value,NULL,10);
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
    #ifdef MEASURE_Z
      else if(!strcasecmp(key,"interval_z")) {
        p->meas_interval_z = strtol(value,NULL,10);
        set_interval_z = 1;
      }
    #endif
    else if(!strcasecmp(key,"checkpoint")) {
      p->checkpoint = strtol(value,NULL,10);
      set_checkpoint = 1;
    } else if(!strcasecmp(key,"measure_local")) {
      p->do_local_meas = strtol(value,NULL,10);
      set_local_meas = 1;
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
    #ifdef SINGLET
      else if(!strcasecmp(key,"update_singlet")) {
        p->update_singlet = strtol(value,NULL,10);
        set_update_singlet = 1;
      }
      else if(!strcasecmp(key,"singlet0")) {
        p->singlet0 = strtod(value,NULL);
        set_singlet0 = 1;
      }
      else if(!strcasecmp(key,"b1_s")) {
        p->b1_s = strtod(value,NULL);
        set_b1s = 1;
      }
      else if(!strcasecmp(key,"msq_s")) {
        p->msq_s = strtod(value,NULL);
        set_msq_s = 1;
      }
      else if(!strcasecmp(key,"b3_s")) {
        p->b3_s = strtod(value,NULL);
        set_b3s = 1;
      }
      else if(!strcasecmp(key,"b4_s")) {
        p->b4_s = strtod(value,NULL);
        set_b4s = 1;
      }
      #if (NHIGGS == 1)
        else if(!strcasecmp(key,"a1_s")) {
          p->a1_s = strtod(value,NULL);
          set_a1s = 1;
        }
        else if(!strcasecmp(key,"a2_s")) {
          p->a2_s = strtod(value,NULL);
          set_a2s = 1;
        }
      #endif
    #endif

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
				if(!l->rank) {
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
      } else if(!strcasecmp(key,"betau1")) {
          p->betau1 = strtod(value,NULL);
          set_betau1 = 1;
      } else if(!strcasecmp(key,"r_u1")) {
          p->r_u1 = strtod(value,NULL);
          set_r_u1 = 1;
      }
    #endif

    #if (NHIGGS > 0)
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
    #ifdef SINGLET
      else if(!strcasecmp(key,"algorithm_singlet")) {
        if(!strcasecmp(value,"metropolis")) {
           p->algorithm_singlet = METROPOLIS;
        } else if(!strcasecmp(value,"overrelax")) {
           p->algorithm_singlet = OVERRELAX;
        } else {
          printf("Unknown algorithm for singlets, using metropolis\n");
          p->algorithm_singlet = METROPOLIS;
        }
        set_singletalg = 1;
      } else if(!strcasecmp(key,"update_singlet")) {
        p->update_singlet = strtol(value,NULL,10);
        set_update_singlet = 1;
      }
    #endif
    // end reading update algorithms

    #ifdef GRADFLOW
    else if(!strcasecmp(key,"do_flow")) {
      p->do_flow = strtol(value,NULL,10);
      set_do_flow = 1;
    } else if(!strcasecmp(key,"flow_meas_interval")) {
      p->flow_meas_interval = strtol(value,NULL,10);
      set_flow_meas_interval = 1;
    } else if(!strcasecmp(key,"flow_interval")) {
      p->flow_interval = strtol(value,NULL,10);
      set_flow_interval = 1;
    } else if(!strcasecmp(key,"flow_dt")) {
      p->flow_dt = strtod(value,NULL);
      set_flow_dt = 1;
    } else if(!strcasecmp(key,"flow_t_max")) {
      p->flow_t_max = strtod(value,NULL);
      set_flow_t_max = 1;
    }
    #endif

    #ifdef CORRELATORS
      else if(!strcasecmp(key,"do_correlators")) {
        p->do_correlators = strtol(value,NULL,10);
        set_do_correlators = 1;
      } else if(!strcasecmp(key,"correlator_interval")) {
        p->correlator_interval = strtol(value,NULL,10);
        set_correlator_interval = 1;
      }
    #endif

    #ifdef BLOCKING
      else if(!strcasecmp(key,"blocks")) {
        p->blocks = strtol(value,NULL,10);
        set_blocks = 1;
      }
    #endif

    #ifdef HB_TRAJECTORY
      else if(!strcasecmp(key,"do_trajectory")) {
        p->do_trajectory = strtol(value,NULL,10);
        set_do_trajectory = 1;
      }
    #endif

    #if (NHIGGS == 2)
      else {
        // for two Higgs potential parameters use the shorthand routine
        set_higgs2_params += read_double(key, "msq_phi2", value, &p->msq_phi2);
        set_higgs2_params += read_double(key, "m12sq_re", value, &p->m12sq.re);
        set_higgs2_params += read_double(key, "m12sq_im", value, &p->m12sq.im);
        set_higgs2_params += read_double(key, "lam2", value, &p->lam2);
        set_higgs2_params += read_double(key, "lam3", value, &p->lam3);
        set_higgs2_params += read_double(key, "lam4", value, &p->lam4);
        set_higgs2_params += read_double(key, "lam5_re", value, &p->lam5.re);
        set_higgs2_params += read_double(key, "lam5_im", value, &p->lam5.im);
        set_higgs2_params += read_double(key, "lam6_re", value, &p->lam6.re);
        set_higgs2_params += read_double(key, "lam6_im", value, &p->lam6.im);
        set_higgs2_params += read_double(key, "lam7_re", value, &p->lam7.re);
        set_higgs2_params += read_double(key, "lam7_im", value, &p->lam7.im);
      }
    #endif

  }

	check_set(set_dim, "dim");

	// Allocate memory for lattice side lengths.
	// This couldn't be done earlier because we didn't know the dimension
	// p.parity is allocated in alloc.c, alloc_neighbors().
	l->L = malloc(l->dim * sizeof(*(l->L)));

	// read the file again for lattice sizes L_i. This could use some polishing
	char istr[255];
	int set = 0;
	for (int i=1; i<=l->dim; i++) {
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
				l->L[i-1] = strtol(value,NULL,10);
				set++;
				continue;
			}
		}
	}

	if (set != l->dim) {
		fprintf(stderr,
		"Failed to set lattice lengths! Only set %d sides.\n", set);
	} else {
		set_L = 1;
	}

	check_set(set_L, "L");


  // calculate total volume now that we know the side lengths
  l->vol = 1;
  for (int dir=0; dir<l->dim; dir++) {
    l->vol *= l->L[dir];
		if (l->L[dir] % 2 != 0) {
			printf0(*l, "WARNING!! Lattice side length L%d is an odd number! Checkerboard updating is not well defined...\n", dir);
		}
  }

	check_set(set_multicanonical, "multicanonical");

  check_set(set_reset, "reset");
  check_set(set_iterations, "iterations");
	check_set(set_interval, "measure interval");
  check_set(set_n_thermalize, "n_thermalize");
  check_set(set_checkpoint, "checkpoint");
  check_set(set_checks, "run_checks");
  check_set(set_local_meas, "measure_local");

  check_set(set_su2alg, "algorithm_su2link");
	check_set(set_betasu2, "betasu2");
  check_set(set_update_links, "update_links");
  check_set(set_scalar_sweeps, "scalar_sweeps");
  #ifdef U1
  check_set(set_betau1, "betau1");
  check_set(set_r_u1, "gammau1");
  check_set(set_u1alg, "algorithm_u1link");
  #endif
  // these parameters are attempted to read but the checks are skipped
  // if preprocessor identifier not defined
	#if (NHIGGS > 0)
  check_set(set_su2DBalg, "algorithm_su2doublet");
	check_set(set_update_doublet, "update_doublet");
	check_set(set_phi0, "phi0");
	check_set(set_lambda_phi, "lambda");
	check_set(set_msq_phi, "msq");
	#endif

  #if (NHIGGS == 2)
    if (set_higgs2_params != higgs2_params) {
      printf("Error setting parameters for second Higgs! Got %d, expected %d\n", set_higgs2_params, higgs2_params);
      die(-1131);
    }
  #endif

	#ifdef TRIPLET
	check_set(set_su2triplet_alg, "algorithm_su2triplet");
	check_set(set_update_triplet, "update_triplet");
	check_set(set_sigma0, "sigma0");
	check_set(set_b4, "b4");
	check_set(set_msq_triplet, "msq_triplet");
	#endif
  #if defined (TRIPLET) && (NHIGGS > 0)
  check_set(set_a2, "a2");
  #endif

  #ifdef SINGLET
  check_set(set_singletalg, "algorithm_singlet");
  check_set(set_update_singlet, "update_singlet");
  check_set(set_singlet0, "singlet0");
  check_set(set_b1s, "b1_s");
  check_set(set_b3s, "b3_s");
  check_set(set_b4s, "b4_s");
  check_set(set_msq_s, "msq_s");
  #if (NHIGGS==1)
    check_set(set_a1s, "a1_s");
    check_set(set_a2s, "a2_s");
  #endif
  #endif

  #ifdef MEASURE_Z
    check_set(set_interval_z, "interval_z");
  #endif

  check_set(set_resultsfile, "resultsfile");
  check_set(set_latticefile, "latticefile");

  #ifdef GRADFLOW
    check_set(set_do_flow, "do_flow");
    check_set(set_flow_meas_interval, "flow_meas_interval");
    check_set(set_flow_interval, "flow_interval");
    check_set(set_flow_dt, "flow_dt");
    check_set(set_flow_t_max, "flow_t_max");
  #endif

  #ifdef CORRELATORS
    check_set(set_do_correlators, "do_correlators");
    check_set(set_correlator_interval, "correlator_interval");
  #endif

  #ifdef BLOCKING
    check_set(set_blocks, "blocks");
  #endif

  #ifdef HB_TRAJECTORY
    check_set(set_do_trajectory, "do_trajectory");
  #endif

  fclose(config);

}


/* Just like get_parameters(), but reads params related to
* multicanonical weighting.
*/
void get_weight_parameters(char *filename, lattice const* l, params *p, weight* w) {

	if (!p->multicanonical) {
		// no multicanonical, so just set dummy values
		w->bins = 0;
		w->min = 0;
		w->max = 0;
		w->mode = 0;
		w->delta = 0;
    w->do_acceptance = 0;
    w->orderparam = -1;
		strcpy(w->weightfile,"weight");
	} else {
		// multicanonical run, read from config

		int set_bins = 0;
		int set_min = 0, set_max = 0;
		int set_muca_delta = 0;
		int set_mode = 0;
		int set_weightfile = 0;
    int set_orderparam = 0;
    int set_checks_per_sweep = 0;

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
			} else if(!strcasecmp(key,"muca_delta")) {
				w->delta = strtod(value,NULL);
				set_muca_delta = 1;
			} else if(!strcasecmp(key,"muca_mode")) {
				w->mode = strtol(value,NULL,10);
				set_mode = 1;
			} else if(!strcasecmp(key,"checks_per_sweep")) {
				w->checks_per_sweep = strtol(value,NULL,10);
				set_checks_per_sweep = 1;
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
				#if (NHIGGS > 0)
        else if(!strcasecmp(value,"phisq")) {
  	       w->orderparam = PHISQ;
           printf0(*l, "Multicanonical order parameter: phi^2\n");
        }
				#endif
        #if (NHIGGS > 1)
        else if(!strcasecmp(value,"phi2sq")) {
  	       w->orderparam = PHI2SQ;
           printf0(*l, "Multicanonical order parameter: phi2^2\n");
        }
				#endif
				#ifdef TRIPLET
				else if (!strcasecmp(value,"Sigmasq")) {
  	       w->orderparam = SIGMASQ;
           printf0(*l, "Multicanonical order parameter: Tr Sigma^2 \n");
        }
          #if (NHIGGS > 0)
          else if (!strcasecmp(value,"phi2minusSigma2")) {
    	       w->orderparam = PHI2MINUSSIGMA2;
             printf0(*l, "Multicanonical order parameter: phi^2 - Tr Sigma^2 \n");
          }
          #endif
				#endif
				else {
          printf0(*l, "Unknown multicanonical order parameter!!\n");
          w->orderparam = 0;
        }
      } // end order param

		}
    if (w->orderparam != 0)
      set_orderparam = 1;

		check_set(set_bins, "bins");
		check_set(set_min, "min");
		check_set(set_max, "max");
		check_set(set_mode, "muca_mode");
		check_set(set_weightfile, "weightfile");
    check_set(set_orderparam, "orderparam");
    check_set(set_muca_delta, "muca_delta");
    check_set(set_checks_per_sweep, "checks_per_sweep");

    if (w->mode != READONLY && w->mode != FAST && w->mode != SLOW) {
      printf0(*l, "Invalid mode specification for multicanonical!! got mode=%d\n", w->mode);
      printf0(*l, "%d = read only, %d = fast weight update, %d = slow weight update\n", READONLY, FAST, SLOW);
      die(551);
    }

		fclose(config);
	}

}

/*
/** Write a summary of the loaded parameters to stderr.
 */

void print_parameters(lattice l, params p) {
	printf("Volume ");
	for (int dir = 0; dir<l.dim; dir++) {
    if (dir>0) printf(" x ");
		printf("%d", l.L[dir]);
	}
	printf( ", total %lu\n", l.vol);

  printf("iterations %lu, measurement interval %lu, checkpoint interval %lu\n",
		p.iterations, p.interval, p.checkpoint);
  printf("Will perform %d gauge sweeps, %d scalar sweeps per iteration\n",
    p.update_links, p.scalar_sweeps);

	printf("-------------------------- Lattice parameters --------------------------\n");
	printf("SU(2) beta %g\n", p.betasu2);
  #ifdef U1
    printf("U(1) beta %g, U(1) r %g\n", p.betau1, p.r_u1);
  #endif

  #if (NHIGGS == 2)
    printf("msq1 %g, msq2 %g, m12sq_re %g, m12sq_im %g \n", p.msq_phi, p.msq_phi2, p.m12sq.re, p.m12sq.im);
    printf("lam1 %g, lam2 %g, lam3 %g, lam4 %g, lam5_re %g, lam5_im %g \n", p.lambda_phi, p.lam2, p.lam3, p.lam4, p.lam5.re, p.lam5.im);
    printf("lam6_re %g, lam6_im %g, lam7_re %g, lam7_im %g \n", p.lam6.re, p.lam6.im, p.lam7.re, p.lam7.im);
    printf("initial phi0 %g, update_su2doublet %d\n",p.phi0, p.update_su2doublet);

  #elif (NHIGGS == 1)
	  printf("msq (Higgs) %g, lambda (Higgs) %g, ", p.msq_phi, p.lambda_phi);
	  printf("initial phi0 %g, update_su2doublet %d\n",p.phi0, p.update_su2doublet);
	#endif

	#ifdef TRIPLET
	printf("msq (triplet) %g, b4 %g, ", p.msq_triplet, p.b4);
    #if (NHIGGS > 0)
    printf("a2 %g, ", p.a2);
    #endif
	printf("initial sigma0 %g, update_su2triplet %d\n",p.sigma0, p.update_su2triplet);
	#endif

  #ifdef SINGLET
  printf("msq (singlet) %g, b1_s %g, b3_s %g, b4_s %g, ", p.msq_s, p.b1_s, p.b3_s, p.b4_s);
    #if (NHIGGS > 0)
    printf("a1_s %g, a2 %g, ", p.a1_s, p.a2_s);
    #endif
	printf("initial singlet0 %g, update_singlet %d\n",p.singlet0, p.update_singlet);
  #endif

	printf("\n");
  //fprintf(stderr, "msq %lf, lambda %lf, g %lf, msigmasq %lf, a2 %lf, b4 %lf\n",
	  //p.msq, p.lambda, p.g, p.msigmasq, p.a2, p.b4);

}

/* Read config file again and update certain values. This is called at every checkpoint to
* see if the iterations limit has been changed by the user.
*/
void read_updated_parameters(char *filename, lattice const* l, params *p, weight* w) {

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
					printf0(*l, "Updated max iterations to %ld\n", new);
				}
			}
			else if(!strcasecmp(key,"interval")) {
				new = strtol(value,NULL,10);
				if (p->interval != new) {
					p->interval = new;
					printf0(*l, "Updated measurement interval to %ld\n", new);
				}
			}
			else if(!strcasecmp(key,"checkpoint")) {
				new = strtol(value,NULL,10);
				if (p->checkpoint != new) {
					p->checkpoint = new;
					printf0(*l, "Updated checkpoint interval to %ld\n", new);
				}
			}
      else if(!strcasecmp(key,"checks_per_sweep") && p->multicanonical) {
				new = strtol(value,NULL,10);
				if (w->checks_per_sweep != new) {
					w->checks_per_sweep = new;
					printf0(*l, "Updated multicanonical checks per sweep to %ld\n", new);
				}
			}

		}
		fclose(config);
	}

}

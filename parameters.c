/** @file parameters.c. Reads parameters from a file.
*/


#include "su2.h"

const int maxLineLen = 200;

/* Open file for reading in root node only. 
* Return value is 1 if access is OK, 0 otherwise.
* Need to pass address of the FILE* pointer that we want to open. */
int OpenRead(char* fname, FILE** file) {
  int success = 1;
  if (!myRank) {
    if(access(fname, R_OK) != 0) {
      printf("!!! Unable to access file %s\n", fname);
      success = 0;
    } else {
      // Note syntax: need to change the real pointer, not the pointer pointing to it
      *file = fopen(fname, "r"); 
    }
  }
  bcast_int(&success, MPI_COMM_WORLD);
  return success;
}

/* Read file line-by-line until a matching label is found.
* Assumes that the lines to be read are of form "label value".
* The resulting value (string) is stored in 'result'.  */
void FindFromFile(FILE* fileIn, char* label, char* result) {

  strcpy(result, "");
  int ok = 1;
  int ret;
  char line[maxLineLen];
  char key[maxLineLen];
  char value[maxLineLen];
  int isEOF = 0;

  // only read if file is already open
  if (!myRank && fileIn == NULL) {
    printf("!!! ERROR: file not open, cannot read...\n");
    ok = 0; 
  } 
  bcast_int(&ok, MPI_COMM_WORLD);
  if (!ok) die(3);

  if (!myRank) {

    while(!feof(fileIn)) {

      // Get next line. If we hit EOF, rewind back to file start,
      // unless we already did that when searching for this parameter.
      if(fgets(line, maxLineLen, fileIn) == NULL) {
        // end of file?
        if (isEOF > 0) {
          break;
        } else {
          isEOF = 1;
          rewind(fileIn);
          continue;
        }
      }

      ret = sscanf(line, "%99s%99s", key, value);
      if(ret == EOF) continue; // No match 
      if(key[0] == '#') continue; // Ignore commented lines
      /* // Check that line is in correct format. 
      if(ret != 2) fprintf("!!! Unable to read input line: %s", line);
      */

      if ( !strcasecmp(key, label) ) {
        // Found matching line
        strcpy(result, value);
        break;
      }

    }

  } else if (myRank == 0) {
    printf("!!! File not open in root node\n");
    ok = 0;
  }

  if ( !strcasecmp(result, "") && myRank == 0) {
    printf("!!! Unable to find parameter %s in file\n", label);
    ok = 0;
  }
  bcast_int(&ok, MPI_COMM_WORLD);
  if (!ok) die(1);

}

// Get integer parameter from file (in root node, and broadcast to others)
int GetInt(FILE* fileIn, char* label) {
  char value[maxLineLen];
  int res;
  FindFromFile(fileIn, label, value);
  res = strtol(value, NULL, 10);
  bcast_int(&res, MPI_COMM_WORLD);
  return res;
}

// Same as getInt(), but for long int type
long GetLong(FILE* fileIn, char* label) {
  char value[maxLineLen];
  long res;
  FindFromFile(fileIn, label, value);
  res = strtol(value, NULL, 10);
  bcast_long(&res, MPI_COMM_WORLD);
  return res;
}

// Same as GetInt(), but for double type
double GetDouble(FILE* fileIn, char* label) {
  char value[maxLineLen];
  double res;
  FindFromFile(fileIn, label, value);
  res = strtod(value, NULL);
  bcast_double(&res, MPI_COMM_WORLD);
  return res;
}

// Same as GetInt(), but for strings (note: no return value)
void GetString(FILE* fileIn, char* label, char* result) {
  FindFromFile(fileIn, label, result);
  bcast_string(result, maxLineLen, MPI_COMM_WORLD);
}

/* Read what update algorithms to use for the fields.
* Return value is integer: METROPOLIS, HEATBATH or OVERRELAX */
int GetUpdateAlgorithm(FILE* fileIn, char* label) {

  int res;
  char alg[maxLineLen];
  GetString(fileIn, label, alg);

  if(!strcasecmp(alg, "metropolis")) {
    res = METROPOLIS;
  } else if(!strcasecmp(alg, "heatbath")) {
    res = HEATBATH;
  } else if(!strcasecmp(alg, "overrelax")) {
    res = OVERRELAX;
  } else {
    if (!myRank) printf("WARNING: Unknown update algorithm for '%s'; using metropolis instead.", label);
    res = METROPOLIS;
  }
  return res;
}


/* Get all parameters from config file. 
* Each parameter is a key-value pair separated by whitespace. A line
* beginning with # is treated as a comment and ignored.
*
* Do NOT call this more than once at runtime!! */
void get_parameters(char *filename, lattice* l, params *p) {

  FILE* config = NULL;
  // Open config file, root node only
  int ok = OpenRead(filename, &config);
  if (!ok) die(1);

  // Start reading params. First lattice dimension:
  l->dim = GetInt(config, "dim");
  // Then alloc the side length array, read in L1, L2 etc and calculate total volume
  l->L = malloc(l->dim * sizeof(*(l->L)));
  l->vol = 1;
  for (int dir=0; dir<l->dim; dir++) {
    char sideLabel[100];
    sprintf(sideLabel, "L%d", dir+1);
    l->L[dir] = GetInt(config, sideLabel);

    l->vol *= l->L[dir];
		if (l->L[dir] % 2 != 0) {
			printf0("!!! WARNING: Lattice side length L%d is an odd number! Checkerboard updating is not well defined\n", dir);
		}
  }

  // Open results file for the root node only
  p->resultsfile = NULL;
  ok = 1;
  char resFileName[maxLineLen];
  GetString(config, "resultsfile", resFileName);
  if(!myRank) p->resultsfile = fopen(resFileName, "a");


  // Read update algorithms to use for fields
  int alg;

  alg = GetUpdateAlgorithm(config, "algorithm_su2link");
  if ( !(alg == METROPOLIS || alg == HEATBATH) ) {
    printf0("Cannot use algorithm %d for SU2 links, changing to metropolis...", alg);
    alg = METROPOLIS;
  }
  p->algorithm_su2link = alg;

  #ifdef U1
    alg = GetUpdateAlgorithm(config, "algorithm_u1link");
    if (alg != METROPOLIS) {
      printf0("Cannot use algorithm %d for U1 links, changing to metropolis...", alg);
      alg = METROPOLIS;
    }
    p->algorithm_u1link = alg;
  #endif

  #if (NHIGGS > 0)
    alg = GetUpdateAlgorithm(config, "algorithm_su2doublet");
    if ( !(alg == METROPOLIS || alg == OVERRELAX) ) {
      printf0("Cannot use algorithm %d for doublets, changing to metropolis...", alg);
      alg = METROPOLIS;
    }
    p->algorithm_su2doublet = alg;
  #endif

  #ifdef TRIPLET
    alg = GetUpdateAlgorithm(config, "algorithm_su2triplet");
    if ( !(alg == METROPOLIS || alg == OVERRELAX) ) {
      printf0("Cannot use algorithm %d for triplets, changing to metropolis...", alg);
      alg = METROPOLIS;
    }
    p->algorithm_su2triplet = alg;
  #endif

  #ifdef SINGLET
    alg = GetUpdateAlgorithm(config, "algorithm_singlet");
    if ( !(alg == METROPOLIS || alg == OVERRELAX) ) {
      printf0("Cannot use algorithm %d for singlets, changing to metropolis...", alg);
      alg = METROPOLIS;
    }
    p->algorithm_singlet = alg;
  #endif
  // end reading update algorithms


  // Everything else from config file:
  GetString(config, "latticefile", p->latticefile); // open only in checkpoint.c

  p->iterations = GetLong(config, "iterations");
  p->interval = GetLong(config, "interval");
  p->checkpoint = GetLong(config, "checkpoint");
  p->n_thermalize = GetLong(config, "n_thermalize");

  p->run_checks = GetInt(config, "run_checks");
  p->reset = GetInt(config, "reset");
  p->multicanonical = GetInt(config, "multicanonical");
  p->do_local_meas = GetInt(config, "measure_local");

  #ifdef MEASURE_Z
    p->meas_interval_z = GetLong(config, "interval_z");
  #endif

  p->update_links = GetInt(config, "update_links");
  p->update_su2doublet = GetInt(config, "update_doublet");
  p->scalar_sweeps = GetInt(config, "scalar_sweeps");

  p->betasu2 = GetDouble(config, "betasu2");
  #if (NHIGGS > 0)
    p->msq_phi = GetDouble(config, "msq");
    p->lambda_phi = GetDouble(config, "lambda");
    p->phi0 = GetDouble(config, "phi0");
  #endif

  #ifdef U1
    p->betau1 = GetDouble(config, "betau1");
    p->r_u1 = GetDouble(config, "r_u1");
  #endif

  // Triplet parameters
  #ifdef TRIPLET
    p->update_su2triplet = GetInt(config, "update_triplet");
    p->sigma0 = GetDouble(config, "sigma0");
    p->msq_triplet = GetDouble(config, "msq_triplet");
    p->b4 = GetDouble(config, "b4");
    #if (NHIGGS > 0)
      p->a2 = GetDouble(config, "a2");
    #endif
  #endif

  // Singlet parameters
  #ifdef SINGLET
    p->update_singlet = GetInt(config, "update_singlet");
    p->singlet0 = GetDouble(config, "singlet0");
    p->b1_s = GetDouble(config, "b1_s");
    p->msq_s = GetDouble(config, "msq_s");
    p->b3_s = GetDouble(config, "b3_s");
    p->b4_s = GetDouble(config, "b4_s");
    #if (NHIGGS > 0)
      p->a1_s = GetDouble(config, "a1_s");
      p->a2_s = GetDouble(config, "a2_s");
    #endif
  #endif

  // 2HDM potential parameters
  #if (NHIGGS == 2)
    p->msq_phi2 = GetDouble(config, "msq_phi2");
    p->m12sq.re = GetDouble(config, "m12sq_re");
    p->m12sq.im = GetDouble(config, "m12sq_im");
    p->lam2 = GetDouble(config, "lam2");
    p->lam3 = GetDouble(config, "lam3");
    p->lam4 = GetDouble(config, "lam4");
    p->lam5.re = GetDouble(config, "lam5_re");
    p->lam5.im = GetDouble(config, "lam5_im");
    p->lam6.re = GetDouble(config, "lam6_re");
    p->lam6.im = GetDouble(config, "lam6_im");
    p->lam7.re = GetDouble(config, "lam7_re");
    p->lam7.im = GetDouble(config, "lam7_im");
  #endif

  #ifdef CORRELATORS
    p->do_correlators = GetInt(config, "do_correlators");
    p->correlator_interval = GetInt(config, "correlator_interval");
  #endif

  #ifdef BLOCKING
    p->blocks = GetInt(config, "blocks");
  #endif

  #ifdef HB_TRAJECTORY
    p->do_trajectory = GetInt(config, "do_trajectory");
  #endif

  #ifdef GRADFLOW
    p->do_flow = GetInt(config, "do_flow");
    p->flow_interval = GetInt(config, "flow_interval");
    p->flow_meas_interval = GetInt(config, "flow_meas_interval");
    p->flow_dt = GetDouble(config, "flow_dt");
    p->flow_t_max = GetDouble(config, "flow_t_max");
  #endif

  #ifdef MEASURE_Z
    p->do_z_meas = GetInt(config, "do_z_meas");
    p->setup_wall = GetInt(config, "setup_wall");
  #endif

  // Reads done!
  if (myRank == 0) fclose(config);
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
		w->mode = 0;
		w->delta = 0;
    w->do_acceptance = 0;
    w->orderparam = -1;
		strcpy(w->weightfile,"weight");
	} else {
		// multicanonical run, read from config
    FILE* config = NULL;
    int ok = OpenRead(filename, &config);
    if (!ok) die(2);

    GetString(config, "weightfile", w->weightfile);

    w->mode = GetInt(config, "muca_mode");
    w->checks_per_sweep = GetInt(config, "checks_per_sweep");
    w->bins = GetInt(config, "bins");

    w->min = GetDouble(config, "min");
    w->max = GetDouble(config, "max");
    w->delta = GetDouble(config, "muca_delta");
    
    // Read order parameter
    char paramName[maxLineLen];
    int res;
    GetString(config, "orderparam", paramName);

    if(!strcasecmp(paramName, "phisq")) {
      res = PHISQ;
    } else if(!strcasecmp(paramName, "Sigmasq")) {
      res = SIGMASQ;
    } else if(!strcasecmp(paramName, "phi2sq")) {
      res = PHI2SQ;
    } else if(!strcasecmp(paramName, "phi2minusSigma2")) {
      res = PHI2MINUSSIGMA2;
    } else {
      printf0("WARNING: Unknown multicanonical order parameter!! using phisq");
      strcpy(paramName, "phisq");
      res = PHISQ;
    }
    w->orderparam = res;


    if (w->mode != READONLY && w->mode != FAST && w->mode != SLOW) {
      printf0("Invalid mode specification for multicanonical!! got mode=%d\n", w->mode);
      printf0("%d = read only, %d = fast weight update, %d = slow weight update\n", READONLY, FAST, SLOW);
      die(551);
    }

		if (!myRank) fclose(config);
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
}


/* Read config file again and update certain values. This is called at every checkpoint to
* see if the iterations limit has been changed by the user. */
void read_updated_parameters(char *filename, lattice const* l, params *p) {
  
	long new;
	char key[maxLineLen];
  char value[maxLineLen];
	char total[maxLineLen];

  FILE* config = NULL;
  // Open config file, root node only
  int ok = OpenRead(filename, &config);
  if (!ok) {
    printf("!!! Warning: Unable to access config file %s at checkpoint; skipping read...\n", filename);
    return;
  }

  new = GetLong(config, "iterations");
  if (p->iterations != new) {
    p->iterations = new;
    printf0("Updated max iterations to %ld\n", new);
  }

  new = GetLong(config, "interval");
  if (p->interval != new) {
    p->interval = new;
    printf0("Updated measurement interval to %ld\n", new);
  }

  new = GetLong(config, "checkpoint");
  if (p->checkpoint != new) {
    p->checkpoint = new;
    printf0("Updated checkpoint interval to %ld\n", new);
  }
  
  if (!myRank) fclose(config);
}
			

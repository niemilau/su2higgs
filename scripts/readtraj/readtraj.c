

#include "stuff.h"

void die(int howbad) {
  exit(howbad);
}


/* 'sample' tells how often to skip measurements,
* with sample = 1 use all t, with sample = 2 use every other t etc.
* blocks = number of jackknife blocks to use */
void realtime(char* fname, int sample, int blocks) {

	if(access(fname, R_OK) != 0) {
		printf("Failed to access file %s !!\n", fname);
		die(1);
	}

	FILE *f = fopen(fname, "r");

	char line[255];

  long traj_tot = 0;

	double min = 0.0, max = 0.0;
	int traj_per_set = 0;
	int read;
  double eps, phi_c;

	// read header
	read = fscanf(f, "%lf %lf %lf %lf %d", &phi_c, &eps, &min, &max, &traj_per_set);
	if (read != 5) {
		printf("Failed to read header!!\n");
		die(1);
	} else {
		printf("Using phi_c = %g, eps = %g, min = %g, max = %g, traj_per_set = %d, sample = %d jackknife blocks = %d\n",
      phi_c, eps, min, max, traj_per_set, sample, blocks);
	}

	// first count the number of trajectory sets
	int nset = 0;
	while (fgets(line, sizeof(line), f)) {
		if (strstr(line, "Begin set ") != NULL) {
			nset++;
		}
	}
	fclose(f);

  // realtime quantities for each set
  double dphidt[nset];
  double dphidt_err[nset];
  double prod[nset]; // d * <|dphi/dt|>
  double prod_err[nset];

  double d[nset]; // "dynamical prefactor"
  double derr[nset];
  for (int i=0; i<nset; i++) {
    derr[i] = 0.0; d[i] = 0.0;
    dphidt_err[i] = 0.0; dphidt[i] = 0.0;
    prod[i] = 0.0; prod_err[i] = 0.0;
  }

	for (int i=0; i<nset; i++) {
    traj_set set;
		make_traj_set(&set, traj_per_set);
		read_traj_set(fname, &set, i+1);

    // now glue the half-trajectories in the set together to form full trajectories
    traj_set set_glue;
    // how many combined trajectories we expect (binomial coefficient):
    long combinations;
    if (set.ntraj < 2) {
      printf("!!! Warning: set contains less than 2 trajectories\n");
      combinations = 1;
    } else {
      //combinations = factorial(set.ntraj) / (factorial(set.ntraj-2) * factorial(2));
      combinations = (set.ntraj * (set.ntraj-1)) / 2;
    }

    make_traj_set(&set_glue, combinations);

    long k = 0;
    for (int j1=0; j1<set.ntraj; j1++) {
      traj* t1 = &set.traj[j1];
      for (int j2=j1+1; j2<set.ntraj; j2++) {
        traj* t2 = &set.traj[j2];

        /* now glue t1 and t2 together:
        * For t1, invert time and skip first measurement,
        * then store in t_new in the inverse order,
        * Depending on the "sampling rate", skip
        * some measurements but always include end and start time */
        traj* t_new = &set_glue.traj[k];

        /* t1->len is the number of time ticks in t1,
        * so the actual "length" is L = t1->len - 1 */

        int L_join = t1->len + t2->len - 2;
        // number of ticks in the glued trajectory:
        int len_new = 1 + L_join / sample; // number of bins in the range + 1 for the last tick
        // one extra bin of different length if the interval didn't divide evenly
        if (L_join % sample != 0) len_new++;

        make_traj(t_new, len_new);

        int cc = 0;
        int l;
        for (l=1; l<t1->len; l += sample) {
          t_new->t[cc] = -1.0*t1->t[t1->len - l];
          t_new->val[cc] = t1->val[t1->len - l];
          cc++;
          //printf("%g %g\n", t_new->t[cc-1], t_new->val[cc-1]);
        }
        for (; l < t1->len + t2->len; l += sample) {
          t_new->t[cc] = t2->t[l - t1->len];
          t_new->val[cc] = t2->val[l - t1->len];
          cc++;
          //printf("%g %g\n", t_new->t[cc-1], t_new->val[cc-1]);
        }
        // add last measurement if it wasn't done above
        if (cc < len_new) {
          t_new->t[len_new-1] = t2->t[t2->len-1];
          t_new->val[len_new-1] = t2->val[t2->len-1];
          cc++;
          //printf("%g %g\n", t_new->t[cc-1], t_new->val[cc-1]);
        }
        if (cc != len_new) printf("!!! Warning: trajectory join (cc=%d, len_new=%d)\n", cc, len_new);

        /*
        for (int l=1; l<t1->len; l++) {
          t_new->t[l-1] = -1.0*t1->t[t1->len - l];
          t_new->val[l-1] = t1->val[t1->len - l];

           //printf("%g %g\n", t_new->t[l-1], t_new->val[l-1]);

        }
        for (int l=0; l<t2->len; l++) {
          t_new->t[l + t1->len-1] = t2->t[l];
          t_new->val[l + t1->len-1] = t2->val[l];

           //printf("%g %g\n", t_new->t[l+t1->len-1], t_new->val[l+t1->len-1]);
        }
        */
        k++;

      }
    } // all combinations done

    if (k != combinations) {
      printf("!!! Error combining trajectories! did %ld, should be %ld\n", k, combinations);
      die(4);
    }

    // realtime quantities for this set:
    /*
    dphidt[i] = dphi_by_dt(&set_glue, phi_c, eps);
    d[i] = dynamical_prefactor(&set_glue, min, max, phi_c);
    prod[i] = d[i] * dphidt[i];
    */

    double res[2];
    dynamics(&set_glue, min, max, phi_c, res);
    dphidt[i] = res[0];
    d[i] = res[1];
    prod[i] = d[i] * dphidt[i];

    traj_tot += k;

    free_traj_set(&set);

    int do_jack = 0;
    if (blocks >= 1 && blocks <= combinations)  {
      /* jackknife */
      do_jack = 1;
      int blen = combinations / blocks;
      int blen_last = blen;
      // if the data doesn't split evenly, last block may have different length from the others
      if (combinations % blocks != 0) {
        blen_last = blen + combinations % blocks;
        //printf("!! Warning: last jackknife block has length %d, others have %d \n", blen_last, blen);
      }

      for (int b_id=0; b_id<blocks; b_id++) {

        // copy all trajectories except those in block 'b_id'
        traj_set set_jack;
        int jack_len;
        if (b_id == blocks-1) jack_len = (blocks-1)*blen;
        else jack_len = (blocks-2)*blen + blen_last;

        make_traj_set(&set_jack, jack_len);

        int cp_tot = 0;
        int skip;
        for (int tr=0; tr<combinations; tr++) {
          // skip trajectories in block 'b_id'
          int end = blen*(b_id+1)-1;
          if (b_id == blocks-1) end = combinations; // special treatment for last block
          if (tr >= blen*b_id && tr <= end) skip = 1;
          else skip = 0;

          if (!skip) {
            traj* cp_tr = &set_glue.traj[tr];
            traj* new_tr = &set_jack.traj[cp_tot];
            make_traj(new_tr, cp_tr->len);
            memcpy(new_tr->t, cp_tr->t, cp_tr->len * sizeof(*cp_tr->t));
            memcpy(new_tr->val, cp_tr->val, cp_tr->len * sizeof(*cp_tr->val));
            cp_tot++;
          }
        }
        if (cp_tot != jack_len) printf("!!! Warning: error making jackknife blocks (jack_len=%d, but collected %d)!\n", jack_len, cp_tot);

        // copy done, now calc jackknife averages for this block

        double res[2];
        dynamics(&set_jack, min, max, phi_c, res);
        double der_jack = res[0];
        double djack = res[1];

        dphidt_err[i] += (dphidt[i] - der_jack) * (dphidt[i] - der_jack);
        derr[i] += (d[i] - djack) * (d[i] - djack);

        double prod_jack = djack * der_jack;
        prod_err[i] += (prod_jack - prod[i]) * (prod_jack - prod[i]);

        free_traj_set(&set_jack);
      }
    } // end jackknife
    else {
      derr[i] = 0.0;
      dphidt_err[i] = 0.0; prod_err[i] = 0.0;
      printf("!!! Warning: no jackknife (requested %d blocks but data length is %ld)\n", blocks, combinations);
    }

    if (do_jack) {
      dphidt_err[i] *= (blocks-1.0) / ((double) blocks);
      dphidt_err[i] = sqrt(dphidt_err[i]);

      derr[i] *= (blocks-1.0) / ((double) blocks);
      derr[i] = sqrt(derr[i]);

      prod_err[i] *= (blocks-1.0) / ((double) blocks);
      prod_err[i] = sqrt(prod_err[i]);
    }

    //printf("d = %g +/- %g\n", d[i], derr[i]);

    free_traj_set(&set_glue);

	} // end set loop

  printf("\nAnalyzed %d sets of half-trajectories, %ld full trajectories in total\n", nset, traj_tot);

  double res[2];
  jackknife_double(dphidt, nset, blocks, res);

  printf("\n --Heatbath <|dphi/dt|> = %g +/- %g\n", res[0], res[1]);
  for (int i=0; i<nset; i++) {
    if (res[1] <= dphidt_err[i])
      printf("!!! Warning: final error smaller than error in set %d (+/- %g)\n", i, dphidt_err[i]);
  }

  jackknife_double(d, nset, blocks, res);

  printf("\n --Dynamical prefactor: d = %g +/- %g\n", res[0], res[1]);
  for (int i=0; i<nset; i++) {
    if (res[1] <= derr[i])
      printf("!!! Warning: final error smaller than error in set %d (+/- %g)\n", i, derr[i]);
  }

  jackknife_double(prod, nset, blocks, res);

  printf("\n -- d * |<dphi/dt>| = %g +/- %g\n", res[0], res[1]);
  for (int i=0; i<nset; i++) {
    if (res[1] <= prod_err[i])
      printf("!!! Warning: final error smaller than error in set %d (+/- %g)\n", i, prod_err[i]);
  }

}

/* jackknife error estimate for array average of doubles.
*Stores in 'res' as res[0] = avg, res[1] = error  */
void jackknife_double(double* arr, long datalen, int blocks, double* res) {

  if (datalen <= 0) {
    printf("!!! Error in jackknife_double: no data??\n");
    die(-21);
  }

  double avg = 0.0;
  for (long k=0; k<datalen; k++) avg += arr[k];
  avg /= (double) datalen;
  res[0] = avg;

  double err = 0.0;
  if (blocks >= 1 && blocks <= datalen)  {

    int blen = datalen / blocks;
    int blen_last = blen;
    // if the data doesn't split evenly, last block may have different length from the others
    if (datalen % blocks != 0) {
      blen_last = blen + datalen % blocks;
    }

    for (int b_id=0; b_id<blocks; b_id++) {
      long jack_len;
      if (b_id == blocks-1) jack_len = (blocks-1)*blen;
      else jack_len = (blocks-2)*blen + blen_last;

      double data_jack[jack_len];
      // take out block 'b_id'
      int cp_tot = 0;
      int skip;
      for (int k=0; k<datalen; k++) {
        int end = blen*(b_id+1)-1;
        if (b_id == blocks-1) end = datalen; // special treatment for last block

        if (k >= blen*b_id && k <= end) skip = 1;
        else skip = 0;

        if (!skip) {
          data_jack[cp_tot] = arr[k];
          cp_tot++;
        }
      }
      if (cp_tot != jack_len) {
        printf("!!! Warning: error making jackknife blocks (jack_len=%ld, but collected %d)!\n", jack_len, cp_tot);
      }

      // copy done, calc jackknife average
      double jack_avg = 0.0;
      for (long l=0; l<jack_len; l++) jack_avg += data_jack[l];
      jack_avg /= (double) jack_len;

      err += (jack_avg - avg) * (jack_avg - avg);

    } // end block loop

  } else {
    printf("!!! Warning: no jackknife (doubles) (requested %d blocks but data length is %ld)\n", blocks, datalen);
    res[1] = 0.0;
    return;
  }

  err *= blocks / ((double) (blocks-1));

  res[1] = sqrt(err);
}


/* Calculate real time stuff: dphi/dt at crossing of critical phi_c, and
* the dynamical prefactor d. Stores res[0] = dphi/dt, phi[1] = dyn. prefactor.
* If a trajectory does not cross phi_c, it still contributes to the statistics of d
* but not for dphi/dt.
* For d they need to be included because it tells you the fraction of configurations
* in the region |phi - phi_c| < eps that actually tunnel */
void dynamics(traj_set* set, double min, double max, double phi_c, double* res) {

  double d = 0.0;
  double dphidt = 0.0;
  traj* t;
  int cross_tot = 0;
  res[0] = 0.0; res[1] = 0.0;

  for (int i=0; i<set->ntraj; i++) {
    t = &set->traj[i];

    double start = t->val[0];
    double end = t->val[t->len-1];

    int tunnel = 1;
    // no tunneling if the trajectory ends on the same side as it started from:
    if (start <= min && end <= min) tunnel = 0;
    else if (start >= max && end >= max) tunnel = 0;

    double der = 0.0;
    int crossings = 0;
    // find crossings of the critical value phi_c and calculate derivative there
    for (int n=1; n<t->len; n++) {
      if ( (t->val[n] - phi_c) * (t->val[n-1] - phi_c) < 0 ) {
        crossings++;
        der += fabs( (t->val[n] - t->val[n-1]) / (t->t[n] - t->t[n-1]) );
      }
    }

    cross_tot += crossings;
    dphidt += der;

    if (crossings > 0) d += tunnel / ((double) crossings);

  } // end traj loop

  // now average. d is averaged over all trajectories and dphi/dt over the number of crossings
  if (cross_tot > 0) {
    res[0] = dphidt / ((double) cross_tot);
    res[1] = d / ((double) set->ntraj);
  }

}


/* calculate the "dynamical prefactor" d from a trajectory set */
double dynamical_prefactor(traj_set const* set, double min, double max, double phi_c) {

  double d = 0.0;
  traj* t;
  int ntraj = 0;

  for (int i=0; i<set->ntraj; i++) {
    t = &set->traj[i];

    double start = t->val[0];
    double end = t->val[t->len-1];

    int tunnel = 1;
    // no tunneling if the trajectory ends on the same side as it started from:
    if (start <= min && end <= min) tunnel = 0;
    else if (start >= max && end >= max) tunnel = 0;

    ntraj++;

    if (!tunnel) continue;

    // calculate crossings of the critical bubble phi_c
    int last_max;
    if (start >= max) last_max = 1;
    else if (start <= min) last_max = 0;
    else {
      printf("!!! Warning: unfinished trajectory, skipping...\n");
      ntraj--;
      continue;
    }

    int crossings = 0;
    for (int n=0; n<t->len; n++) {
      if (last_max && t->val[n] < phi_c) {
        crossings++;
        last_max = 0;
      } else if (!last_max && t->val[n] > phi_c) {
        crossings++;
        last_max = 1;
      }
    }

    if (crossings > 0.0) {
      d += 1.0 / ((double) crossings);
    } else {
      // no crossings, so something went wrong
      printf("!!! Warning: no crossing of phi_c in traj %d !!\n", i);
      ntraj--;
    }
  }

  return d / ((double) ntraj);
}

void read_traj_set(char* fname, traj_set* set, int id) {

	FILE *f = fopen(fname, "r");

  char line[255];
	char header[255];
	char nextheader[255];
	snprintf(header, 255, "Begin set %d, start", id);
	snprintf(nextheader, 255, "Begin set %d, start", id+1);

	int ntraj = 0;
	int ok = 0;
	int init_max = 500;

	while (fgets(line, sizeof(line), f)) {
		if (strstr(line, header) != NULL) {
			// found set 'id'
			ok = 1;
			continue;
		} else if (strstr(line, nextheader) != NULL) {
			// reached next set of trajectories, exit
			ok = 0;
			break;
		} else if (ok && strstr(line, "--- Trajectory ") != NULL) {

			// begin and read new trajectory
			make_traj(&set->traj[ntraj], init_max);
			read_trajectory(f, &set->traj[ntraj]);
			ntraj++;
		}
	}

	fclose(f);
	realloc_traj_set(set, ntraj);
}

void read_trajectory(FILE* f, traj* t) {
	fpos_t fpos; // position in the file
	int read = 2;
	int n = 0;
	double time_read, val_read;

	// read the trajectory, will stop once fscanf does not found 2 numbers
	while (1==1) {

		if (fgetpos(f, &fpos) != 0) printf("Warning: fpos error!!\n");

		read = fscanf(f, "%lf %lf\n", &time_read, &val_read);
		if (read != 2) {
			// rewind to the line that was just read
			fsetpos(f, &fpos);
			break;
		}

		if (n >= t->len) {
			// need to realloc
			realloc_traj(t, t->len + 100);
		}
		t->t[n] = time_read;
		t->val[n] = val_read;

		// printf("%g %g\n", t->t[n], t->val[n]); // print trajectory, for debugging
		n++;
	}
	// read done, now just realloc
	realloc_traj(t, n);
}


void make_traj(traj* t, int len) {
	t->t = calloc(len, sizeof(*t->t));
	t->val = calloc(len, sizeof(*t->val));
	t->len = len;
}

void realloc_traj(traj* t, int newlen) {
	if (t->len == newlen) return;
	t->len = newlen;
	t->t = realloc(t->t, newlen * sizeof(*t->t));
	t->val = realloc(t->val, newlen * sizeof(*t->val));

	if (t->t == NULL || t->val == NULL) {
		printf("!!! realloc error!! in readtraj.c\n");
		die(-21);
	}
}

void free_traj(traj* t) {
	free(t->t);
	free(t->val);
};


void make_traj_set(traj_set* set, int ntrajs) {
	set->ntraj = ntrajs;
	set->traj = malloc(ntrajs * sizeof(*set->traj));
}


void realloc_traj_set(traj_set* set, int newlen) {
	if (newlen == set->ntraj) return;
	set->traj = realloc(set->traj, newlen * sizeof(*set->traj));
	set->ntraj = newlen;
}

void free_traj_set(traj_set* set) {
	for (int i=0; i<set->ntraj; i++) free_traj(&set->traj[i]);
  free(set->traj);
}

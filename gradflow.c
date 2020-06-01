
/** @file gradflow.c
*
* Routines for performing gradient (Wilson) flow on the fields.
* References: 0907.5491, 1006.4518.
*
* Measurements are are stored in "measure_flow". First column is the (dimensionless)
* time and the other columns follow the same pattern as in "labels".
*
* Plan is to first calculate the gradient force at time t on all fields everywhere,
* then update all fields to time t + dt. This way the ordering of updates does not matter.
*
* For gauge links, we flow is
*   (d/dt) V_t(x,mu) = -g^2 [\partial_{x,mu} S(V_t)] V_t(x,mu)
* with V_t(x,mu) = U_mu(x) at t=0. Here S(V_t) is the action due to the gauge
* link V_t(x,mu), and the derivative can be calculated as a projection as in
* 0907.5491. Here I use the "Euler" scheme for integrating the flow, so the derivative
* is assumed to be a constant in the interval [t, t+dt]
*
* TODO: Add Higgs!!
*
*/

#ifdef GRADFLOW // do nothing if compiler flag is not set

#include "su2.h"


/* Do a "real time" gradient flow of the fields and measure stuff
* as a function of the time.
* flow_id is an identifier for the current flow and is used as a "header"
* in the measurement file.
* Requires a "weight" struct to be compatible with measure().
* This does not modify the original field config.
*/
void grad_flow(params const* p, fields const* f, comlist_struct* comlist,
                  weight* w, double t_max, double dt, int flow_id) {

  fields flow; // flowing fields
  fields forces; // gradient forces for each field

  alloc_fields(p, &flow);
  alloc_fields(p, &forces);

  // copy starting configuration to "flow"
  copy_fields(p, f, &flow);
  sync_halos(&flow, p, comlist); // initialize halos in "flow"


  FILE* file;
  if (!p->rank) {
    file = fopen("measure_flow", "a");
    // write header for the current set of measurements
    fprintf(file, "\n =========== Flow id: %d ===========\n", flow_id);
  }

  double t = 0.0;
  // initial measurements at flow time t = 0
  if (!p->rank) {
    fprintf(file, "%g ", t); // first column is time, rest come from measure()
  }
  measure(file, &flow, p, w);

  /* "time" loop */
  while (t < t_max) {

    // if time after updating surpasses t_max, calculate only up to t_max
    if (t + dt > t_max) {
      dt = t_max - t;
    }

    calc_gradient(p, &flow, &forces);

    /* update everything. Halos can be synced afterwards,
     * because the force is known already. */
    flow_gauge(p, &flow, &forces, dt); // SU(2) gauge links
    #ifdef HIGGS

    #endif

    #ifdef TRIPLET

    #endif


    sync_halos(&flow, p, comlist);

    t += dt;

    if (!p->rank) {
      fprintf(file, "%g ", t);
    }
    measure(file, &flow, p, w);

  }

  if (!p->rank) {
    fclose(file);
  }

  free_fields(p, &forces);
  free_fields(p, &flow);

}


/* Calculate gradient force for the SU(2) link at a given site and direction
* and store in "force". "f" contains the field configuration at current flow time.
*
* If the flow equation is:    d/dt V_\mu(x,t) = Z[V] V_\mu(x,t),
* then the "force" here refers to the Lie-algebra values coefficient Z
* which I write, in terms of Pauli matrices sigma^a, as
*   Z[V] = i F^a sigma^a .
*
* So this routine calculates the components F^a.
* NOTE: the real force is an adjoint vector with 3 components, but here I take
* "force" to have 4 components with the 0. component not being used for anything.
* This is convenient because then I can make a "fields" struct and store the link
* force in force->su2link[x][mu]
*
*/
void grad_force_link(params const* p, fields const* f, double* force, long i, int dir) {

  /* Pure gauge part: the Wilson force is
  *     Z = -0.5 [U_mu(x) S_mu(x) - U^+_mu(x) S^+_mu(x)
  *          - 1/N * Tr(U_mu(x) S_mu(x) - U^+_mu(x) S^+_mu(x) )]
  * where S_mu(x) is the staple calculated by su2staple_wilson():
  *
  * S_mu(x) = \sum_{nu != mu} (U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
  *																	+ U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu) )
  */

  double s[SU2LINK];
  su2staple_wilson(f, p, i, dir, s);

  /* staple s is a non-unitary matrix that can still be parametrized as
    s = I s[0] + i s[a]*sigma_a , with real s[a], simply because it is a sum of such matrices.
  */

  #ifdef TRIPLET
    /* Add triplet contribution to the "staple", which actually depends on U^+_mu(x):
      s <- s + g^2 Sigma(x+mu) U^+_mu(x) Sigma(x)
    */
    double* u = f->su2link[i][dir];
    double* a1 = f->su2triplet[i];
    long next = p->next[i][dir];
    double* a2 = f->su2triplet[next];

    s[0] += (1.0/p->betasu2) * (a1[0]*a2[0]*u[0] + a1[1]*a2[1]*u[0] + a1[2]*a2[2]*u[0]
          - a1[2]*a2[1]*u[1] + a1[1]*a2[2]*u[1] + a1[2]*a2[0]*u[2] - a1[0]*a2[2]*u[2]
          - a1[1]*a2[0]*u[3] + a1[0]*a2[1]*u[3])/4.0;

    s[1] += (1.0/p->betasu2) * (a1[2]*a2[1]*u[0] - a1[1]*a2[2]*u[0] - a1[0]*a2[0]*u[1]
          + a1[1]*a2[1]*u[1] + a1[2]*a2[2]*u[1] - a1[1]*a2[0]*u[2] - a1[0]*a2[1]*u[2]
          - a1[2]*a2[0]*u[3] - a1[0]*a2[2]*u[3])/4.0;

    s[2] += (1.0/p->betasu2) * (-(a1[2]*a2[0]*u[0]) + a1[0]*a2[2]*u[0] - a1[1]*a2[0]*u[1]
          - a1[0]*a2[1]*u[1] + a1[0]*a2[0]*u[2] - a1[1]*a2[1]*u[2] + a1[2]*a2[2]*u[2]
          - a1[2]*a2[1]*u[3] - a1[1]*a2[2]*u[3])/4.0;

    s[3] += (1.0/p->betasu2) * (a1[1]*a2[0]*u[0] - a1[0]*a2[1]*u[0] - a1[2]*a2[0]*u[1]
          - a1[0]*a2[2]*u[1] - a1[2]*a2[1]*u[2] - a1[1]*a2[2]*u[2] + a1[0]*a2[0]*u[3]
          + a1[1]*a2[1]*u[3] - a1[2]*a2[2]*u[3])/4.0;

  #endif

  // staple done, multiply by the link
  double m[SU2LINK];
  memcpy(m, f->su2link[i][dir], SU2LINK*sizeof(m[0]));
  su2rot(m, s); // m <- U_mu(x) S_mu(x)

  /* project onto su(2) algebra:
  *   Z = -0.5 * [m - m^+ - 0.5 * Tr(m - m^+) ] , which is of the form
  *   Z = i z^a sigma^a,   with real z^a.
  * For SU(2) the trace part vanishes, and m - m^+ = 2 i m[a] sigma[a]
  * because all our m[i] are real.
  */
  // "force" is an adjoint vector while force has 4 components,
  // but I can just set force[0] to vanish.
  force[0] = 0.0;
  for(int a=1; a<SU2TRIP+1; a++) {
    force[a] = -1.0 * m[a];
  }
  // link force done

}

/* Calculate gradient force for each field at all sites (not halos)
* and store in the "fields" array "forces".
*/
void calc_gradient(params const* p, fields const* f, fields* forces) {

  for (long i=0; i<p->sites; i++) {

    for (int dir=0; dir<p->dim; dir++) {
      /* store force as a "link" matrix. In reality it is an adjoint vector,
      * so the 0. component is not used */
      grad_force_link(p, f, forces->su2link[i][dir], i, dir);
    }

  #ifdef HIGGS

  #endif

  #ifdef TRIPLET

  #endif
  } // end site loop

}


/* Flow the SU(2) gauge field one step forward in time.
* force[i][mu] should give the gradient force on the link U_mu at site i.
* Here again force[i][mu] is assumed to have 4 components, with the 0. component not used.
*/
void flow_gauge(params const* p, fields* flow, fields const* forces, double dt) {


  for (long i=0; i<p->sites; i++) {
    for (int dir=0; dir<p->dim; dir++) {
      /* d/dt V = (i z_a sigma^a) V
      * => approximate V(t+dt) = exp[i dt*z_a sigma^a ] V(t)
      * and z_a are in force[i][mu]
      */
      double z[SU2TRIP];
      double mod = 0.0;
      for (int a=0; a<SU2TRIP; a++) {
        z[a] = dt * forces->su2link[i][dir][a+1]; // offset force by 1 since 0-component is placeholder
        mod += z[a]*z[a];
      }

      if (mod <= 0.0) {
        /* force is zero, so the update does nothing.
        This happens if the matrix U_mu(x).S_mu(x) is Hermitian, or if dt=0 */
        return;
      }

      mod = sqrt(mod); // length of the force vector

      /* exp[i z_a sigma^a] = exp[i mod n_a sigma^a ]
      *  with unit vector n_a = z_a / mod */
      double u[SU2LINK];
      u[0] = cos(mod);
      for (int a=0; a<SU2TRIP; a++) {
        u[a+1] = sin(mod) * z[a] / mod;
      }
      // now u = exp[i dt*z_a sigma^a]
      su2rot(u, flow->su2link[i][dir]); // u <- exp[...].U_mu(x)

      memcpy(flow->su2link[i][dir], u, SU2LINK*sizeof(u[0])); // update the flow link


      /*
      // test:
      double det = su2sqr(flow->su2link[i][dir]);
      if (det - 1.0 > 10e-6) {
        printf("Error in flow_gauge()! Got det = %lf\n", det);
      }
      */

    } // end dir

  } // end i
  // gauge update done
}



#endif // end #ifdef GRADFLOW

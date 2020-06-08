#ifndef SU2_H
#define SU2_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>


#include "comms.h"

// degrees of freedom per site for different fields
#define SU2DB 4
#define SU2LINK 4
#define SU2TRIP 3

// update algorithms
#define METROPOLIS 1 // Metropolis
#define HEATBATH 2 // Heatbath
#define OVERRELAX 3 // Overrelaxation

// parity identifiers
#define EVEN 0
#define ODD 1

// multicanonical order parameters
#define PHISQ 1
#define SIGMASQ 2
#define PHI2MINUSSIGMA2 3
#define PHI2SIGMA2 4

// nasty globals for keeping track of evaluation time
double waittime;
double Global_comms_time, Global_total_time;
double Global_current_action; // for debugging gradient flows etc

// general multipurpose struct params.
typedef struct {
	// layout parameters for MPI
	int rank, size;

	// lattice dimensions, set in get_parameters()
	int dim;
	int *L;
	long vol;
	// slicing the lattice for MPI, set in layout.c
	int *nslices; // how many slices in each direction
	long *sliceL; // how many sites per slice in each direction
	long sites; // how many sites in total in one hypercubic slice
	long halos; // how many artificial sites are needed for haloing
	// how many actual sites we have. This can be less than sites + halos,
	// because we remove halos that correspond to real sites in my the same node.
	long sites_total;
	// neighboring sites: next[i][dir] gives the index of the next site after i in direction dir
	long **next;
	long **prev;
	/* coords[i][dir] = x_dir coordinate on the full lattice of site i.
	* Typically not used after layout(), so could free to save memory... */
	long **coords;
	long *offset; // coordinate offset in my node wrt. the full lattice
	long firstsite; // index of (x,y,z)=(0,0,0), with (x,y,z) being coordinates on the node

	// parity of a site is EVEN if the physical index x+y+z+... is even, ODD otherwise
	char *parity;
	long evensites, oddsites;
	long evenhalos, oddhalos;
	// in layout.c we reorder lattice sites so that EVEN sites come first.
	int reorder_parity; // for debugging purposes

	#ifdef MEASURE_Z
		/* Single out a particular direction ("z coordinate") to study e.g. wall profile, correlation lengths.
		* These are initialized in init_z_coord() which is called by layout() */
		int z_dir; // specify the z direction, default = longest direction
		long sites_per_z, offset_z; // offset = physical value of z at the "first" site in my node
		long area; // how many sites per z coordinate on the whole lattice; area = vol / L[z_dir]
		// list of all real sites at fixed z coord in no particular order
		long** site_at_z; // sites_at_z[z][j], 0<=j<sites_per_z
		int n_meas_z; // how many quantities to measure along z
		int meas_interval_z; // read in from config file in get_parameters()
	#endif

	#ifdef GRADFLOW
		// these are all read from config (in parameters.c)
		int do_flow; // if 0, will not do gradient flows
		double flow_dt; // timestep
		double flow_t_max; // max time for a single flow (start from t=0)
		int flow_interval; // how many non-flow iterations between flows
		int flow_meas_interval; // how often to measure during flowing (in units of dt)
	#endif

	// max iterations etc
	int reset;
	long iterations;
	long checkpoint;
	long n_thermalize;
	long interval;
	FILE *resultsfile;
	char latticefile[100];
	int run_checks;
	int do_local_meas;

	int multicanonical;
	/* randomize the ordering of updates in update_lattice() or not.
	* Typically 0 (non-random), but for heatbath trajectories this is set to 1.
	*/
	int random_sweeps;

	// Parameters in the action
	// possible generalization: make separate structures
	// for each field in the theory that contain their respective params
	double betasu2; // 4/(g^2)
	double betau1; // 2/(g'^2)
	// doublet
	double msq_phi; // Higgs mass
	double lambda_phi; // Higgs quartic
	// triplet
	double msq_triplet;
	double a2; // Higgs portal
	double b4; // self quartic

	// initial values for fields
	double phi0; // doublet
	double sigma0; // triplet

	// How to update the fields
	short algorithm_su2link;
	short algorithm_u1link;
	short algorithm_su2doublet;
	short algorithm_su2triplet;

	// How many times to update a field per sweep
	short update_links;
	short scalar_sweeps; // update all scalars n times per iteration
	// additional sweeps on top of scalar_sweeps
	short update_su2doublet;
	short update_su2triplet;

} params;


// Field structure. Note that this can safely be passed by value to updating functions
// after the contents have been alloc'd
typedef struct {

	//double *su2singlet;
	double ***su2link;
	double **su2doublet;
	double **su2triplet;
	double **u1link;

	// backup arrays. these are used in global multicanonical steps in case the update sweep needs to be undone
	double **backup_doublet;
	double **backup_triplet;
} fields;

typedef struct {

	long iter;

	// keep track of when a metropolis sweep should be forced
	int higgs_sweeps, triplet_sweeps;

	// count metropolis updates
	long total_su2link, accepted_su2link;
	long total_u1link, accepted_u1link;
	long total_doublet, accepted_doublet;
	long total_triplet, accepted_triplet;
	// count overrelax updates
	long total_overrelax_doublet, acc_overrelax_doublet;
	long total_overrelax_triplet, acc_overrelax_triplet;

	long total_muca, accepted_muca;
} counters;


// multicanonical weight
typedef struct {
	int orderparam; // specify which order parameter to use

	/* current value of the order parameter, separated by parity.
	* e.g. param_value[0] is the full contribution from EVEN sites
	* and param_value[1] is the contribution from ODD sites */
	double param_value[2];

	/* 1 if muca accept/reject is to be performed after update sweeps, 0 otherwise.
	* This is used by e.g. realtime trajectory routines to temporarily disable weighting.
	* Will not update weight either if set to 0. Read in update sweep routines.
	*/
	int do_acceptance;

	long bins;
	double min, max; // range of orderparam values where weight is to be updated
	double min_abs, max_abs; // weighting range. can be larger than the range where weight is actually updated
	long min_bin, max_bin; // indices of the bins containing w.min and w.max
	double* pos; // weight "position", i.e. values of the order param in the given range
	double* W; // value of the weight function at each position
	double increment; // how much the weight is increased after visiting a bin
	double reduction_factor; // after probing the full order param range, w.increment is multiplied by this number

	int* hits; // keep track of which bins we have visited
	int update_interval; // how many measurements until weight is updated
	int m; // current number of measurements (resets after weight update)

 	int last_max; // 1 if system recently visited the bin containing w.max (keep track of tunneling)
	char readonly; // 1 if weight is to be updated, 0 otherwise
	char weightfile[100]; // file name
} weight;


/* inlines */
inline char otherparity(char parity) {
	if (parity == EVEN)
		return ODD;
	else
		return EVEN;
}

// comms.c (move elsewhere later?)
void make_comlists(params *p, comlist_struct *comlist);
void reorder_sitelist(params* p, sendrecv_struct* sr);
void reorder_comlist(params* p, comlist_struct* comlist);
double reduce_sum(double res);
double allreduce(double res);
void bcast_int(int *res);
void bcast_long (long *res);
void bcast_double(double *res);
void bcast_int_array(int *arr, int size);
void barrier();
// gauge links:
void update_gaugehalo(comlist_struct* comlist, char parity, double*** field, int dofs, int dir);
#ifdef MPI
void send_gaugefield(sendrecv_struct* send, MPI_Request* req, char parity, double*** field, int dofs, int dir);
void recv_gaugefield(sendrecv_struct* recv, char parity, double*** field, int dofs, int dir);
#endif
// non-gauge fields:
void update_halo(comlist_struct* comlist, char parity, double** field, int dofs);
#ifdef MPI
void recv_field(sendrecv_struct* recv, char parity, double** field, int dofs);
void send_field(sendrecv_struct* send, MPI_Request* req, char parity, double** field, int dofs);
#endif
void test_comms(params p, comlist_struct comlist);
void test_comms_individual(params p, comlist_struct comlist);

// layout.c
void layout(params *p, comlist_struct *comlist);
void make_slices(params *p);
void sitemap(params *p);
void set_parity(params *p);
void paritymap(params* p, long* newindex);
void remap_latticetable(params* p, long** arr, long* newindex, long maxindex);
void remap_neighbor_table(params* p, long** arr, long* newindex, long maxindex);
void remap_lattice_arrays(params* p, long* newindex, long maxindex);
long findsite(params* p, long* x, int include_halos);
void test_coords(params p);
void test_neighbors(params p);
void indexToCoords(short dim, int* L, long i, long* x);
long coordsToIndex(short dim, int* L, long* x);
int coordsToRank(params p, long* x);
void printf0(params p, char *msg, ...);
void die(int howbad);
void print_lattice_2D(params *p);

// alloc.c
double *make_singletfield(long sites);
double **make_field(long sites, int dofs);
double ***make_gaugefield(long sites, int dim, int dofs);
void free_singletfield(double *field);
void free_field(double **field);
void free_gaugefield(long sites, double ***field);
void alloc_fields(params const* p, fields *f);
void free_fields(params const* p, fields *f);
void alloc_lattice_arrays(params *p, long sites);
long **alloc_latticetable(int dim, long sites);
long **realloc_latticetable(long** arr, int dim, long oldsites, long newsites);
void realloc_lattice_arrays(params *p, long oldsites, long newsites);
void free_latticetable(long** list);
void free_lattice_arrays(params *p);
void free_comlist(comlist_struct* comlist);

// su2u1.c
double su2sqr(double *u);
void su2rot(double *u1, double *u2);
double su2ptrace(fields const* f, params const* p, long i, int dir1, int dir2);
long double local_su2wilson(fields const* f, params const* p, long i);
double localact_su2link(fields const* f, params const* p, long i, int dir);
void su2staple_wilson(fields const* f, params const* p, long i, int dir, double* V);
void su2link_staple(fields const* f, params const* p, long i, int dir, double* V);
double su2trace4(double *u1, double *u2, double *u3, double *u4);
void su2staple_counterwise(double* V, double* u1, double* u2, double* u3);
void su2staple_clockwise(double* V, double* u1, double* u2, double* u3);
double hopping_trace(double* phi1, double* u, double* phi2);
double hopping_trace_su2u1(double* phi1, double* u, double* phi2, double a);
double hopping_trace_triplet(double* a1, double* u, double* a2);
// U(1) routines
double u1ptrace(fields const* f, params const* p, long i, int dir1, int dir2);
long double local_u1wilson(fields const* f, params const* p, long i);
double localact_u1link(fields const* f, params const* p, long i, int dir);
// doublet routines
double doubletsq(double* a);
long double avg_doublet2(fields const* f, params const* p);
long double avg_doublet4(fields const* f, params const* p);
double hopping_doublet_forward(fields const* f, params const* p, long i, int dir);
double hopping_doublet_backward(fields const* f, params const* p, long i, int dir);
double covariant_doublet(fields const* f, params const* p, long i);
double localact_doublet(fields const* f, params const* p, long i);
double higgspotential(fields const* f, params const* p, long i);
// triplet routines
double tripletsq(double* a);
double hopping_triplet_forward(fields const* f, params const* p, long i, int dir);
double hopping_triplet_backward(fields const* f, params const* p, long i, int dir);
double covariant_triplet(fields const* f, params const* p, long i);
double localact_triplet(fields const* f, params const* p, long i);

// metropolis.c
int metro_su2link(fields* f, params const* p, long i, int dir);
int metro_u1link(fields* f, params const* p, long i, int dir);
int metro_doublet(fields* f, params const* p, long i);
int metro_triplet(fields* f, params const* p, long i);

// heatbath.c
int heatbath_su2link(fields* f, params const* p, long i, int dir);

// overrelax.c
double polysolve3(long double a, long double b, long double c, long double d);
int overrelax_doublet(fields* f, params const* p, long i);
int overrelax_triplet(fields* f, params const* p, long i);

// update.c
void update_lattice(fields* f, params const* p, comlist_struct* comlist, counters* c, weight* w);
void checkerboard_sweep_su2link(fields* f, params const* p, counters* c, int parity, int dir);
void checkerboard_sweep_u1link(fields* f, params const* p, counters* c, int parity, int dir);
int checkerboard_sweep_su2doublet(fields* f, params const* p, counters* c, weight* w, int parity, int metro);
int checkerboard_sweep_su2triplet(fields* f, params const* p, counters* c, weight* w, int parity, int metro);
void sync_halos(fields* f, params const* p, comlist_struct* comlist);
void shuffle(int *arr, int len);

// init.c
void setsu2(fields f, params p);
void random_su2link(double *su2);
void setu1(fields f, params p);
void setfields(fields f, params p);
void setdoublets(fields f, params p);
void settriplets(fields f, params p);
void copy_fields(params const* p, fields const* f_old, fields* f_new);
void init_counters(counters* c);

// checkpoint.c
void print_acceptance(params p, counters c);
void read_field(params p, FILE *file, double *field, int size);
void write_field(params p, FILE *file, double *field, int size);
void save_lattice(params p, fields f, counters c);
void load_lattice(params const* p, fields* f, counters* c, comlist_struct* comlist);

// parameters.c
void check_set(int set, char *name); // helper routine
void get_parameters(char *filename, params *p);
void get_weight_parameters(char *filename, params *p, weight* w);
void print_parameters(params p);
void read_updated_parameters(char *filename, params *p);


// measure.c
void measure(FILE* file, fields const* f, params const* p, weight* w);
// weight not necessarily constant because measure() can recalculate the order parameter
double action_local(fields const* f, params const* p, long i);
void print_labels();
void measure_local(char* fname, params const* p, fields const* f);
void print_labels_local(params* p, char* fname);

// multicanonical.c
void load_weight(params const* p, weight *w);
void save_weight(params const* p, weight const* w);
double get_weight(weight const* w, double val);
void update_weight(params const* p, weight* w);
int multicanonical_acceptance(params const* p, weight* w, double oldval, double newval);
long whichbin(weight const* w, double val);
double calc_orderparam(params const* p, fields const* f, weight* w, char par);
void check_tunnel(params const* p, weight *w);
void store_muca_fields(params const* p, fields* f, weight* w);
void reset_muca_fields(params const* p, fields* f, weight* w, char par);
void alloc_backup_arrays(params const* p, fields* f, weight const* w);
void free_muca_arrays(fields* f, weight *w);

#ifdef MEASURE_Z
	// z_coord.c
	void init_z_coord(params* p);
	void print_z_labels(params* p);
	void measure_along_z(fields const* f, params const* p, long id);
	void prepare_wall(fields* f, params const* p, comlist_struct* comlist);
#endif


#ifdef TRIPLET
	// magfield.c
	void matmat(double *in1, double *in2, int dag);
	void projector(double *proj, double *adjoint);
	void project_u1(params const* p, fields const* f, long i, int dir, double* pro);
	double alpha_proj(params const* p, fields const* f, long i, int dir1, int dir2);
	double magfield(params const* p, fields const* f, long i, int dir);
	double magcharge_cube(params const* p, fields const* f, long i);
#endif


#ifdef GRADFLOW
	// gradflow.c
	void grad_flow(params const* p, fields const* f, comlist_struct* comlist, weight* w, double t_max, double dt, int flow_id);
	void grad_force_link(params const* p, fields const* f, double* force, long i, int dir);
	void grad_force_triplet(params const* p, fields const* f, double* force, long i);
	void calc_gradient(params const* p, fields const* f, fields* forces);
	void flow_gauge(params const* p, fields* flow, fields const* forces, double dt);
	void flow_triplet(params const* p, fields* flow, fields const* forces, double dt);
#endif


#endif // end #ifndef SU2_H

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


typedef unsigned int uint;
typedef unsigned short ushort;

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

// nasty global...
double waittime;

// general multipurpose struct params.
// Multicanonical can modify the lattice parameters (redefine to include weight),
// so this should be passed by reference. TODO make a separate lattice structure
// that contains the layouting
typedef struct {
	// layout parameters for MPI
	int rank, size;

	// lattice dimensions
	int dim;
	int *L;
	long vol;
	// slicing the lattice for MPI
	int *nslices; // how many slices in each direction
	long *sliceL; // how many sites per slice in each direction
	long sites; // how many sites in total in one hypercubic slice
	long halos; // how many artificial sites are needed for haloing
	// how many actual sites we have. This can be less than sites + halos,
	// because we remove halos that correspond to real sites in my the same node.
	long sites_total;
	// neighboring sites
	//next[i][dir] gives the index of the next site after i in direction dir
	long **next;
	long **prev;

	// parity of a site is EVEN if the physical index x+y+z+... is even, ODD otherwise
	char *parity;
	long evensites, oddsites;
	long evenhalos, oddhalos;
	// in layout.c we reorder lattice sites so that EVEN sites come first.
	int reorder_parity; // for debugging purposes

	// max iterations etc
	int reset;
	long iterations;
	long checkpoint;
	long interval;
	FILE *resultsfile;
	char latticefile[100];
	int run_checks;

	int multicanonical;

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

	// keep track of evaluation time
	double comms_time;
	double total_time;
	long iter;

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
	int orderparam;
	/* current value of the order parameter, separated by parity.
	* e.g. param_value[0] is the full contribution from EVEN sites
	* and param_value[1] is the contribution from ODD sites */
	double param_value[2];

	long bins;
	double min, max;
	double dbin; // size (width) of one bin
	double* pos; // weight "position", i.e. values of the order param in the given range
	double* W; // value of the weight function at each position
	double increment; // how much the weight is increased after visiting a bin
	double outsideW_min, outsideW_max; // weight values outside binning range

	int* hits; // keep track of which bins we have visited
	int update_interval; // how many measurements until weight is updated
	int m; // current number of measurements (resets after weight update)

 	int last_max; // 1 if system was recently in one of the last bins (keep track of tunneling)
	char readonly;
	char absolute_bounds; // if 1, uses infinite weight outside binning range
	char weightfile[100];
} weight;


/* inlines */
inline char otherparity(char parity) {
	if (parity == EVEN)
		return ODD;
	else
		return EVEN;
}

// comms.c (move elsewhere later?)
void make_comlists(params *p, comlist_struct *comlist, long** xphys);
void reorder_sitelist(params* p, sendrecv_struct* sr);
void reorder_comlist(params* p, comlist_struct* comlist);
double reduce_sum(double res);
double allreduce(double res);
void bcast_int(int *res);
void bcast_long (long *res);
void bcast_double(double *res);
void barrier();
// gauge links:
double update_gaugehalo(comlist_struct* comlist, char parity, double*** field, int dofs, int dir);
#ifdef MPI
void send_gaugefield(sendrecv_struct* send, MPI_Request* req, char parity, double*** field, int dofs, int dir);
void recv_gaugefield(sendrecv_struct* recv, char parity, double*** field, int dofs, int dir);
#endif
// non-gauge fields:
double update_halo(comlist_struct* comlist, char parity, double** field, int dofs);
#ifdef MPI
void recv_field(sendrecv_struct* recv, char parity, double** field, int dofs);
void send_field(sendrecv_struct* send, MPI_Request* req, char parity, double** field, int dofs);
#endif
void test_comms(params p, comlist_struct comlist);
void test_comms_individual(params p, comlist_struct comlist, long** xphys);

// layout.c
void layout(params *p, comlist_struct *comlist);
void make_slices(params *p);
void sitemap(params *p, long** xphys, long* newindex);
void calculate_neighbors(params *p, long** xphys, long* newindex);
void set_parity(params *p, long** xphys);
void paritymap(params* p, long* newindex);
void remap_lattice_arrays(params* p, long** xphys, long* newindex);
long findsite(params p, long* x, long** xphys, int include_halos);
void test_xphys(params p, long** xphys);
void test_neighbors(params p, long** xphys);
void indexToCoords(ushort dim, uint* L, long i, long* x);
long coordsToIndex(ushort dim, uint* L, long* x);
int coordsToRank(params p, long* xphys);
void printf0(params p, char *msg, ...);
void die(int howbad);
void print_lattice_2D(params *p, long** xphys);

// alloc.c
double *make_singletfield(long sites);
double **make_field(long sites, int dofs);
double ***make_gaugefield(long sites, int dim, int dofs);
void free_singletfield(double *field);
void free_field(double **field);
void free_gaugefield(long sites, double ***field);
void alloc_fields(params p, fields *f);
void free_fields(params p, fields *f);
void alloc_lattice_arrays(params *p);
long **alloc_latticetable(int dim, long sites);
void free_latticetable(long** list);
void free_lattice_arrays(params *p);
void free_comlist(comlist_struct* comlist);

// su2u1.c
double su2sqr(double *u);
void su2rot(double *u1, double *u2);
double su2ptrace(fields f, params p, long i, int dir1, int dir2);
long double local_su2wilson(fields f, params p, long i);
double localact_su2link(fields f, params p, long i, int dir);
void su2staple_wilson(fields f, params p, long i, int dir, double* V);
void su2link_staple(fields f, params p, long i, int dir, double* V);
// U(1) routines
double u1ptrace(fields f, params p, long i, int dir1, int dir2);
long double local_u1wilson(fields f, params p, long i);
double localact_u1link(fields f, params p, long i, int dir);
// doublet routines
double doubletsq(double* a);
long double avg_doublet2(fields f, params p);
long double avg_doublet4(fields f, params p);
double hopping_doublet_forward(fields f, params p, long i, int dir);
double hopping_doublet_backward(fields f, params p, long i, int dir);
double covariant_doublet(fields f, params p, long i);
double localact_doublet(fields f, params p, long i);
double higgspotential(fields f, params p, long i);
// triplet routines
double tripletsq(double* a);
double hopping_triplet_forward(fields f, params p, long i, int dir);
double hopping_triplet_backward(fields f, params p, long i, int dir);
double covariant_triplet(fields f, params p, long i);
double localact_triplet(fields f, params p, long i);

// metropolis.c
int metro_su2link(fields f, params p, long i, int dir);
int metro_u1link(fields f, params p, long i, int dir);
int metro_doublet(fields f, params p, long i);
int metro_triplet(fields f, params p, long i);

// heatbath.c
int heatbath_su2link(fields f, params p, long i, int dir);

// overrelax.c
double polysolve3(long double a, long double b, long double c, long double d);
int overrelax_doublet(fields f, params p, long i);
int overrelax_triplet(fields f, params p, long i);

// update.c
void update_lattice(fields* f, params p, comlist_struct* comlist, counters* c, char metro);
void update_lattice_muca(fields* f, params p, comlist_struct* comlist, weight* w, counters* c, char metro);
void checkerboard_sweep_su2link(fields f, params p, counters* c, char parity, int dir);
void checkerboard_sweep_u1link(fields f, params p, counters* c, char parity, int dir);
void checkerboard_sweep_su2doublet(fields f, params p, counters* c, char parity, char metro);
void checkerboard_sweep_su2triplet(fields f, params p, counters* c, char parity, char metro);
void sync_halos(fields* f, params p, comlist_struct* comlist);

// init.c
void setsu2(fields f, params p);
void random_su2link(double *su2);
void setu1(fields f, params p);
void setfields(fields f, params p);
void setdoublets(fields f, params p);
void settriplets(fields f, params p);
void init_counters(counters* c);

// checkpoint.c
void print_acceptance(params p, counters c);
void read_field(params p, FILE *file, double *field, int size);
void write_field(params p, FILE *file, double *field, int size);
void save_lattice(params p, fields f, counters c);
void load_lattice(params p, fields f, counters* c);

// parameters.c
void get_parameters(char *filename, params *p);
void get_weight_parameters(char *filename, params *p, weight* w);
void print_parameters(params p);
void read_updated_parameters(char *filename, params *p);


// measure.c
void measure(fields f, params p, counters* c, weight* w);
double action_local(fields f, params p, long i);
void print_labels();

// multicanonical.c
void load_weight(params p, weight *w);
void save_weight(params p, weight w);
double get_weight(weight w, double val);
void update_weight(params p, weight* w);
int multicanonical_acceptance(params p, weight* w, double oldval, double newval);
long whichbin(weight w, double val);
double calc_orderparam(params p, fields f, weight* w, char par);
void measure_muca(params p, fields f, weight* w);
void check_tunnel(params p, weight *w);
void store_muca_fields(params p, fields* f, weight* w);
void reset_muca_fields(params p, fields* f, weight* w, char par);
void alloc_backup_arrays(params p, fields* f, weight w);
void free_muca_arrays(fields* f, weight *w);


#endif

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
	ushort dim;
	uint *L;
	long vol;
	// slicing the lattice for MPI
	int *nslices; // how many slices in each direction
	uint *sliceL; // how many sites per slice in each direction
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

	// max iterations etc
	long iterations;
	long checkpoint;
	long interval;
	FILE *resultsfile;
	int run_checks;

	int multicanonical;

	// Parameters in the action
	// possible generalization: make separate structures
	// for each field in the theory that contain their respective params
	double betasu2; // 4/(g^2)
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
	short algorithm_su2doublet;
	short algorithm_su2triplet;

	// How many times to update a field per sweep
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

} fields;

typedef struct {

	// keep track of evaluation time
	double comms_time;

	// count metropolis updates
	long total_su2link;
	long total_doublet;
	long total_triplet;
	long accepted_su2link;
	long accepted_doublet;
	long accepted_triplet;
	// count overrelax updates
	long acc_overrelax_doublet;
	long total_overrelax_doublet;
	long acc_overrelax_triplet;
	long total_overrelax_triplet;

	long accepted_muca;
	long total_muca;
} counters;


// multicanonical weight: piecewise linear weight[i] = slope[i] * x + W[i]
typedef struct {
	long bins;
	double min, max;
	double dbin; // size (width) of one bin
	double* pos; // weight "position", i.e. values of the order param in the given range
	double* W; // value of the weight function at each position
	double* slope;
	double increment; // how much the weight is increased after visiting a bin
	int* visited; // keep track of which bins we have visited
	long visited_total;
	char readonly;
	char weightfile[100];
} weight;

// comms.c (move elsewhere later?)
void make_comlists(params *p, comlist_struct *comlist, long** xphys);
void reorder_sitelist(params* p, sendrecv_struct* sr);
void find_max_sendrecv(comlist_struct* comlist);
void reorder_comlist(params* p, comlist_struct* comlist);
double reduce_sum(double res);
double allreduce(double res);
void bcast_int(int *res);
void bcast_long (long *res);
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
long **alloc_latticetable(ushort dim, long sites);
void free_latticetable(long** list);
void free_lattice_arrays(params *p);
void free_comlist(comlist_struct* comlist);

// su2.c
double su2sqr(double *u);
void su2rot(double *u1, double *u2);
double su2ptrace(fields f, params p, long i, int dir1, int dir2);
long double local_su2wilson(fields f, params p, long i);
double localact_su2link(fields f, params p, long i, int dir);
void su2staple_wilson(fields f, params p, long i, int dir, double* V);
void su2link_staple(fields f, params p, long i, int dir, double* V);
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
int metro_doublet(fields f, params p, long i, weight* w);
int metro_doublet_transverse(fields f, params p, long i);
int metro_triplet(fields f, params p, long i);


// heatbath.c
int heatbath_su2link(fields f, params p, long i, int dir);

// overrelax.c
int overrelax_doublet(fields f, params p, long i);
int overrelax_triplet(fields f, params p, long i);

// update.c
void update_lattice(fields f, params* p, comlist_struct* comlist, weight* w, counters* c, char metro);
void checkerboard_sweep_su2link(fields f, params p, counters* c, char parity, int dir);
void checkerboard_sweep_su2doublet(fields f, params* p, weight* w, counters* c, char parity, char metro);
void checkerboard_sweep_su2triplet(fields f, params* p, counters* c, char parity, char metro);
void transverse_sweep_doublet(fields f, params p, counters* c, char parity);

// init.c
void setsu2(fields f, params p);
void random_su2link(double *su2);
void setfields(fields f, params p);
void setdoublets(fields f, params p);
void settriplets(fields f, params p);
void init_counters(counters* c);

// parameters.c
void get_parameters(char *filename, params *p);
void get_weight_parameters(char *filename, params *p, weight* w);
void print_parameters(params p);


// measure.c
void measure(fields f, params p, counters* c);
double action_local(fields f, params p, long i);
void print_labels();

// multicanonical.c
void load_weight(params p, weight *w);
void save_weight(params p, weight w);
double get_weight(weight w, double val);
void update_weight(params* p, weight* w, double val);
int multicanonical_acceptance(params* p, weight* w, double oldval, double newval);
long whichbin(weight w, double val);
double calc_orderparam(params p, fields f);
double calc_weight_slope(weight w, long bin);
void muca_param_shift(params* p, weight w, long bin, double oldslope);

#endif

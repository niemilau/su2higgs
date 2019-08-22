#ifndef SU2_H
#define SU2_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

// did we compile with MPI?
#ifdef MPI
#include <mpi.h>
#endif

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;

// degrees of freedom per site for different fields
#define SU2DB 4
#define SU2LINK 4
#define SU2TRIP 3

// update algorithms
#define METROPOLIS 1 // Metropolis
#define HEATBATH 2 // Heatbath
#define OVERRELAX 3 // Overrelaxation

typedef struct {

	// lattice dimensions
	ushort dim;
	uint *L;
	ulong vol;
	// neighboring sites
	//next[i][dir] gives the index of the next site after i in direction dir
	ulong **next;
	ulong **prev;
	// parity of a site is 0 if the physical index x+y+z+... is even, 1 otherwise
	char *parity;

	// max iterations etc
	ulong iterations;
	ulong checkpoint;
	ulong interval;
	FILE *resultsfile;

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

	// Stuff to measure
	double current_action;
	double current_wilson;
	double current_doublet2;
	double current_doublet4;
	double current_hopping_doublet;

	// How to update the fields
	short algorithm_su2link;
	short algorithm_su2doublet;
	short algorithm_su2triplet;

	// How many times to update a field per sweep
	short update_su2doublet;
	short update_su2triplet;


} params;

typedef struct {

	//double *su2singlet;
	double ***su2link;
	double **su2doublet;
	double **su2triplet;

} fields;

typedef struct {

	// count metropolis updates
	ulong total_su2link;
	ulong total_doublet;
	ulong total_triplet;
	ulong accepted_su2link;
	ulong accepted_doublet;
	ulong accepted_triplet;
	// count overrelax updates
	ulong acc_overrelax_doublet;
	ulong total_overrelax_doublet;
	ulong acc_overrelax_triplet;
	ulong total_overrelax_triplet;
} counters;


// alloc.c
double *make_singletfield(params p);
double **make_field(params p, int dofs);
double ***make_gaugefield(params p, int dofs);
void free_singletfield(params p, double *field);
void free_field(params p, double **field);
void free_gaugefield(params p, double ***field);
void alloc_fields(params p, fields *f);
void free_fields(params p, fields *f);
void alloc_neighbors(params *p);
ulong **alloc_neighborList(params p);
void free_neighbors(params *p);


// su2.c
double su2sqr(double *u);
void su2rot(double *u1, double *u2);
double su2ptrace(fields f, params p, ulong i, int dir1, int dir2);
long double local_su2wilson(fields f, params p, ulong i);
double localact_su2link(fields f, params p, ulong i, int dir);
void su2staple_wilson(fields f, params p, ulong i, int dir, double* V);
void su2link_staple(fields f, params p, ulong i, int dir, double* V);
// doublet routines
double doubletsq(double* a);
long double avg_doublet2(fields f, params p);
long double avg_doublet4(fields f, params p);
double hopping_doublet_forward(fields f, params p, ulong i, int dir);
double hopping_doublet_backward(fields f, params p, ulong i, int dir);
double covariant_doublet(fields f, params p, ulong i);
double localact_doublet(fields f, params p, ulong i);
double higgspotential(fields f, params p, ulong i);
// triplet routines
double tripletsq(double* a);
double hopping_triplet_forward(fields f, params p, ulong i, int dir);
double hopping_triplet_backward(fields f, params p, ulong i, int dir);
double covariant_triplet(fields f, params p, ulong i);
double localact_triplet(fields f, params p, ulong i);

// metropolis.c
int metro_su2link(fields f, params p, ulong i, int dir);
int metro_doublet(fields f, params p, ulong i);
int metro_triplet(fields f, params p, ulong i);

// heatbath.c
int heatbath_su2link(fields f, params p, ulong i, int dir);

// overrelax.c
int overrelax_doublet(fields f, params p, ulong i);
int overrelax_triplet(fields f, params p, ulong i);

// update.c
void checkerboard_sweep(fields f, params p, counters* c, int dir, char parity);
void update_lattice(fields f, params p, counters* c, char metro);
void checkerboard_sweep_su2link(fields f, params p, counters* c, char parity, int dir);
void checkerboard_sweep_su2doublet(fields f, params p, counters* c, char parity, char metro);
void checkerboard_sweep_su2triplet(fields f, params p, counters* c, char parity, char metro);

// init.c
void indexToCoords(params p, ulong i, uint* x);
ulong coordsToIndex(params p, uint* x);
void calculate_neighbors(params *p);
void setsu2(fields f, params p);
void random_su2link(double *su2);
void setfields(fields f, params p);
void setdoublets(fields f, params p);
void settriplets(fields f, params p);

// parameters.c
void get_parameters(char *filename, params *p);
void print_parameters(params p);


// measure.c
void measure(fields f, params p);
double action_local(fields f, params p, ulong i);
void print_labels();

#endif

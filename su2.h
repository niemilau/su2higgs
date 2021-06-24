#ifndef SU2_H
#define SU2_H

#include "comms.h"
#include "generic/stddefs.h"

#ifndef NHIGGS // specify in makefile, otherwise assume no Higgs doublets
	#define NHIGGS 0
#endif

#ifdef U1
	// nasty global for Higgs hypercharge
	extern const double higgs_Y; // set Y=1/2 in parameters.c
#endif

/* Struct "lattice": contains info on lattice dimensions, lookup tables for
* sites and everything related to parallelization. */
typedef struct {
	// layout parameters for MPI
	int rank, size;

	MPI_Comm comm;

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
	// coords[i][dir] = x_dir coordinate on the full lattice of site i.
	long **coords;
	long *offset; // coordinate offset in my node wrt. the full lattice
	long firstsite; // index of (x,y,z)=(0,0,0), with (x,y,z) being coordinates on the node

	// parity of a site is EVEN if the physical index x+y+z+... is even, ODD otherwise
	char *parity;
	long evensites, oddsites;
	long evenhalos, oddhalos;
	comlist_struct comlist;
	// in layout.c we reorder lattice sites so that EVEN sites come first.
	int reorder_parity; // for debugging purposes

	/* miscellaneous info about the lattice, used for example in correlation.c.
	* Alloc'd and filled in by make_misc_tables() in layout.c */
	int longest_dir; // longest direction on the full lattice
	long* sites_per_coord; // how many sites for each coordinate x_dir
	long*** sites_at_coord; // list of sites for a fixed x_dir (sites_at_coord[dir][x][i])

	#ifdef BLOCKING
		// communications between the blocked lattice and the original
		comlist_struct blocklist;
		// some MPI nodes may not fit on the blocked lattice and have to standby
		int standby;
		int blocking_level;
	#endif

	#ifdef MEASURE_Z
		/* Single out a particular direction ("z coordinate") to study e.g. wall profile, correlation lengths.
		* These are initialized in init_z_coord() which is called by layout() */
		int z_dir; // specify the z direction, default = longest direction
		long sites_per_z, offset_z; // offset = physical value of z at the "first" site in my node
		long area; // how many sites per z coordinate on the whole lattice; area = vol / L[z_dir]
		// list of all real sites at fixed z coord in no particular order
		long** site_at_z; // sites_at_z[z][j], 0<=j<sites_per_z
	#endif
} lattice;


/* struct "params": contains control parameters such as number of iterations,
* and also values for model-dependent parameters in the action */
typedef struct {

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

	/* Parameters in the action */
	double betasu2; // 4 / (a g^2)


	// doublet
	double msq_phi; // Higgs mass
	double lambda_phi; // Higgs quartic

	// triplet
	double msq_triplet;
	double a2; // Higgs portal
	double b4; // self quartic

	#if (NHIGGS == 2)
		/* with 2 Higgs doublets, the couplings are assumed to be in a doublet basis
		* where the kinetic terms are diagonal and canonically normalized */
		double msq_phi2;
		complex m12sq;
		double lam2, lam3, lam4;
		complex lam5, lam6, lam7;
	#endif

	// initial values for fields
	double phi0; // doublet
	double sigma0; // triplet

	// How to update the fields
	short algorithm_su2link;
	short algorithm_su2doublet;
	short algorithm_su2triplet;

	// How many times to update a field per sweep
	short update_links;
	short scalar_sweeps; // update all scalars n times per iteration
	// additional sweeps on top of scalar_sweeps
	short update_su2doublet;
	short update_su2triplet;

	// U(1) hypercharge
	#ifdef U1
		// Wilson action is S = gamma^2 * betau1 * sum_{x, i<j} (1 - Re p_{ij}^{1/gamma})
		double betau1; // 1 / (a g'^2) in this code
		double gammau1; // my gamma = 2 * gamma of hep-lat/9705003 and hep-lat/9612006
		short algorithm_u1link;
	#endif

	#ifdef SINGLET
		int update_singlet;
		short algorithm_singlet;
		double singlet0;
		double b1_s, msq_s, b3_s, b4_s;
		#if (NHIGGS == 1)
			double a1_s, a2_s;
		#endif
	#endif

	#ifdef CORRELATORS
		int do_correlators;
		int correlator_interval; // how often to measure correlators
	#endif

	#ifdef BLOCKING
		int blocks; // how many times to block the lattice for correlators
	#endif

	#ifdef MEASURE_Z
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

	#ifdef HB_TRAJECTORY
		int do_trajectory; // measure realtime trajectories if nonzero, specify in config file
	#endif

} params;


// Field structure. Note that this can safely be passed by value to updating functions
// after the contents have been alloc'd
typedef struct {

	//double *su2singlet;
	double ***su2link;
	double **su2triplet;

	#ifdef U1
		double **u1link; // actually the exponent
	#endif

	// NHIGGS copies of Higgs doublets, accessed as su2doublet[id][site][component]
	#if (NHIGGS > 0)
		double** su2doublet[NHIGGS];
	#endif

	#ifdef SINGLET
		double** singlet; // access as singlet[site][0], using 2d table so that this works with comms routines
	#endif

} fields;


typedef struct {

	long iter;

	// keep track of when a metropolis sweep should be forced
	int higgs_sweeps, triplet_sweeps;
	#ifdef SINGLET
		int singlet_sweeps;
	#endif

	// count metropolis updates
	long total_su2link, accepted_su2link;
	long total_u1link, accepted_u1link;
	#if (NHIGGS > 0)
		long total_doublet[NHIGGS], accepted_doublet[NHIGGS];
		long total_overrelax_doublet[NHIGGS], acc_overrelax_doublet[NHIGGS];
	#endif

	#ifdef SINGLET
		long accepted_singlet;
		long total_singlet;
		long acc_overrelax_singlet;
		long total_overrelax_singlet;
	#endif

	long total_triplet, accepted_triplet;
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
	* Will not update weight either if set to 0. Used in update sweep routines. */
	int do_acceptance;

	// how many times is multicanonical_acceptance() called per update sweep (separately for each field, parity)
	int checks_per_sweep;

	int bins;
	double min, max; // weighting range
	double wrk_min, wrk_max; // range of orderparam values where weight is to be updated
	//long min_bin, max_bin; // indices of the bins containing w.wrk_min and w.wrk_max
	double* pos; // weight "position", i.e. values of the order param in the given range
	double* W; // value of the weight function at the beginning of each bin
	double* slope; double* b; // linearized weight: W(R) = w_i + (R - R_i)*slope[i] = b + slope*R

	double delta; // how much the weight is increased in update_weight()
	int* hits; // keep track of which bins we have visited
	int muca_count; // how many muca acc/rej steps performed (resets after weight update)
	int update_interval; // how many muca acc/rej steps until weight is updated

 	int last_max; // 1 if system recently visited the bin containing w.wrk_max (keep track of tunneling)
	int mode; // one of the multicanonical "modes" defined above, affects weight recursion
	char weightfile[100]; // file name

	fields fbu; // backup fields for undoing a rejected multicanonical update sweep

	// additional data arrays used in slow update mode only
	long* gsum;
	long* nsum;
	double* hgram;
} weight;


// comms.c
void make_comlists(lattice *l, comlist_struct *comlist);
int addto_comlist(comlist_struct* comlist, int rank, long i, int sendrecv, char evenodd, long init_max);
void reorder_sitelist(lattice* l, sendrecv_struct* sr);
void reorder_comlist(lattice* l, comlist_struct* comlist);
double reduce_sum(double res, MPI_Comm comm);
double allreduce(double res, MPI_Comm comm);
long reduce_sum_long(long res, MPI_Comm comm);
void bcast_int(int *res, MPI_Comm comm);
void bcast_long (long *res, MPI_Comm comm);
void bcast_double(double *res, MPI_Comm comm);
void bcast_int_array(int *arr, int size, MPI_Comm comm);
void barrier(MPI_Comm comm);
// gauge links:
void update_gaugehalo(lattice* l, char parity, double*** field, int dofs, int dir);
#ifdef MPI
void send_gaugefield(sendrecv_struct* send, MPI_Comm comm, MPI_Request* req, char parity, double*** field, int dofs, int dir);
void recv_gaugefield(sendrecv_struct* recv, MPI_Comm comm, char parity, double*** field, int dofs, int dir);
#endif
// non-gauge fields:
void update_halo(lattice* l, char parity, double** field, int dofs);
#ifdef MPI
void recv_field(sendrecv_struct* recv, MPI_Comm comm, char parity, double** field, int dofs);
void send_field(sendrecv_struct* send, MPI_Comm comm, MPI_Request* req, char parity, double** field, int dofs);
#endif
void test_comms(lattice* l);
void test_comms_individual(lattice* l);

// layout.c
void layout(lattice *l, int do_prints, int run_checks);
void make_slices(lattice *l, int do_prints);
void sitemap(lattice *l);
void set_parity(lattice *l);
void paritymap(lattice* l, long* newindex);
void make_misc_tables(lattice* l);
void remap_latticetable(lattice* l, long** arr, long* newindex, long maxindex);
void remap_neighbor_table(lattice* l, long** arr, long* newindex, long maxindex);
void remap_lattice_arrays(lattice* l, long* newindex, long maxindex);
long findsite(lattice const* l, long* x, int include_halos);
void test_coords(lattice const* l);
void test_neighbors(lattice const* l);
void indexToCoords(short dim, int* L, long i, long* x);
long coordsToIndex(short dim, int* L, long* x);
int coordsToRank(lattice const* l, long* coords);
void printf0(lattice l, char *msg, ...);
void die(int howbad);
void print_lattice_2D(lattice *l);

// alloc.c
double *make_singletfield(long sites);
double **make_field(long sites, int dofs);
double ***make_gaugefield(long sites, int dim, int dofs);
void free_singletfield(double *field);
void free_field(double **field);
void free_gaugefield(long sites, double ***field);
void alloc_fields(lattice const* l, fields *f);
void free_fields(lattice const* l, fields *f);
void alloc_lattice_arrays(lattice *l, long sites);
long **alloc_latticetable(int dim, long sites);
long **realloc_latticetable(long** arr, int dim, long oldsites, long newsites);
void realloc_lattice_arrays(lattice *l, long oldsites, long newsites);
void alloc_comlist(comlist_struct* comlist, int nodes);
void realloc_comlist(comlist_struct* comlist, int sendrecv);
void free_latticetable(long** list);
void free_lattice(lattice *l);
void free_comlist(comlist_struct* comlist);

// su2u1.c
double su2sqr(double *u);
void su2rot(double *u1, double *u2);
void su2plaquette(lattice const* l, fields const* f, long i, int dir1, int dir2, double* u1);
double su2ptrace(lattice const* l, fields const* f, long i, int dir1, int dir2);
long double local_su2wilson(lattice const* l, fields const* f, params const* p, long i);
double localact_su2link(lattice const* l, fields const* f, params const* p, long i, int dir);
double su2trace4(double *u1, double *u2, double *u3, double *u4);
void clover_su2(lattice const* l, fields const* f, long i, int d1, int d2, double* clover);
double hopping_trace(double* phi1, double* u, double* phi2);
double hopping_trace_su2u1(double* phi1, double* u, double* phi2, double a);
double hopping_trace_triplet(double* a1, double* u, double* a2);
#ifdef U1
// U(1) routines
double u1ptrace(lattice const* l, fields const* f, long i, int dir1, int dir2);
double local_u1wilson(lattice const* l, fields const* f, params const* p, long i);
double localact_u1link(lattice const* l, fields const* f, params const* p, long i, int dir);
#endif
// doublet routines
double doubletsq(double* a);
void phiproduct(double* f1, double const* f2, int conj);
complex get_phi12(double const* h1, double const* h2);
double hopping_doublet_forward(lattice const* l, fields const* f, long i, int dir, int higgs_id);
double hopping_doublet_backward(lattice const* l, fields const* f, long i, int dir, int higgs_id);
double covariant_doublet(lattice const* l, fields const* f, long i, int higgs_id);
double localact_doublet(lattice const* l, fields const* f, params const* p, long i, int higgs_id);
double higgspotential(fields const* f, params const* p, long i);
// triplet routines
double tripletsq(double* a);
double hopping_triplet_forward(lattice const* l, fields const* f, params const* p, long i, int dir);
double hopping_triplet_backward(lattice const* l, fields const* f, params const* p, long i, int dir);
double covariant_triplet(lattice const* l, fields const* f, params const* p, long i);
double localact_triplet(lattice const* l, fields const* f, params const* p, long i);
#ifdef SINGLET
double localact_singlet(lattice const* l, fields const* f, params const* p, long i);
double potential_singlet(fields const* f, params const* p, long i);
#endif
#ifdef BLOCKING
// smearing
void smear_link(lattice const* l, fields const* f, int const* smear_dir, double* res, long i, int dir);
void smear_triplet(lattice const* l, fields const* f, int const* smear_dir, double* res, long i);
void smear_fields(lattice const* l, fields const* f, fields* f_b, int const* block_dir);
#endif


// staples.c
void su2staple_counterwise(double* V, double* u1, double* u2, double* u3);
void su2staple_clockwise(double* V, double* u1, double* u2, double* u3);
void su2staple_wilson(lattice const* l, fields const* f, long i, int dir, double* V);
void su2staple_wilson_onedir(lattice const* l, fields const* f, long i, int mu, int nu, int dagger, double* res);
void su2link_staple(lattice const* l, fields const* f, params const* p, long i, int dir, double* V);
void staple_doublet(double* res, lattice const* l, fields const* f, params const* p, long i, int higgs_id);


// metropolis.c
int metro_su2link(lattice const* l, fields* f, params const* p, long i, int dir);
int metro_u1link(lattice const* l, fields* f, params const* p, long i, int dir);
int metro_doublet(lattice const* l, fields* f, params const* p, long i, int higgs_id);
int metro_triplet(lattice const* l, fields* f, params const* p, long i);
#ifdef SINGLET
int metro_singlet(lattice const* l, fields* f, params const* p, long i);
#endif

// heatbath.c
int heatbath_su2link(lattice const* l, fields* f, params const* p, long i, int dir);

// overrelax.c
double polysolve3(long double a, long double b, long double c, long double d);
#if (NHIGGS > 0)
int overrelax_doublet(lattice const* l, fields* f, params const* p, long i); // for N=1 Higgs potentials
int overrelax_higgs2(lattice const* l, fields* f, params const* p, long i, int higgs_id); // for N>1 Higgs potentials
#endif
int overrelax_triplet(lattice const* l, fields* f, params const* p, long i);
#ifdef SINGLET
int overrelax_singlet(lattice const* l, fields* f, params const* p, long i);
#endif


// update.c
void update_lattice(lattice* l, fields* f, params const* p, counters* c, weight* w);
void checkerboard_sweep_su2link(lattice const* l, fields* f, params const* p, counters* c, int parity, int dir);
void checkerboard_sweep_u1link(lattice const* l, fields* f, params const* p, counters* c, int parity, int dir);
int checkerboard_sweep_su2doublet(lattice const* l, fields* f, params const* p, counters* c,
			weight* w, int parity, int metro, int higgs_id);
int checkerboard_sweep_su2triplet(lattice const* l, fields* f, params const* p, counters* c, weight* w, int parity, int metro);
#ifdef SINGLET
int checkerboard_sweep_singlet(lattice const* l, fields* f, params const* p, counters* c,
			weight* w, int parity, int metro);
#endif
void sync_halos(lattice* l, fields* f);
int muca_check(lattice const* l, fields* f, params const* p, counters* c, weight* w, int parity);
void shuffle(int *arr, int len);

// init.c
void setsu2(fields* f, lattice const* l);
void random_su2link(double *su2);
void setu1(fields* f, lattice const* l);
#ifdef SINGLET
void set_singlets(fields* f, lattice const* l, params const* p);
#endif
void setfields(fields* f, lattice* l, params const* p);
void setdoublets(fields* f, lattice const* l, params const* p);
void settriplets(fields* f, lattice const* l, params const* p);
void cp_field(lattice const* l, double** field, double** new, int dofs, int parity);
void copy_fields(lattice const* l, fields const* f_old, fields* f_new);
void init_counters(counters* c);

// checkpoint.c
void print_acceptance(params p, counters c);
void write_field(lattice const* l, FILE *file, double *field, int size);
void read_field(lattice const* l, FILE *file, double *field, int size);
void save_lattice(lattice const* l, fields f, counters c, char* fname);
void load_lattice(lattice* l, fields* f, counters* c, char* fname);

// parameters.c
void check_set(int set, char *name); // helper routine
void get_parameters(char *filename, lattice* l, params *p);
void get_weight_parameters(char *filename, lattice const* l, params *p, weight* w);
void print_parameters(lattice l, params p);
void read_updated_parameters(char *filename, lattice const* l, params *p, weight* w);


// measure.c
void measure(FILE* file, lattice const* l, fields const* f, params const* p, weight* w);
// weight not necessarily constant because measure() can recalculate the order parameter
double action_local(lattice const* l, fields const* f, params const* p, long i);
void print_labels();
void measure_local(char* fname, lattice const* l, fields const* f, params const* p);
void print_labels_local(lattice const* l, char* fname);

// multicanonical.c
void load_weight(lattice const* l, weight *w);
void load_weight_params(int rank, char* fname, weight* w);
void save_weight(lattice const* l, weight const* w);
void linearize_weight(weight* w);
double get_weight(weight const* w, double val);
void muca_accumulate_hits(weight* w, double val);
int update_weight(weight* w);
int multicanonical_acceptance(lattice const* l, weight* w, double oldval, double newval);
int whichbin(weight const* w, double val);
double calc_orderparam(lattice const* l, fields const* f, params const* p, weight* w, char par);
void alloc_muca_backups(lattice const* l, weight* w);
void free_muca_arrays(fields* f, weight *w);
void init_last_max(weight* w);
void update_weight_slow(weight* w);

#ifdef CORRELATORS
// correlation.c
	#if (NHIGGS > 0)
		double higgs_correlator(lattice* l, fields const* f, int x, int dir, int higgs_id);
	#endif
	#ifdef TRIPLET
		double triplet_correlator(lattice* l, fields const* f, int d, int dir);
		complex projected_photon_operator(lattice* l, fields const* f, int z, int dir, int* mom);
		complex projected_photon_correlator(lattice* l, fields const* f, int d, int dir);
		double projected_photon_operator_old(lattice* l, fields const* f, params const* p, int z, int dir,
		      int* mom, double* res_re, double* res_im);
		void projected_photon_correlator_old(lattice* l, fields const* f, params const* p, int d,
		  		int dir, double* res_re, double* res_im);
	#endif
	double plane_sum(double (*funct)(double*), double** field, lattice* l, long x, int dir);
	void measure_correlators(char* fname, lattice* l, fields const* f, params const* p, int dir, int meas_id);
	void print_labels_correlators();
	#ifdef BLOCKING
	void measure_blocked_correlators(lattice* l, lattice* b, fields const* f, fields* f_b, params const* p,
				int const* block_dir, int dir, int id);
	#endif
#endif

#ifdef MEASURE_Z
	// z_coord.c
	void init_z_coord(lattice* l);
	void print_z_labels(lattice const* l, params* p);
	void measure_along_z(lattice const* l, fields const* f, params const* p, long id);
	void prepare_wall(lattice* l, fields* f, params const* p);
#endif


#ifdef TRIPLET
	// magfield.c
	void matmat(double *in1, double *in2, int dag);
	void projector(double *proj, double *adjoint);
	void project_u1(lattice const* l, fields const* f, long i, int dir, double* pro);
	double alpha_proj(lattice const* l, fields const* f, params const* p, long i, int dir1, int dir2);
	double magfield(lattice const* l, fields const* f, params const* p, long i, int dir);
	double magcharge_cube(lattice const* l, fields const* f, params const* p, long i);
#endif


#ifdef GRADFLOW
	// gradflow.c
	void grad_flow(lattice* l, fields const* f, params* p, weight* w, double t_max, double dt, int flow_id);
	void grad_force_link(lattice const* l, fields const* f, params const* p, double* force, long i, int dir);
	void grad_force_triplet(lattice const* l, fields const* f, params const* p, double* force, long i);
	void calc_gradient(lattice const* l, fields const* f, params const* p, fields* forces);
	void flow_gauge(lattice const* l, fields* flow, fields const* forces, double dt);
	void flow_triplet(lattice const* l, fields* flow, fields const* forces, double dt);
	void remove_counterterms(params* p);
#endif

#ifdef BLOCKING
	// blocking.c
	int max_block_level(lattice const* l, int const* block_dir);
	void block_lattice(lattice* l, lattice* b, int const* block_dir);
	void make_blocklists(lattice* l, lattice* b, int const* block_dir);
	void standby_layout(lattice* l);
	void transfer_blocked_gaugefield(lattice* l, lattice* b, double*** field, double*** field_b, int dofs, int dir);
	void transfer_blocked_field(lattice* l, lattice* b, double** field, double** field_b, int dofs);
	void make_blocked_fields(lattice* l, lattice* b, fields const* f, fields* f_blocked);
	void block_fields_ownnode(lattice const* l, lattice const* b, fields const* f_smeared, fields* f_blocked);
	void test_blocking(lattice* l, lattice* b, int const* block_dir);
#endif


#ifdef HB_TRAJECTORY
	// realtime trajectories with heatbath algorithm

	typedef struct {
	  FILE *trajectoryfile;
	  int n_traj; // how many trajectories to generate from each initial configuration
	  int interval; // how often to measure in trajectory mode
	  int mode_interval; // how often to switch to trajectory mode
	  double min, max; // min and max values of the order parameter before trajectory is considered complete
	} trajectory;

	// hb_trajectory.c
	void make_realtime_trajectories(lattice* l, fields const* f, params* p, counters* c,
		 	weight* w, trajectory* traj, int id);
	void read_realtime_config(char *filename, trajectory* traj);
	void print_realtime_params(params p, trajectory traj);

#endif

#endif // end #ifndef SU2_H

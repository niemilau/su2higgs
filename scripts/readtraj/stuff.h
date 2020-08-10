
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

typedef struct {
  double* t; // time
  double* val;
  int len;
} traj;

typedef struct {
  int ntraj;
  traj* traj;
} traj_set;

inline long factorial(int n) {
  if (n==0) return 1;
  else if (n < 0) {
    printf("!!! Factorial error !\n");
    return 1;
  }
  long res = 1;
  for (int i=2; i<n; i++) res *= i;
  return n*res;
}

void die(int howbad);
void realtime(char* fname, int sample, int blocks);
void read_trajectory(FILE* f, traj* t);
void read_traj_set(char* fname, traj_set* set, int id);
double dynamical_prefactor(traj_set const* set, double min, double max, double phi_c);
void jackknife_double(double* arr, long datalen, int blocks, double* res);
void dynamics(traj_set* set, double min, double max, double phi_c, double* res);

void make_traj(traj* t, int len);
void realloc_traj(traj* t, int newlen);
void free_traj(traj* t);
void realloc_traj_set(traj_set* set, int newlen);
void free_traj_set(traj_set* set);
void make_traj_set(traj_set* set, int ntrajs);

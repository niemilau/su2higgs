/* A C-program for MT19937: Real number version  (1998/4/6)    */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

/** Modified somewhat by L. Niemi 2021: removed global definitions of N and M
* and moved some definitions to mersenne.h. Included integer generator too. **/

#include<stdlib.h>
#include<stdio.h>
#include "stddefs.h"


static unsigned int mt[MERSENNE_N]; /* the array for the state vector  */
int mersenne_i = -1; /*  < 0 means mt[N] is not initialized */
int mersenne_i_int = -1;
double mersenne_array[MERSENNE_N];
unsigned int mersenne_array_int[MERSENNE_N];

/* initializing the array with a NONZERO seed */
void
seed_mersenne(long seed)
{
  int mti;
  mt[0]= seed & 0xffffffffUL;
  for (mti=1; mti<MERSENNE_N; mti++) {
    mt[mti] =
      (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
  mersenne_i = 0;
  mersenne_i_int = 0;
}

double /* generating reals */
/* unsigned int */ /* for integer generation */
mersenne_generate(int *dummy)
{
  register unsigned int y;
  register int kk;
  static unsigned int mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mersenne_i < 0) {  /* if sgenrand() has not been called, */
    printf("!!! Mersenne generator not seeded!\n");
    exit(0);
  }

  /* generate MERSENNE_N words at one time */

  for (kk=0;kk<MERSENNE_N - MERSENNE_M;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk + MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
  }
  for (;kk<MERSENNE_N - 1;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+(MERSENNE_M - MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
  }
  y = (mt[MERSENNE_N - 1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
  mt[MERSENNE_N - 1] = mt[MERSENNE_M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

  for (kk=0; kk<MERSENNE_N; kk++) {
    y = mt[kk];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    mersenne_array[kk] = (double)y * 2.3283064365386963e-10;  /* reals: interval [0,1) */
  }

  mersenne_i = MERSENNE_N;
  return ( mersenne_array[--mersenne_i] );
    /* return y; */ /* for integer generation */
}

/* Sme as above but generates integers instead. Range: [0, 2^32−1] */
unsigned int mersenne_generate_int(int *dummy)
{
  register unsigned int y;
  register int kk;
  static unsigned int mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mersenne_i_int < 0) {  /* if sgenrand() has not been called, */
    printf("!!! Mersenne generator not seeded!\n");
    exit(0);
  }

  /* generate MERSENNE_N words at one time */

  for (kk=0;kk<MERSENNE_N - MERSENNE_M;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk + MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
  }
  for (;kk<MERSENNE_N - 1;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+(MERSENNE_M - MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
  }
  y = (mt[MERSENNE_N - 1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
  mt[MERSENNE_N - 1] = mt[MERSENNE_M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

  for (kk=0; kk<MERSENNE_N; kk++) {
    y = mt[kk];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    mersenne_array_int[kk] = y; /* integers: interval [0, 2^32−1] */
  }

  mersenne_i_int = MERSENNE_N;
  return ( mersenne_array_int[--mersenne_i_int] );
}

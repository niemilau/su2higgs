#ifndef MERSENNE_H
#define MERSENNE_H

/***************************************************************
 *  mersenne.h
 * This header just declares the Mersenne Twister functions and defines macros.
 */

double genrand64_real2(void);
void init_genrand64(unsigned long long seed);
unsigned long long genrand64_int64(void);

#define seed_mersenne(seed) init_genrand64(seed)
#define dran() genrand64_real2() /* random double on interval [0,1) */
#define iran() genrand64_int64() /* random unsigned long long integer on [0, 2^64-1]-interval */
/* The integer generator is total overkill for my purposes, but whatever */

#endif

#ifndef MERSENNE_H
#define MERSENNE_H

/***************************************************************
 *  mersenne.h
 *  for the inline version of the mersenne generator
 */

 /* Period parameters */
 #define MERSENNE_N 624
 #define MERSENNE_M 397
 #define MATRIX_A 0x9908b0df   /* constant vector a */
 #define UPPER_MASK 0x80000000 /* most significant w-r bits */
 #define LOWER_MASK 0x7fffffff /* least significant r bits */

 /* Tempering parameters */
 #define TEMPERING_MASK_B 0x9d2c5680
 #define TEMPERING_MASK_C 0xefc60000
 #define TEMPERING_SHIFT_U(y)  (y >> 11)
 #define TEMPERING_SHIFT_S(y)  (y << 7)
 #define TEMPERING_SHIFT_T(y)  (y << 15)
 #define TEMPERING_SHIFT_L(y)  (y >> 18)

extern int mersenne_i;
extern int mersenne_i_int;
extern double mersenne_array[MERSENNE_N];
extern unsigned int mersenne_array_int[MERSENNE_N];

/* mersenne_generate() produces MERSENNE_N numbers and stores them in mersenne_array. Here I take one such number
* from the array, unless I already used up all of them, in which case generate new numbers and proceed */
#define mersenne() ( mersenne_i > 0 ? mersenne_array[--mersenne_i] : mersenne_generate(&mersenne_i) )
#define mersenne_int() ( mersenne_i_int > 0 ? mersenne_array_int[--mersenne_i_int] : mersenne_generate_int(&mersenne_i_int) )

#define dran() mersenne() /* random double in interval [0,1) */
#define iran() mersenne_int() /* random integer in interval [0, 2^32−1] */

void seed_mersenne(long a);
double mersenne_generate(int *);
unsigned int mersenne_generate_int(int *);


#endif

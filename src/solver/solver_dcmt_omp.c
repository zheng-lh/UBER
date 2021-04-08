/* 
      Gaussian pseudo-random number generator from Dynamic Creator of Mersenne
      Twisters (DCMT) with OpenMP parallelization
      
      Generating multiple independent streams of double precision pseudo-random
      numbers in Gaussian distribution. This source file also serves as a wrapper
      of the C library DCMT for Fortran applications.
      
      Author: Liheng Zheng
              Department of Physics and Astronomy
              Rice University
      2014

      Exponential random number generator added, to be used in implementing Gobet
      [2001] half-space reflection algorithm.

      Liheng Zheng, 2016
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>
#include "dc.h"

/* 
 Exponent of the period is a Mersenne prime, it could be 521, 607, 1279, 2203,
 2281,... The period of the PRN's is then 2^p - 1. For large p, the time needed to
 find an mt_struct increases as O(p^3) (See README in the DCMT library).
 */   
static int w = 32;
static int p = 2203;
static uint32_t mt_seed = 863;
static mt_struct *mts;
#pragma omp threadprivate(mts)

void get_mt_parameter_id_st_f_(int *id);
void sgenrand_mt_f_(int *seed);
void free_mt_struct_f_(void);
double unidev(mt_struct *mts);
double gasdev(mt_struct *mts);
double expdev(mt_struct *mts);
double unirand_(void);
double gausrand_(void);
double exprand_(void);

/*
The functions whose names are suffixed by underscores will be accessed by Fortran.
*/

void get_mt_parameter_id_st_f_(int *id)
{   /* Find the mt_struct parameters for thread *id. */ 
    mts = get_mt_parameter_id_st(w, p, *id, mt_seed);
    if( mts==NULL ){
       printf(" Failed to get mt_struct parameter for id = %d\n",*id);
       printf("  w = %d, p = %d, mt_seed = %d\n", w, p, mt_seed);
       exit(1);
    }
}



void sgenrand_mt_f_(int *seed)
{   /* initialize the psuedo-random number generator */
    uint32_t useed;

    if( *seed<0 ){
       printf(" Warning: seed = %d < 0 passed to sgenrand_mt_f.\n",*seed);
    }
    useed = *seed;
    sgenrand_mt(useed, mts);
}



void free_mt_struct_f_(void)
{   /* Release the memory claimed by mts. */
    free_mt_struct(mts);
}



double unidev(mt_struct *mts)
{   /*
    Obtaining double precision (53-bit resolution) uniform deviates in [0,1] from
    32-bit unsigned random integers. Code was adapted from FORTRAN function
    mt_genrand_double1 in Ken-Ichi Ishikawa's Multiple Stream Mersenne Twister PRNG:
    <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>
    */
    double r;
    double a, b;
    uint32_t ia, ib;
    
    ia = genrand_mt(mts);    /* ia in [0,0xFFFFFFFF] */
    ib = genrand_mt(mts);    /* ib in [0,0xFFFFFFFF] */
    ia >>= 5;                /* ia in [0,2^27-1] */
    ib >>= 6;                /* ib in [0,2^26-1] */
    a = (double) ia;
    b = (double) ib;
    /* ( a*2^26 + b ) in [0,2^53-1]
       r = ( a*2^26 + b )/(2^53-1) */
    r= (a*67108864.0 + b)*(1.0/9007199254740991.0);
    return r;
}



double gasdev(mt_struct *mts)
{   /* 
    Normal distribution random number generator modified from the one on p.280 in 
    "Numerical Recipes in FORTRAN", second edition, by W. H. Press et al., Cambridge
    University Press, 1992.
    */
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    
    if (iset == 0){
       /* we don't have an extra deviate handy */
       do{
          /* so pick two uniform uumbers in the square extending from -1 to 1 in
             each direction */
          v1 = 2.0*unidev(mts) - 1.0;
          v2 = 2.0*unidev(mts) - 1.0;
          rsq = v1*v1 + v2*v2;
          /* see if they are in the unit circle, and if not, try again */
       } while (rsq>=1.0 || rsq==0.0);
       
       /* now make the Box-Muller transformation to get two normal deviates. return
          one and save the other for the next time. */
       fac = sqrt( -2.0*log(rsq)/rsq );       
       gset = v1*fac;
       /* set flag */
       iset = 1;
       return v2*fac;
    }
    else{
       /* we have an extra deviate handy, so return it, and unset the flag */
       iset = 0;
       return gset;
    }
}



double expdev(mt_struct *mts)
{   /* Exponential distribution random number with p(x) = exp(-x). */
    double x, y;

    do{
       x = unidev(mts);
    } while (x==0.0);

    y = -log(x);
    return y;
}



double unirand_(void)
{   /* Returns a double precision uniform deviate in [0,1] upon every call. */
    double unirand = unidev(mts);
    return unirand;
}



double gausrand_(void)
{   /* Returns a double precision normal deviate upon every call.  */
    double gausrand = gasdev(mts);
    return gausrand;
}



double exprand_(void)
{   /* Returns a double precision exponential deviate upon every call. */
    double exprand = expdev(mts);
    return exprand;
}


#ifdef TEST
void dcuint32_(uint32_t intarr[])
{   /* Test function to return count number of unsigned 32-bit random integers. */
    int i;    
    for( i=0; i<count; i++ )
        intarr[i] = genrand_mt(mtss[i]);
}



void dcunidev_(double uniarr[])
{   /* Test function to return count number of double precision uniform deviates. */
    int i;
    for( i=0; i<count; i++ )
        uniarr[i] = unidev(mtss[i]);
}
#endif











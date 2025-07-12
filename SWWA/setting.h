#ifndef _SETTING_H
#define _SETTING_H

#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include<gsl/gsl_math.h>

//#define AUX_INACTIVE
//#define BOS_INACTIVE

/******************************************************************************
 * Begin macros
 ******************************************************************************/
 enum BOOL {FALSE, TRUE};
 enum PARITY {EVEN, ODD, BOTH};
 #define MAXLENGTH 120
 #define THRESHOLD_OMP 64
/******************************************************************************
  * new data types to help with readibility and permit
  * easy swapping between single and double precision
 ******************************************************************************/
typedef size_t Uint;
#ifndef PRECISION
    #define PRECISION 2
#endif
#define L1ALIGN 64
#if 0
//#if ( PRECISION == 1 )
    typedef float Scalar __attribute__ ((aligned (L1ALIGN)));
    typedef float complex Complex __attribute__ ((aligned (L1ALIGN)));
    typedef float* Vector __attribute__ ((aligned (L1ALIGN)));
    typedef float* Sfield __attribute__ ((aligned (L1ALIGN)));
    typedef float** Cfield __attribute__ ((aligned (L1ALIGN)));
    typedef float** Vfield __attribute__ ((aligned (L1ALIGN)));
    typedef float** V_of_fields __attribute__ ((aligned (L1ALIGN)));
    #define MACH_EPS FLT_EPSILON
#else
    typedef double Scalar __attribute__ ((aligned (L1ALIGN)));
    typedef double complex Complex __attribute__ ((aligned (L1ALIGN)));
    typedef double* Vector __attribute__ ((aligned (L1ALIGN)));
    typedef double* Sfield __attribute__ ((aligned (L1ALIGN)));
    typedef double** Cfield __attribute__ ((aligned (L1ALIGN)));
    typedef double** Vfield __attribute__ ((aligned (L1ALIGN)));
    typedef double** V_of_fields __attribute__ ((aligned (L1ALIGN)));
    #define MACH_EPS DBL_EPSILON
#endif
/** The following macro VAS(Vfield) is conditional on the chosen contiguous memory alignment for vector fields. Hence, they can traversed using this macro just as a single large one-dimensional scalar field. **/
#define VAS(field) (*field)
/** The following macros make complex algebra a little bit more transparent.
  **/
#define CABS(z) (z[REAL]*z[REAL]+z[IMAG]*z[IMAG])
#define CMULS(x,y,s) { x[REAL] = y[REAL] * s; \
                       x[IMAG] = y[IMAG] * s; }
#define CMULC(x,y,z) { x[REAL] = y[REAL] * z[REAL] - y[IMAG] * z[IMAG]; \
                       x[IMAG] = y[IMAG] * z[REAL] + y[REAL] * z[IMAG]; }
#define CMULA(x,y,z) { x[REAL] = y[REAL] * z[REAL] + y[IMAG] * z[IMAG]; \
                       x[IMAG] = y[IMAG] * z[REAL] - y[REAL] * z[IMAG]; }

/******************************************************************************
  * making directions and loops more transparent throughout the code
  *****************************************************************************/
enum CMPLX {REAL, IMAG,NCMPLX};
#define NDIM 2
#define NDIR (NDIM * 2)
#define XUP 0
#define YUP 1
#define YDN 2
#define XDN 3
#define OPP_DIR(DIR) (NDIR-1 - DIR)
#define UP_DIR(DIR) ( DIR < NDIM ? DIR : OPP_DIR(DIR) )
#define DN_DIR(DIR) ( DIR < NDIM ? OPP_DIR(DIR) : DIR )
/* alternative scheme closer to what Ilaria had
#define XUP 0
#define XDN 1
#define YUP 2
#define YDN 3
#define OPP_DIR(DIR) (DIR+(1-2*(DIR%2)))
*/


#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

/******************************************************************************
  * list of global variables
  *****************************************************************************/

/* Physics parameters */
EXTERN Uint L; /* array for dimensions, will supersede L0, L1 */
EXTERN Uint VOL; /* volume */
EXTERN Uint NBOS; /* number of bosons can be a parameter */
EXTERN Uint NAUX; /* number of auxiliary fields depends on NBOS */
EXTERN Uint NFER; /* number of fermions could be a parameter */
EXTERN Scalar g; /* coupling strength */


EXTERN Uint RNG_seed;
EXTERN Uint save_each;

/* counter for warnings */
EXTERN Uint run_status;
/* counter for HMC acceptance */
EXTERN Uint total_accept;
/* counter for CG steps (manual accumulation through CG return value) */
EXTERN Uint total_iters;

/* files for I/O */
EXTERN char *Field_Savefile;
EXTERN char *State_Savefile;


EXTERN double time_algebra;
EXTERN double time_cg;
EXTERN double time_master;
EXTERN double time_random;

#endif /* _SETTING_H */

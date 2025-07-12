#ifndef _RNG_H
#define _RNG_H

#include <stddef.h>
#include<gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "setting.h"

enum RNG_STATE_FLAG {RNG_UNDEFINED, RNG_DEFINED};
enum RNG_IO_FLAG {RNG_IO_SUCCESS=0, RNG_IO_FAIL=GSL_EFAILED};

void rng_fill_array(Scalar *array, Uint n);
void rng_free();
Scalar rng_get();
int rng_read_state(FILE *fp);
int rng_save_state(FILE *fp);
void rng_set_seed(Uint seed);

#endif /* _RNG_H */

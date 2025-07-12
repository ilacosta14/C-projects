#include "rng.h"

#include <omp.h>
#include <gsl/gsl_rng.h>

static gsl_rng *r;

// Fills given array with n random numbers in the range [0,1).
void rng_fill_array(Scalar *array, Uint n) {
    if (!r) {
        const gsl_rng_type *T;

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }

    for (Uint i = 0; i < n; i++) {
        array[i] = gsl_rng_uniform(r);
    }
}

// Frees the random number generator.
void rng_free() {
    if (r) {
        gsl_rng_free(r);
        r = 0;
    }
}

// Returns a single random number.
Scalar rng_get() {
    if (!r) {
        const gsl_rng_type *T;

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }

    return gsl_rng_uniform(r);
}

// Writes the state of the random number generator.
int rng_read_state(FILE *fp) {
    return( gsl_rng_fread(fp, r) );
}

// Writes the state of the random number generator.
int rng_save_state(FILE *fp) {
    return( gsl_rng_fwrite(fp, r) );
}

// Set the seed for the random number generator.
void rng_set_seed(unsigned long int seed) {
    if (!r) {
        const gsl_rng_type *T;

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }

    gsl_rng_set(r, seed);
}

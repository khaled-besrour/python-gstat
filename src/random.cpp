#ifdef CPP_STANDALONE
    #include <R.h>
    #include <Rinternals.h>

    double r_uniform(void) {
        return(unif_rand());
    }

    double r_normal(void) {
        return(norm_rand());
    }
#else
    #include "cpp_utils/random_generator.hpp"
    extern "C" double r_uniform(void) {
        return(random_generator.unif_rand());
    }

    extern "C" double r_normal(void) {
        return(random_generator.norm_rand());
    }
#endif
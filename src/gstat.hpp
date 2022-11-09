//
//  s_python.h
//  kriging
//
//  Created by khaled besrour on 20/10/2022.
//

#ifndef gstat_hpp
#define gstat_hpp

#ifdef CPP_STANDALONE
#include "cpp_utils/random_generator.hpp"
#endif

#include <stdio.h>
#include <string>

#ifdef __cplusplus
extern "C"
{
#endif
#include "defs.h"
#include "data.h"
#include "select.h"
#include "utils.h"
#include "userio.h"
#include "vario.h"
#include "fit.h"
#include "sem.h"
#include "glvars.h"
#include "debug.h"
#include "mapio.h"
#include "msim.h"
#include "getest.h"
#ifdef __cplusplus
}
#endif

std::vector<std::string> gstat_get_variogram_models(int do_long);
int gstat_init(int s_debug_level, int random_seed);
void gstat_new_data(std::vector<double> & sy, std::vector<double> & slocs,
                    std::vector<double> & sX, long has_intercept,
                    std::vector<double> & beta, int nmax, int nmin,
                    double maxdist, int force, int vfn, std::vector<double> & sw,
                    std::vector<double> & grid, int degree, int is_projected,
                    int vdist, double lambda, int omax);

void gstat_new_dummy_data(int loc_dim, long has_intercept,
                          std::vector<double> & beta,
                          int nmax, int nmin, double maxdist, int vfn, int is_projected,
                          int vdist);

static void gstat_set_block(long i, std::vector<double> & block, std::vector<int> & block_cols, DPOINT *current);

static DATA_GRIDMAP *gstat_S_fillgrid(std::vector<double> & gridparams);

void gstat_load_variogram(std::vector<long> & s_ids, std::vector<std::string> & s_model, std::vector<double> & s_sills,
                          std::vector<double> & s_ranges, std::vector<double> & s_kappas,
                          std::vector<double> & s_anis_all, std::vector<double> & s_table, std::vector<double> & s_max_val);

std::vector<std::vector<double> > get_covariance_list(std::vector<long> & ids , int covariance, std::vector<int> & dist_list);

void gstat_load_ev(std::vector<unsigned long> & np, std::vector<double> & dist, std::vector<double> & gamma);

std::vector<std::vector<double> > gstat_variogram_values(std::vector<long> & ids, std::vector<double> & pars, int covariance, std::vector<int> & dist_values);


void gstat_predict(std::vector<double> & retvector, std::vector<int> & retvector_dim, long sn,  std::vector<double> & slocs,
                   std::vector<double> & sX, std::vector<int> & block_cols,
                   std::vector<double> & block, std::vector<double> & weights,
                   int nsim, long blue);

std::vector<std::vector<double> > gstat_variogram(std::vector<long> & s_ids, std::vector<double> & cutoff,  std::vector<double> & width,  std::vector<double> & direction,
                                                  int cressie, std::vector<double> & dX, std::vector<double> & boundaries, std::vector<double> & grid, int cov,
                                                  int pseudo) ;
std::vector<std::vector<double> > gstat_fit_variogram(int fit, std::vector<int> & fit_sill, std::vector<int> & fit_range);

std::string gstat_set_method(std::string to);
int gstat_set_merge(int id1, int col1, int id2, int col2);
int gstat_exit(int x);


#endif /* gstat_hpp */

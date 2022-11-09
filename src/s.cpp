
/*
 * all functions exposed to R
 */
#include <R.h>
#include <Rinternals.h>
/* #include <Rdefines.h> */

#include "gstat.hpp"
#include "s.hpp"
#include <vector>

#if defined(__cplusplus)
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
#if defined(__cplusplus)
}
#endif

static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current);
int do_print_progress = 0;
#define NAME_SIZE 20 /* buffer size for name */

extern unsigned int n_pred_locs; /* msim.c */

SEXP gstat_init(SEXP s_debug_level) {
	int level = INTEGER(s_debug_level)[0];
    int seed = 0; // Will take R random inside gstat_init

	gstat_init(level, seed);
    return s_debug_level;
}

SEXP gstat_exit(SEXP x) {
    int a = INTEGER(x)[0];
    gstat_exit(a);
	return(x);
}

void copy_vector(SEXP & output_vect, std::vector<double> & input_vect) {
    for (int j = 0 ; j < input_vect.size(); j++)
        REAL(output_vect)[j] = input_vect[j];
}

void copy_vector(SEXP & output_vect, std::vector<int> & input_vect) {
    for (int j = 0 ; j < input_vect.size(); j++)
        INTEGER(output_vect)[j] = input_vect[j];
}

SEXP gstat_new_data(SEXP sy, SEXP slocs, SEXP sX, SEXP has_intercept, 
			SEXP beta, SEXP nmax, SEXP nmin, SEXP maxdist, SEXP force,
			SEXP vfn, SEXP sw, SEXP grid, SEXP degree, SEXP is_projected,
			SEXP vdist, SEXP lambda, SEXP omax) {
    
    std::vector<double> cpp_sy  {REAL(sy), REAL(sy) + LENGTH(sy)};
    std::vector<double> cpp_slocs  {REAL(slocs), REAL(slocs) + LENGTH(slocs)};
    std::vector<double> cpp_sX  {REAL(sX), REAL(sX) + LENGTH(sX)};
    std::vector<double> cpp_beta  {REAL(beta), REAL(beta) + LENGTH(beta)};
    std::vector<double> cpp_sw  {REAL(sw), REAL(sw) + LENGTH(sw)};
    std::vector<double> cpp_grid  {REAL(grid), REAL(grid) + LENGTH(grid)};

    long cpp_has_intercept = INTEGER(has_intercept)[0];
    int cpp_nmax = INTEGER(nmax)[0];
    int cpp_nmin = INTEGER(nmin)[0];
    int cpp_force = INTEGER(force)[0];
    int cpp_vfn = INTEGER(vfn)[0];
    int cpp_degree = INTEGER(degree)[0];
    int cpp_is_projected = INTEGER(is_projected)[0];
    int cpp_vdist = INTEGER(vdist)[0];
    int cpp_omax = INTEGER(omax)[0];

    double cpp_lambda = REAL(lambda)[0];
    double cpp_maxdist = REAL(maxdist)[0];
    
    gstat_new_data(cpp_sy,
                   cpp_slocs,
                   cpp_sX,
                   cpp_has_intercept,
                   cpp_beta,
                   cpp_nmax,
                   cpp_nmin,
                   cpp_maxdist,
                   cpp_force,
                   cpp_vfn,
                   cpp_sw,
                   cpp_grid,
                   cpp_degree,
                   cpp_is_projected,
                   cpp_vdist,
                   cpp_lambda,
                   cpp_omax);
    
	return(sy);
}

SEXP gstat_new_dummy_data(SEXP loc_dim, SEXP has_intercept, SEXP beta, 
		SEXP nmax, SEXP nmin, SEXP maxdist, SEXP vfn, SEXP is_projected,
		SEXP vdist) {
    
    std::vector<double> cpp_beta  {REAL(beta), REAL(beta) + LENGTH(beta)};
    
    int cpp_loc_dim = INTEGER(loc_dim)[0];
    long cpp_has_intercept = INTEGER(has_intercept)[0];
    int cpp_nmax = INTEGER(nmax)[0];
    int cpp_nmin = INTEGER(nmin)[0];
    double cpp_maxdist = REAL(maxdist)[0];
    int cpp_vfn = INTEGER(vfn)[0];
    int cpp_is_projected = INTEGER(is_projected)[0];
    int cpp_vdist = INTEGER(vdist)[0];
    

    
    gstat_new_dummy_data(cpp_loc_dim,
                              cpp_has_intercept,
                              cpp_beta,
                              cpp_nmax,
                              cpp_nmin,
                              cpp_maxdist,
                              cpp_vfn,
                              cpp_is_projected,
                              cpp_vdist);

	return(loc_dim);
}

SEXP gstat_predict(SEXP sn, SEXP slocs, SEXP sX, SEXP block_cols, SEXP block, 
			SEXP weights, SEXP nsim, SEXP blue) {
    
    std::vector<double> cpp_retvector ;
    std::vector<int> cpp_retvector_dim;
    
    std::vector<double> cpp_slocs {REAL(slocs), REAL(slocs)+ LENGTH(slocs)};
    std::vector<double> cpp_sX {REAL(sX), REAL(sX)+ LENGTH(sX)};
    std::vector<double> cpp_block {REAL(block), REAL(block)+ LENGTH(block)};
    std::vector<double> cpp_weights {REAL(weights), REAL(weights)+ LENGTH(weights)};

    std::vector<int> cpp_block_cols {INTEGER(block_cols), INTEGER(block_cols) + LENGTH(block_cols)};

    int cpp_nsim = INTEGER(nsim)[0];
    long cpp_blue = INTEGER(blue)[0];
    long cpp_sn = INTEGER(sn)[0];

    gstat_predict(
                  cpp_retvector,
                  cpp_retvector_dim,
                  cpp_sn,
                  cpp_slocs,
                  cpp_sX,
                  cpp_block_cols,
                  cpp_block,
                  cpp_weights,
                  cpp_nsim,
                  cpp_blue);
    
	SEXP retvector;
	SEXP retvector_dim;
    SEXP ret;

    
	PROTECT(ret = allocVector(VECSXP, 1));
	PROTECT(retvector_dim = allocVector(REALSXP, 2));
    PROTECT(retvector = allocVector(REALSXP, cpp_retvector.size()));
    

    copy_vector(retvector_dim,cpp_retvector_dim);
    copy_vector(retvector, cpp_retvector);

	/* SET_DIM(retvector, retvector_dim); */
	setAttrib(retvector, R_DimSymbol, retvector_dim);
	SET_VECTOR_ELT(ret, 0, retvector);
    
	UNPROTECT(3);
	return(ret);
}

static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current) {
    std::vector<double> cpp_block {REAL(block), REAL(block) + LENGTH(block)};
    std::vector<int>    cpp_block_cols {INTEGER(block_cols), INTEGER(block_cols) + LENGTH(block_cols)};
    
    gstat_set_block(i, cpp_block,cpp_block_cols, current);
}

SEXP gstat_variogram(SEXP s_ids, SEXP cutoff, SEXP width, SEXP direction, 
		SEXP cressie, SEXP dX, SEXP boundaries, SEXP grid, SEXP cov,
		SEXP pseudo) {
    
    std::vector<double> cpp_cutoff {REAL(cutoff), REAL(cutoff) + LENGTH(cutoff)};
    std::vector<double> cpp_width {REAL(width), REAL(width) + LENGTH(width)};
    std::vector<double> cpp_direction {REAL(direction), REAL(direction) + LENGTH(direction)};
    std::vector<double> cpp_dX {REAL(dX), REAL(dX) + LENGTH(dX)};
    std::vector<double> cpp_boundaries {REAL(boundaries), REAL(boundaries) + LENGTH(boundaries)};
    std::vector<double> cpp_grid {REAL(grid), REAL(grid) + LENGTH(grid)};
    
    int cpp_cov = INTEGER(cov)[0];
    int cpp_pseudo = INTEGER(pseudo)[0];
    int cpp_cressie = INTEGER(cressie)[0];

    std::vector<long> cpp_s_ids {INTEGER(s_ids), INTEGER(s_ids) + LENGTH(s_ids)};
    
    std::vector<std::vector<double> > res = gstat_variogram(cpp_s_ids,
                                                            cpp_cutoff,
                                                            cpp_width,
                                                            cpp_direction,
                                                            cpp_cressie,
                                                            cpp_dX,
                                                            cpp_boundaries,
                                                            cpp_grid,
                                                            cpp_cov,
                                                            cpp_pseudo);
    
    SEXP ret;
    SEXP r1; // sx OR np
    SEXP r2; // sy OR dist
    SEXP r3; // np OR gamma
    SEXP r4; // gamma OR ev_parameters

    PROTECT(ret = allocVector(VECSXP, 4));
    
    PROTECT(r1 = allocVector(REALSXP, res[0].size()));
    PROTECT(r2 = allocVector(REALSXP, res[1].size()));
    PROTECT(r3 = allocVector(REALSXP, res[2].size()));
    PROTECT(r4 = allocVector(REALSXP, res[3].size()));
    
    copy_vector(r1, res[0]);
    copy_vector(r2, res[1]);
    copy_vector(r3, res[2]);
    copy_vector(r4, res[3]);
    
    SET_VECTOR_ELT(ret, 0, r1);
    SET_VECTOR_ELT(ret, 1, r2);
    SET_VECTOR_ELT(ret, 2, r3);
    SET_VECTOR_ELT(ret, 3, r4);
    UNPROTECT(5);
        
	return(ret);
}

SEXP gstat_load_variogram(SEXP s_ids, SEXP s_model, SEXP s_sills, SEXP s_ranges, 
		SEXP s_kappas, SEXP s_anis_all, SEXP s_table, SEXP s_max_val) 
{
    
    std::vector<long> cpp_s_ids;

    std::vector<std::string> cpp_s_model;

    std::vector<double> cpp_s_sills {REAL(s_sills), REAL(s_sills) + LENGTH(s_sills)};
    std::vector<double> cpp_s_ranges {REAL(s_ranges), REAL(s_ranges) + LENGTH(s_ranges)};
    std::vector<double> cpp_s_kappas {REAL(s_kappas), REAL(s_kappas) + LENGTH(s_kappas)};
    std::vector<double> cpp_s_anis_all {REAL(s_anis_all), REAL(s_anis_all) + LENGTH(s_anis_all)};
    std::vector<double> cpp_s_table {REAL(s_table), REAL(s_table) + LENGTH(s_table)};
    std::vector<double> cpp_s_max_val {REAL(s_max_val), REAL(s_max_val) + LENGTH(s_max_val)};
    
    long n = LENGTH(s_sills);
    for (int i = 0; i < n; i++)
        cpp_s_model.push_back(CHAR(STRING_ELT(s_model, i)));
    
    gstat_load_variogram(
                cpp_s_ids,
                cpp_s_model,
                cpp_s_sills,
                cpp_s_ranges,
                cpp_s_kappas,
                cpp_s_anis_all,
                cpp_s_table,
                cpp_s_max_val
                         );

	return(s_model);
}

SEXP gstat_variogram_values(SEXP ids, SEXP pars, SEXP covariance, SEXP dist_values) {
    
    std::vector<long>  cpp_ids  {INTEGER(ids ), INTEGER(ids ) + LENGTH(ids )};
    std::vector<double>  cpp_pars {REAL(pars), REAL(pars) + LENGTH(pars)};
    std::vector<int>  cpp_dist_values {INTEGER(dist_values), INTEGER(dist_values) + LENGTH(dist_values)};

    int cpp_covariance = INTEGER(covariance)[0];
    
    std::vector<std::vector<double> > res = gstat_variogram_values(cpp_ids, cpp_pars, cpp_covariance, cpp_dist_values);
    

	SEXP dist;
	SEXP gamma;
	SEXP ret;

    PROTECT(dist = allocVector(REALSXP, res[0].size()));
    PROTECT(gamma = allocVector(REALSXP, res[1].size()));
	PROTECT(ret = allocVector(VECSXP, 2));
    
    copy_vector(dist, res[0]);
    copy_vector(gamma, res[1]);
    
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

// Added by Paul Hiemstra, 30-06-2008
SEXP get_covariance_list(SEXP ids, SEXP covariance, SEXP dist_list) {
    SEXP dist, gamma, ret;
    
    std::vector<long> cpp_ids {INTEGER(ids), INTEGER(ids) + LENGTH(ids)};
    std::vector<int> cpp_dist_list {INTEGER(dist_list), INTEGER(dist_list) + LENGTH(dist_list)};

    int cpp_covariance = INTEGER(covariance)[0];
    
    std::vector<std::vector<double> > res = get_covariance_list(cpp_ids , cpp_covariance, cpp_dist_list);
    
    
    PROTECT(dist = allocVector(REALSXP, res[0].size()));
    PROTECT(gamma = allocVector(REALSXP, res[1].size()));
    PROTECT(ret = allocVector(VECSXP, 2));
    
    copy_vector(dist, res[0]);
    copy_vector(gamma, res[1]);
    
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

SEXP gstat_get_variogram_models(SEXP dolong) {
	SEXP ret;
	int i, n = 0, do_long;
	
	for (i = 1; v_models[i].model != NOT_SP; i++)
		n++;

	do_long = INTEGER(dolong)[0];
	PROTECT(ret = allocVector(STRSXP, n));
	for (i = 1; v_models[i].model != NOT_SP; i++)
		SET_STRING_ELT(ret, i-1, 
				mkChar(do_long ? v_models[i].name_long : v_models[i].name));
	UNPROTECT(1);
	return(ret);
}

SEXP gstat_load_ev(SEXP np, SEXP dist, SEXP gamma) {
    
    std::vector<unsigned long> cpp_np {REAL(np), REAL(np) + LENGTH(np)};
    std::vector<double> cpp_dist {REAL(dist), REAL(dist) + LENGTH(dist)};
    std::vector<double> cpp_gamma {REAL(gamma), REAL(gamma) + LENGTH(gamma)};
    
    gstat_load_ev(cpp_np , cpp_dist, cpp_gamma);

	return(np);
}

SEXP gstat_fit_variogram(SEXP fit, SEXP fit_sill, SEXP fit_range) {
    
    int cpp_fit = INTEGER(fit)[0];
    std::vector<int> cpp_fit_sill {INTEGER(fit_sill), INTEGER(fit_sill) + LENGTH(fit_sill)};
    std::vector<int> cpp_fit_range {INTEGER(fit_range), INTEGER(fit_range) + LENGTH(fit_range)};
    
    std::vector<std::vector<double> > res = gstat_fit_variogram(cpp_fit, cpp_fit_sill, cpp_fit_range);
    
	SEXP ret;
	SEXP sills;
	SEXP ranges;
	SEXP SSErr;
	SEXP fit_is_singular;

	PROTECT(sills = allocVector(REALSXP, res[0].size()));
	PROTECT(ranges = allocVector(REALSXP, res[1].size()));
    PROTECT(fit_is_singular = allocVector(REALSXP, res[2].size()));
    PROTECT(SSErr = allocVector(REALSXP, res[3].size()));
    PROTECT(ret = allocVector(VECSXP, 4));
    
    copy_vector(sills, res[0]);
    copy_vector(ranges, res[1]);
    copy_vector(fit_is_singular, res[2]);
    copy_vector(SSErr, res[3]);
    
	SET_VECTOR_ELT(ret, 0, sills);
	SET_VECTOR_ELT(ret, 1, ranges);
    SET_VECTOR_ELT(ret, 2, fit_is_singular);
    SET_VECTOR_ELT(ret, 3, SSErr);

	UNPROTECT(5);
	return(ret);
}

SEXP gstat_debug_level(SEXP level) {
	debug_level = INTEGER(level)[0];
	return(level);
}

SEXP gstat_set_method(SEXP to) {
    std::string what = CHAR(STRING_ELT(to, 0));
    gstat_set_method(what);
	return(to);
}

SEXP gstat_set_set(SEXP arg, SEXP val) {
	const char *name;
	int i;
    enum TYPE_VAR { UNKNOWN, IS_INT, IS_UINT, IS_REAL, IS_STRING, IS_D_VECTOR, NO_ARG } ;
    enum LIMIT { NOLIMIT, GEZERO, GTZERO } ;
    
	typedef struct {
		const char *name;
		void *ptr;
		TYPE_VAR what;
        LIMIT limit;
	} GSTAT_EXPR;
	const GSTAT_EXPR set_options[] = {
	{ "alpha",          &gl_alpha,        IS_REAL, GEZERO  },
	{ "beta",           &gl_beta,         IS_REAL, GEZERO  },
	{ "blas",           &gl_blas,         IS_INT,  GEZERO  },
	{ "choleski",       &gl_choleski,     IS_INT,  GEZERO  },
	{ "co$incide",      &gl_coincide,     IS_INT,  GEZERO  },
	{ "Cr$essie",       &gl_cressie,      IS_INT,  GEZERO  },
	{ "cutoff",         &gl_cutoff,       IS_REAL, GTZERO  },
	{ "de$bug",         &debug_level,     IS_INT,  GEZERO  },
	{ "fit",            &gl_fit,          IS_INT,  GEZERO  },
	{ "fit_l$imit",     &gl_fit_limit,    IS_REAL, GTZERO  },
	{ "fr$action",      &gl_fraction,     IS_REAL, GTZERO  },
	/* { "display",        &gl_display,      IS_STRING, NOLIMIT }, */
	{ "gls$_residuals", &gl_gls_residuals, IS_INT, GEZERO  },
	{ "id$p",           &gl_idp,          IS_REAL, GEZERO  },
	{ "in$tervals",     &gl_n_intervals,  IS_INT,  GTZERO  },
	{ "it$er",          &gl_iter,         IS_INT,  GEZERO  },
	{ "lhs",            &gl_lhs,          IS_INT,  GEZERO  },
	{ "longlat",        &gl_longlat,      IS_INT, GEZERO },
	{ "sim_beta",  		&gl_sim_beta,     IS_INT, GEZERO },
	{ "n_uk",           &gl_n_uk,         IS_INT,  GEZERO  },
	{ "numbers",        &gl_numbers,      IS_INT,  GEZERO  },
	{ "nb$lockdiscr",   &gl_nblockdiscr,  IS_INT,  GTZERO  },
	{ "no$check",       &gl_nocheck,      IS_INT,  GEZERO  },
	{ "ns$im",          &gl_nsim,         IS_INT,  GTZERO  },
	{ "or$der",         &gl_order,        IS_INT,  GEZERO },
	{ "q$uantile",      &gl_quantile,     IS_REAL, GEZERO  },
	{ "rowwise",        &gl_rowwise,      IS_INT,  GEZERO  },
	{ "rp",             &gl_rp,           IS_INT,  GEZERO  },
	{ "see$d",          &gl_seed,         IS_INT,  GTZERO  },
	{ "useed",          &gl_seed,         IS_UINT,  GEZERO  },
	{ "spa$rse",        &gl_sparse,       IS_INT,  GEZERO  },
	{ "spi$ral",        &gl_spiral,       IS_INT,  GEZERO  },
	{ "spl$it",         &gl_split,        IS_INT,  GTZERO  },
	{ "sy$mmetric",     &gl_sym_ev,       IS_INT,  GEZERO  },
	{ "tol_h$or",       &gl_tol_hor,      IS_REAL, GEZERO  },
	{ "tol_v$er",       &gl_tol_ver,      IS_REAL, GEZERO  },
	{ "v$erbose",       &debug_level,     IS_INT,  GEZERO  },
	{ "w$idth",         &gl_iwidth,       IS_REAL, GEZERO  },
	{ "x$valid",        &gl_xvalid,       IS_INT,  GEZERO  },
	{ "zero_di$st",     &gl_zero_est,     IS_INT,  GEZERO  },
	{ "zero",           &gl_zero,         IS_REAL, GEZERO  },
	{ "zm$ap",          &gl_zmap,         IS_REAL, NOLIMIT },
	{ NULL, NULL, TYPE_VAR(0), LIMIT(0) }
	};
	name = CHAR(STRING_ELT(arg, 0));
	for (i = 0; set_options[i].name; i++)
		if (almost_equals(name, set_options[i].name))
			break; /* break out i-loop */
	if (set_options[i].name == NULL)
		ErrMsg(ER_SYNTAX, name);

	if (almost_equals((const char *)name, "nb$lockdiscr"))
		gl_gauss = 0; /* side effect */

	switch (set_options[i].what) {
		case IS_INT: 
			*((int *) set_options[i].ptr) = asInteger(val);
			/* Rprintf("int arg: %s val %d\n", name, asInteger(val)); */
			break;
		case IS_UINT: 
			*((unsigned int *) set_options[i].ptr) = (unsigned int) asInteger(val);
			/* Rprintf("uint arg: %s val %d\n", name, asInteger(val)); */
			break;
		case IS_REAL: 
			*((double *) set_options[i].ptr) = asReal(val);
			/* Rprintf("real arg: %s val %d\n", name, asReal(val)); */
			break; 
		case IS_STRING: 
			*((const char **) set_options[i].ptr) = CHAR(STRING_ELT(val, 0));
			break;
		default:
			ErrMsg(ER_SYNTAX, name);
			break;
	}
	return val;
}

SEXP gstat_set_merge(SEXP a, SEXP b, SEXP c, SEXP d) { /* merge a(b) with c(d); */
    int id1 = INTEGER(a)[0];
    int col1 = INTEGER(b)[0];
    int id2 = INTEGER(c)[0];
    int col2 = INTEGER(d)[0];
    gstat_set_merge(id1, col1, id2, col2);
	return(a);
}

extern "C" double r_uniform(void) {
	return(unif_rand());
}

extern "C" double r_normal(void) {
	return(norm_rand());
}

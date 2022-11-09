#ifndef CPP_STANDALONE
    #include <R.h>
    #include <Rinternals.h>
#else
    #include <limits.h>
    #define NA_REAL ' '
    #define NA_INT std::numeric_limits<int>::min()
#endif


#include "gstat.hpp"
#include "s.hpp"
#include <vector>

#ifdef CPP_STANDALONE
    int do_print_progress = 0;
#endif

extern unsigned int n_pred_locs; /* msim.c */

int gstat_init(int s_debug_level, int random_seed) {
    
    do_print_progress = 0;
    remove_all();
    init_global_variables();
    init_data_minmax();
    
    #ifndef CPP_STANDALONE
        GetRNGstate();
    #else
        random_generator.set_seed(random_seed);
    #endif
    
    debug_level = s_debug_level;
    if (debug_level < 0) {
        debug_level = -debug_level;
        do_print_progress = 1;
    }
    return(s_debug_level);
}

int gstat_exit(int x) {
    /* write seed back to R/S engine */
    remove_all();
    return(x);
}

// x position puis y position
void gstat_new_data(std::vector<double> & sy, std::vector<double> & slocs,
                    std::vector<double> & sX, long has_intercept,
                    std::vector<double> & beta, int nmax, int nmin,
                    double maxdist, int force, int vfn, std::vector<double> & sw,
                    std::vector<double> & grid, int degree, int is_projected,
                    int vdist, double lambda, int omax) {
    
    double *y, *locs, *X, *w = NULL;
    long i, j, id, n, dim, n_X, has_int;
    DPOINT current;
    DATA **d;
    char name[NAME_SIZE];
    
    n = sy.size();
    y = &sy[0];
    
    if (n == 0)
        ErrMsg(ER_IMPOSVAL, "no data read");
    
    long slocs_size = slocs.size();
    long sX_size = sX.size();
    
    if (slocs_size % n != 0)
        error("dimensions do not match: locations %d and data %ld",
              slocs_size, n);
    dim = slocs_size / n;
    if (dim <= 0)
        error("too few spatial dimensions: %ld", dim);
    if (dim > 3)
        error("too many spatial dimensions: %ld", dim);
    locs = &slocs[0];
    
    if (sw.size() == n)
        w = &sw[0];
    
    if (sX_size % n != 0)
        error("dimensions do not match: X %d and data %ld: missing values in data?",
              sX_size, n);
    n_X = sX_size / n;
    X = &sX[0];
    
    assert(n_X > 0);
    current.z = 0.0;
    current.bitfield = 0;
    
    id = get_n_vars();
    snprintf(name, NAME_SIZE, "var%ld", id);
    which_identifier(name);
    d = get_gstat_data();
    d[id]->id = id;
    
    d[id]->n_list = d[id]->n_max = 0;
    d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
    d[id]->x_coord = "x";
    d[id]->y_coord = "y";
    d[id]->z_coord = "z";
    d[id]->variable = "R data";
    d[id]->fname = "R data";
    d[id]->lambda = lambda;
    has_int = has_intercept;
    /* increase d[id]->n_X and set d[id]->colX[i]: */
    for (i = d[id]->n_X = 0; i < n_X; i++)
        data_add_X(d[id], i + (has_int ? 0 : 1));
    assert(d[id]->n_X == n_X);
    for (i = 0; i < beta.size(); i++) /* do nothing if beta is numeric(0) */
        d[id]->beta = push_d_vector(beta[i], d[id]->beta);
    if (nmax > 0) /* leave default (large) if < 0 */
        d[id]->sel_max = nmax;
    if (omax > 0) /* leave default (0) if <= 0 */
        d[id]->oct_max = omax;
    if (nmin > 0) /* leave default (0) if <= 0 */
        d[id]->sel_min = nmin;
    if (maxdist > 0.0)
        d[id]->sel_rad = maxdist;
    if (force > 0)
        d[id]->force = force;
    switch(vfn) {
        case 1: /* d[id]->variance_fn = v_identity; == leave NULL */ break;
        case 2: d[id]->variance_fn = v_mu; break;
        case 3: d[id]->variance_fn = v_bin; break;
        case 4: d[id]->variance_fn = v_mu2; break;
        case 5: d[id]->variance_fn = v_mu3; break;
        default: error("unknown variance function %d", vfn);
    }
    gl_longlat = (is_projected == 0);
    d[id]->mode = X_BIT_SET | V_BIT_SET;
    if (dim > 1)
        d[id]->mode = d[id]->mode | Y_BIT_SET;
    if (dim > 2)
        d[id]->mode = d[id]->mode | Z_BIT_SET;
    set_norm_fns(d[id]); /* so gstat can calculate distances */
    if (w != NULL)
        d[id]->colnvariance = 1; /* it is set */
    switch(grid.size()) {
        case 0: case 1: break; /* empty, i.e., numeric(0) */
        case 6: d[id]->grid = gstat_S_fillgrid(grid); break;
        default: error("length of grid topology %d unrecognized", (int) grid.size());
    }
    d[id]->polynomial_degree = degree;
    if (d[id]->polynomial_degree < 0 || d[id]->polynomial_degree > 3) {
        error("polynomial degree should be 0, 1, 2 or 3");
    }
    if (d[id]->polynomial_degree > 0) {
        /* we're doing polynomials through degree: */
        if (id > 0) {
            error("polynomial degree will only work for a single variable");
        } if (n_X > 1) {
            error("polynomial degree only works when no other predictors are given");
        }
        setup_polynomial_X(d[id]); /* standardized coordinate polynomials */
    }
    d[id]->vdist = vdist;
    assert(n_X <= d[id]->n_X);
    current.X = (double *) emalloc(d[id]->n_X * sizeof(double));
    
    SET_POINT(&current);
    current.u.stratum = 0;
    current.attr = current.x = current.y = current.z = 0.0;
    for (i = 0; i < n; i++) { /* loop over points */
        current.attr = y[i];
        current.x = locs[i];
        if (dim >= 2)
            current.y = locs[n + i];
        if (dim >= 3)
            current.z = locs[2 * n + i];
        /* track min/max coordinates, also for z, for the qtree bbox */
        if (i == 0) {
            d[id]->maxX = d[id]->minX = current.x;
            d[id]->maxY = d[id]->minY = current.y;
            d[id]->maxZ = d[id]->minZ = current.z;
        } else {
            d[id]->minX = MIN(d[id]->minX, current.x);
            d[id]->maxX = MAX(d[id]->maxX, current.x);
            d[id]->minY = MIN(d[id]->minY, current.y);
            d[id]->maxY = MAX(d[id]->maxY, current.y);
            d[id]->minZ = MIN(d[id]->minZ, current.z);
            d[id]->minZ = MIN(d[id]->minZ, current.z);
        }
        for (j = 0; j < n_X; j++)
            current.X[j] = X[j * n + i];
        if (w != NULL)
            current.variance = 1.0/(w[i]);
        push_point(d[id], &current);
    }
    check_global_variables();
    d[id]->n_original = d[id]->n_list;
    efree(current.X);
}



void gstat_new_dummy_data(int loc_dim, long has_intercept,
                          std::vector<double> & beta,
                          int nmax, int nmin, double maxdist, int vfn, int is_projected,
                          int vdist) {
    int i, id, dim, has_int;
    char name[NAME_SIZE];
    DATA **d = NULL;
    long beta_size = beta.size();

    dim = loc_dim;
    if (dim <= 0)
        error("dimension value impossible: %d", dim);
    if (dim > 3)
        error("too many dimensions: %d", dim);
    assert(beta_size > 0);

    id = get_n_vars();
    snprintf(name, NAME_SIZE, "var%d", id);
    which_identifier(name);
    d = get_gstat_data();
    d[id]->id = id;

    d[id]->n_list = d[id]->n_max = 0;
    d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
    d[id]->x_coord = "x";
    d[id]->y_coord = "y";
    d[id]->z_coord = "z";
    d[id]->variable = "R data";
    d[id]->fname = "R data";
    has_int = has_intercept;
    for (i = d[id]->n_X = 0; i < beta_size; i++)
        data_add_X(d[id], i + (has_int ? 0 : 1));
    assert(d[id]->n_X == beta_size);
    d[id]->dummy = 1;
    for (i = 0; i < beta_size; i++)
        d[id]->beta = push_d_vector(beta[i], d[id]->beta);
    if (nmax > 0) /* leave default (large) if < 0 */
        d[id]->sel_max = nmax;
/* I doubt whether using nmin for dummy data _ever_ can have a
 * meaning, but hey, let's add it anyway. */
    if (nmin > 0) /* leave default (0) if <= 0 */
        d[id]->sel_min = nmin;
    if (maxdist > 0.0)
        d[id]->sel_rad = maxdist;
    switch(vfn) {
        case 1: /* d[id]->variance_fn = v_identity; -> leave NULL */ break;
        case 2: d[id]->variance_fn = v_mu; break;
        case 3: d[id]->variance_fn = v_bin; break;
        case 4: d[id]->variance_fn = v_mu2; break;
        case 5: d[id]->variance_fn = v_mu3; break;
        default: error("unknown variance function %d", vfn);
    }
    gl_longlat = (is_projected == 0);
    d[id]->vdist = vdist;
    d[id]->mode = X_BIT_SET | V_BIT_SET;
    if (dim > 1)
        d[id]->mode = d[id]->mode | Y_BIT_SET;
    if (dim > 2)
        d[id]->mode = d[id]->mode | Z_BIT_SET;
    set_norm_fns(d[id]); /* so gstat can calculate distances */
    check_global_variables();
    d[id]->n_original = d[id]->n_list;
}


void gstat_predict(std::vector<double> & retvector, std::vector<int> & retvector_dim, long sn,  std::vector<double> & slocs,
                                  std::vector<double> & sX, std::vector<int> & block_cols,
                                  std::vector<double> & block, std::vector<double> & weights,
                                  int nsim, long blue) {
    
    double *locs, **est_all, *X;
    long i, j, k, n, nvars, nest, dim, n_X, ncols_block,
        nrows_block, pos;
    DPOINT current, *bp = NULL;
    DATA **d = NULL, *vd = NULL, *area = NULL;

    extern unsigned int n_pred_locs; /* predict.c, used in msim.c */
    float ***msim = NULL;

    nvars = get_n_vars();
    nest = nvars + (nvars * (nvars + 1))/2;
    n = sn;
    long slocs_size = slocs.size();
    long sx_size = sX.size();
    long block_size = block.size();
    
    if (n <= 0 || slocs_size == 0 || sx_size == 0)
        ErrMsg(ER_IMPOSVAL, "newdata empty or only NA's");
    if (slocs_size % n != 0)
        error("dimensions do not match: locations %d, nrows in X %ld",
              slocs_size, n);
    dim = slocs_size / n;
    if (dim > 3)
        error("too many spatial dimensions: %ld", dim);
    if (dim <= 0)
        error("too few spatial dimensions: %ld", dim);
    locs = &slocs[0];
    if (sx_size % n != 0)
        error("dimensions do not match: X %d and data %ld", sx_size, n);
    n_X = sx_size / n;

    current.attr = current.x = current.y = current.z = 0.0;
    current.bitfield = 0;
    /* assuming points ... */
    SET_POINT(&current);
    /* and then do the block thing: */
    if (block_cols.size() == 0) {
        bp = get_block_p();
        bp->x = bp->y = bp->z = 0.0; /* obsolete, I'd guess */
        if (block_size >= 1) {
            bp->x = block[0];
            SET_BLOCK(&current);
        }
        if (block_size >= 2)
            bp->y = block[1];
        if (block_size >= 3)
            bp->z = block[2];
        if (block_size > 3)
            pr_warning("block dimension can only be 3; using the first 3");
    } else if (block_cols.size() == 1) { /* if > 1, block contains multiple 2D blocks */
        ncols_block = block_cols[0];
        if (ncols_block < 1 || ncols_block > 3)
            ErrMsg(ER_IMPOSVAL, "block dimensions should be in [1..3]");
        nrows_block = block_size / ncols_block; /* nr of rows */
        if (nrows_block > 0) {
            area = create_data_area();
            area->colnvariance = 0;
            area->n_list = area->n_max = 0;
            area->id = ID_OF_AREA;
            area->mode = X_BIT_SET;
            if (ncols_block > 1)
                area->mode = area->mode & Y_BIT_SET;
            if (ncols_block > 2)
                area->mode = area->mode & Z_BIT_SET;
            for (i = 0; i < nrows_block; i++) {
                current.x = block[i];
                if (ncols_block > 1)
                    current.y = block[nrows_block + i];
                if (ncols_block > 2)
                    current.z = block[2 * nrows_block + i];
                if (weights.size() > 0) {
                    area->colnvariance = 1;
                    current.variance = weights[i];
                }
                push_point(area, &current);
            }
            SET_BLOCK(&current);
        }
        if (DEBUG_FORCE)
            print_data_list(area);
    }

    X = &sX[0];
    assert(n_X > 0);
    current.X = (double *) emalloc(n_X * sizeof(double));
    current.u.stratum = 0;
    d = get_gstat_data();
    est_all = (double **) emalloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
        est_all[i] = (double *) emalloc(nest * sizeof(double));
    /*
     * the following is to fake gstat's default method handling:
     * we got to suggest that we'll go through a list of prediction
     * locations, a la the gstat ``data(): ... ;'' command.
     */
    vd = get_dataval();
    vd->id = ID_OF_VALDATA;
    vd->mode = d[0]->mode;
    /* set min/max[XYZ] */
    vd->minY = vd->maxY = vd->minZ = vd->maxZ = 0.0;
    vd->minX = vd->maxX = locs[0];
    for (i = 1; i < n; i++) {
        vd->minX = MIN(vd->minX, locs[i]);
        vd->maxX = MAX(vd->maxX, locs[i]);
    }
    if (dim >= 2) {
        vd->minY = vd->maxY = locs[n];
        for (i = 1; i < n; i++) {
            vd->minY = MIN(vd->minY, locs[n + i]);
            vd->maxY = MAX(vd->maxY, locs[n + i]);
        }
    }
    if (dim >= 3) {
        vd->minZ = vd->maxZ = locs[2 * n];
        for (i = 1; i < n; i++) {
            vd->minZ = MIN(vd->minZ, locs[2 * n + i]);
            vd->maxZ = MAX(vd->maxZ, locs[2 * n + i]);
        }
    }

    /* fill, and standardize coordinate predictors from degree = x */
    for (i = 0; i < nvars; i++)
        setup_data_minmax(d[i]);
    setup_data_minmax(vd);
    for (i = 0; i < nvars; i++)
        calc_polynomials(d[i]);
    /* calc_polynomials(vd); */ /* still no data in fake vd */

    vd->polynomial_degree = d[0]->polynomial_degree;
    if (vd->polynomial_degree > 0) {
        setup_polynomial_X(vd); /* standardized coordinate polynomials */
        current.X = (double *) erealloc(current.X, vd->n_X * sizeof(double));
    }

    /* so far for the faking; now let's see what gstat makes out of this: */
    if (nsim == 0) {
        if (blue == 0) { /* FALSE */
            if (get_method() == NSP) /* choose default */
                set_method(get_default_method());
        } else
            set_method(LSLM);
    } else {
        if (nsim < 0) {
            gl_nsim = -(nsim);
            set_method(ISI);
        } else {
            gl_nsim = nsim;
            set_method(GSI);
        }
        n_pred_locs = n;
        if (gl_nsim > 1)
            init_simulations(d);
        if (get_n_beta_set() != get_n_vars())
            setup_beta(d, get_n_vars(), gl_nsim);
    }
    set_mode();  /* simple, stratified, multivariable? */
    check_global_variables(); /* it's there, better do it now */
    if (debug_level)
        Rprintf("[%s]\n", method_string(get_method()));

    #ifndef CPP_STANDALONE
        #ifdef WIN32
            R_FlushConsole();
            R_ProcessEvents();
        #endif
    #endif
    
    for (i = 0; i < n; i++) {
        print_progress(i, n);
        if (block_cols.size() > 1)
            gstat_set_block(i, block, block_cols, &current);
        current.x = locs[i];
        if (dim >= 2)
            current.y = locs[n + i];
        if (dim >= 3)
            current.z = locs[2 * n + i];
        for (j = 0; j < n_X; j++)
            current.X[j] = X[j * n + i];
        /* transform standardized coordinate polynomial here */
        if (vd->polynomial_degree)
            calc_polynomial_point(vd, &current);
        for (j = 0; j < get_n_vars(); j++)
            select_at(d[j], &current);
        get_est(d, get_method(), &current, est_all[i]);
        #ifndef CPP_STANDALONE
            #ifdef WIN32
                R_FlushConsole();
                R_ProcessEvents();
            #endif
            R_CheckUserInterrupt();
        #endif
    }
    
    print_progress(100, 100);
    retvector_dim.resize(2);
    retvector_dim[0] = n; /* nrows */

    if (gl_nsim > 1) {
        retvector.resize(gl_nsim * nvars * n);
        msim = get_msim();
        for (i = pos = 0; i < nvars; i++)
            for (j = 0; j < gl_nsim; j++)
                for (k = 0; k < n; k++) {
                    if (is_mv_float(&(msim[i][k][j])))
                        retvector[pos++] = NA_REAL;
                    else
                        retvector[pos++] = msim[i][k][j];
                }
        retvector_dim[1] = nvars * gl_nsim; /* ncols */
    } else {
        retvector.resize(n * nest);
        for (j = pos = 0; j < nest; j++) {
            for (i = 0; i < n; i++) {
                if (is_mv_double(&(est_all[i][j])))
                    retvector[pos] =  NA_REAL;
                else
                    retvector[pos] = est_all[i][j];
                pos++;
            }
        }
        retvector_dim[1] = nest; /* ncols */
    }
    if (gl_nsim > 0)
        free_simulations();
    
    for (i = 0; i < n; i++)
        efree(est_all[i]);
    efree(est_all);
    efree(current.X);
}


std::vector<std::vector<double> > gstat_variogram(std::vector<long> & s_ids, std::vector<double> & cutoff,  std::vector<double> & width,  std::vector<double> & direction,
        int cressie, std::vector<double> & dX, std::vector<double> & boundaries, std::vector<double> & grid, int cov,
        int pseudo) {
    /* SEXP y; */
    long i, id1, id2, nest;
    VARIOGRAM *vgm;
    DATA **d;
    
    std::vector<std::vector<double> > ret;
    
    GRIDMAP *m;
    unsigned int row, col, n;

    id1 = s_ids[0];
    if (s_ids.size() > 1)
        id2 = s_ids[1];
    else
        id2 = id1;
    vgm = get_vgm(LTI(id1,id2));
    vgm->id = LTI(id1,id2);
    vgm->id1 = id1;
    vgm->id2 = id2;
    if (cov == 0)
        vgm->ev->evt = (id1 == id2 ? SEMIVARIOGRAM : CROSSVARIOGRAM);
    else if (cov == 1)
        vgm->ev->evt = (id1 == id2 ? COVARIOGRAM : CROSSCOVARIOGRAM);
    else {
        if (id1 != id2)
            ErrMsg(ER_IMPOSVAL,
            "cannot compute pairwise relative cross semivariogram");
        if (cov == 2)
            vgm->ev->evt = PRSEMIVARIOGRAM;
    }
    /* vgm->ev->is_asym = INTEGER(asym)[0]; */
    vgm->ev->pseudo = pseudo;
    vgm->ev->recalc = 1;
    if (cutoff.size() > 0)
        gl_cutoff = cutoff[0];
    if (width.size() > 0)
        gl_iwidth = width[0];
    gl_alpha = direction[0];
    gl_beta = direction[1];
    gl_tol_hor = direction[2];
    gl_tol_ver = direction[3];
    gl_cressie = cressie;
    if (dX.size() > 0) {
        d = get_gstat_data();
        d[id1]->dX = dX[0];
        d[id2]->dX = dX[0];
    }
    for (i = 0; i < boundaries.size(); i++) /* does nothing if LENGTH is 0 */
        push_bound(boundaries[i]);
    switch (grid.size()) {
        case 0: case 1: break;
        case 6: vgm->ev->S_grid = gstat_S_fillgrid(grid); break;
        default: error("unrecognized grid length in gstat_variogram");
            break;
    }

    calc_variogram(vgm, NULL);

    if (vgm->ev->S_grid != NULL) {
        m =(gridmap *) vgm->ev->map;
        n = m->rows * m->cols;
        std::vector<double> np(n);
        std::vector<double> gamma(n);
        std::vector<double> sx(n);
        std::vector<double> sy(n);
        
        for (row = i = 0; row < m->rows; row++) {
            for (col = 0; col < m->cols; col++) {
                map_rowcol2xy(m, row, col, &(sx[i]),
                                &(sy[i]));
                np[i] = vgm->ev->nh[i];
                if (vgm->ev->nh[i] > 0)
                    gamma[i] = vgm->ev->gamma[i];
                else
                    gamma[i] = NA_REAL;
                i++;
            }
        }
        ret.push_back(sx);
        ret.push_back(sy);
        ret.push_back(np);
        ret.push_back(gamma);
        free_data_gridmap((DATA_GRIDMAP *) vgm->ev->S_grid);
    } else {
        if (vgm->ev->cloud)
            nest = vgm->ev->n_est;
        else {
            if (vgm->ev->zero == ZERO_SPECIAL)
                nest = vgm->ev->n_est;
            else
                nest = vgm->ev->n_est - 1;
        }
        if (nest <= 0) {
            return(ret);
        }
        std::vector<double> np(nest);
        std::vector<double> dist(nest);
        std::vector<double> gamma(nest);
        std::vector<double> ev_parameters(4);
        
        ev_parameters[0] = vgm->ev->cutoff;
        ev_parameters[1] = vgm->ev->iwidth;
        ev_parameters[2] = vgm->ev->pseudo;
        ev_parameters[3] = vgm->ev->is_asym;
        
        for (i = 0; i < nest; i++) {
            np[i] = vgm->ev->nh[i];
            dist[i] = vgm->ev->dist[i];
            gamma[i] = vgm->ev->gamma[i];
        }
        ret.push_back(np);
        ret.push_back(dist);
        ret.push_back(gamma);
        ret.push_back(ev_parameters);

    }
    return(ret);
}


void gstat_load_variogram(std::vector<long> & s_ids, std::vector<std::string> & s_model, std::vector<double> & s_sills,
                          std::vector<double> & s_ranges, std::vector<double> & s_kappas,
                          std::vector<double> & s_anis_all, std::vector<double> & s_table, std::vector<double> & s_max_val)
{
    VARIOGRAM *vgm;
    long i, n, id1, id2, max_id;
    double anis[5] = {0.0, 0.0, 0.0, 1.0, 1.0}, rpars[2], *sills, *ranges,
        *kappas, *anis_all;
    const char *model;

    sills = &s_sills[0];
    ranges = &s_ranges[0];
    kappas = &s_kappas[0];
    anis_all = &s_anis_all[0];

    id1 = s_ids[0];
    id2 = s_ids[1];
    max_id = MAX(id1, id2);

    if (get_n_vars() == 0)
        which_identifier("xx"); /* at least "load" one dummy var */
    if (max_id >= get_n_vars())
        ErrMsg(ER_IMPOSVAL,
            "gstat_load_variogram has been called with max_id >= n_vars");

    vgm = get_vgm(LTI(id1,id2));
    assert(vgm != NULL);

    vgm->id = LTI(id1,id2);
    vgm->id1 = id1;
    vgm->id2 = id2;
    vgm->n_models = vgm->n_fit = 0;

    n = s_sills.size();
    for (i = 0; i < n; i++) { /* loop over sub models */
        model = s_model[i].c_str();
        anis[0] = anis_all[0 * n + i];
        anis[1] = anis_all[1 * n + i];
        anis[2] = anis_all[2 * n + i];
        anis[3] = anis_all[3 * n + i];
        anis[4] = anis_all[4 * n + i];
        rpars[0] = ranges[i];
        rpars[1] = kappas[i];
        if (s_table.size() > 0)
            push_to_v_table(vgm, rpars[0],
                    s_table.size(), &s_table[0],
                    (anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis);
        else
            push_to_v(vgm, model, sills[i], rpars, 2,
                (anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis, 1, 1);
    }
    update_variogram(vgm);
    if (s_max_val[0] > 0.0 || s_max_val[1] > 0.0 || s_max_val[2] > 0.0)
        vgm->max_val = get_semivariance(vgm,
                s_max_val[0], s_max_val[1], s_max_val[2]);
    if (DEBUG_DUMP)
        logprint_variogram(vgm, 1);
}


std::vector<std::vector<double> > gstat_variogram_values(std::vector<long> & ids, std::vector<double> & pars, int covariance, std::vector<int> & dist_values) {
    double from, to, n, d, x = 1.0, y = 0.0, z = 0.0;
    int i, id1, id2, cov = 0, ndist = 0;
    VARIOGRAM *vgm;
    
    std::vector<double> dist;
    std::vector<double> gamma;
    std::vector<std::vector<double> > ret;
    
    long pars_size = pars.size();

    if (pars_size != 3 && pars_size != 6)
        error("supply three or six distance parameters");
    from = pars[0];
    to = pars[1];
    n = pars[2];
    ndist = dist_values.size();
    cov = covariance;
    if (pars_size == 6) {
        x = pars[3];
        y = pars[4];
        z = pars[5];
    }

    id1 = ids[0];
    id2 = ids[1];
    vgm = get_vgm(LTI(id1,id2));

    if (ndist > 0) {
        dist.resize(ndist);
        gamma.resize(ndist);
        for (i = 0; i < ndist; i++) {
            d = dist_values[i];
            dist[i] = d;
            gamma[i] = (cov ?
                get_covariance(vgm, d * x, d * y, d * z) :
                get_semivariance(vgm, d * x, d * y, d * z));
        }
    } else {
        dist.resize(n);
        gamma.resize(n);
        for (i = 0; i < n; i++) {
            d = from;
            if (i > 0) /* implies n > 1 */
                d += (i/(n-1))*(to-from);
            dist[i] = d;
            gamma[i] = (cov ?
                get_covariance(vgm, d * x, d * y, d * z) :
                get_semivariance(vgm, d * x, d * y, d * z));
        }
    }
    ret.push_back(dist);
    ret.push_back(gamma);
    return(ret);
}

// Added by Paul Hiemstra, 30-06-2008 edited by Khaled Besrour in 20-10-2022
std::vector<std::vector<double> > get_covariance_list(std::vector<long> & ids , int covariance, std::vector<int> & dist_list) {
    double d, x = 1.0, y = 0.0, z = 0.0;
    int i, id1, id2, cov = 0;
    VARIOGRAM *vgm;
    std::vector<double> dist;
    std::vector<double> gamma;
    std::vector<std::vector<double> > ret;
    int length_list = dist_list.size();

    cov = covariance;

    id1 = ids[0];
    id2 = ids[1];
    vgm = get_vgm(LTI(id1,id2));
    
    dist.resize(length_list);
    gamma.resize(length_list);

    for (i = 0; i < length_list; i++) {
        d = dist_list[i];
        dist[i] = d;
        gamma[i] = (cov ?
            get_covariance(vgm, d * x, d * y, d * z) :
            get_semivariance(vgm, d * x, d * y, d * z));
    }
    
    ret.push_back(dist);
    ret.push_back(gamma);
    return(ret);
}

std::vector<std::string> gstat_get_variogram_models(int do_long) {
    std::vector<std::string> ret;
    int i, n = 0;
    
    for (i = 1; v_models[i].model != NOT_SP; i++)
        n++;
    
    for (i = 1; v_models[i].model != NOT_SP; i++)
        ret.push_back(do_long ? v_models[i].name_long : v_models[i].name);

    return(ret);
}

void gstat_load_ev(std::vector<unsigned long> & np, std::vector<double> & dist, std::vector<double> & gamma) {

    int i, cloud = 1;
    VARIOGRAM *vgm;

    which_identifier("xx");
    /*
     * vgm = get_vgm(LTI(INTEGER(id)[0], INTEGER(id)[1]));
     * */
    vgm = get_vgm(LTI(0, 0));
    if (vgm->ev == NULL)
        vgm->ev = init_ev();
    vgm->ev->evt = SEMIVARIOGRAM;
    vgm->ev->n_est = np.size();
    vgm->ev->n_max = np.size();
    vgm->ev->gamma = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
    vgm->ev->dist = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
    vgm->ev->nh = (unsigned long *) emalloc (sizeof(long) * vgm->ev->n_max);
    for (i = 0; i < vgm->ev->n_est; i++) {
        vgm->ev->nh[i] = np[i];
        vgm->ev->dist[i] = dist[i];
        vgm->ev->gamma[i] = gamma[i];
        if (cloud && vgm->ev->nh[i] > 1)
            cloud = 0;
    }
    vgm->ev->cloud = cloud;
    if (DEBUG_VGMFIT)
        fprint_sample_vgm(vgm->ev);
}


std::vector<std::vector<double> > gstat_fit_variogram(int fit, std::vector<int> & fit_sill, std::vector<int> & fit_range) {
    int i;
    VARIOGRAM *vgm;
    
    std::vector<std::vector<double> > ret;
    
    std::vector<double> sills;
    std::vector<double> ranges;
    std::vector<double> SSErr;
    std::vector<double> fit_is_singular;
    
    vgm = get_vgm(LTI(0, 0));
    vgm->ev->fit = FIT_TYPE(fit);
    for (i = 0; i < vgm->n_models; i++) {
        vgm->part[i].fit_sill = fit_sill[i];
        vgm->part[i].fit_range = fit_range[i];
    }
    update_variogram(vgm);
    if (DEBUG_VGMFIT)
        logprint_variogram(vgm, 1);
    fit_variogram(vgm);
    if (DEBUG_VGMFIT)
        logprint_variogram(vgm, 1);
    
    sills.resize(vgm->n_models);
    ranges.resize(vgm->n_models);
    
    for (i = 0; i < vgm->n_models; i++) {
        sills[i] = vgm->part[i].sill;
        ranges[i] = vgm->part[i].range[0];
    }

    
    fit_is_singular.push_back(vgm->fit_is_singular);
    SSErr.push_back(vgm->SSErr);
    
    ret.push_back(sills);
    ret.push_back(ranges);
    ret.push_back(fit_is_singular);
    ret.push_back(SSErr);
    
    return(ret);
}

//-------------------------------------------------------
static void gstat_set_block(long i, std::vector<double> & block, std::vector<int> & block_cols, DPOINT *current) {
    DATA *area;
    VARIOGRAM *v;
    long nrows_block, start, end, j;
    
    if (i >= block_cols.size() || i < 0)
        ErrMsg(ER_IMPOSVAL, "block_cols length less than nr of prediction locations");
    nrows_block = block.size() / 2; /* nr of rows */
    start = block_cols[i];
    if (i == block_cols.size() - 1)
        end = nrows_block;
    else
        end = block_cols[i+1] - 1;
    area = get_data_area();
    if (area != NULL)
        free_data(area);
    area = create_data_area();
    area->n_list = area->n_max = 0;
    area->id = ID_OF_AREA;
    area->mode = X_BIT_SET & Y_BIT_SET;
    for (j = start - 1; j < end; j++) {
        current->x = block[j];
        current->y = block[nrows_block + j];
        push_point(area, current);
    }
    SET_BLOCK(current);
    if (DEBUG_FORCE)
        print_data_list(area);
    for (j = 0; j < get_n_vgms(); j++) {
        v = get_vgm(j);
        if (v != NULL)
            v->block_semivariance_set = v->block_covariance_set = 0; /* don't store under these circumstances! */
    }
    return;
}

#ifdef CPP_STANDALONE
    extern "C" double r_uniform(void) {
        return(random_generator.unif_rand());
    }

    extern "C" double r_normal(void) {
        return(random_generator.norm_rand());
    }
#endif

int gstat_debug_level(int level) {
    debug_level = level;
    return(level);
}

std::string gstat_set_method(std::string to) {
    const char *what;

    what = to.c_str();
    for (int id = 1; methods[id].name != NULL; id++) {
        if (almost_equals(what, methods[id].name)) {
            set_method(methods[id].m);
            break; /* id-loop */
        }
    }
    return(to);
}

int gstat_set_merge(int id1, int col1, int id2, int col2) { /* merge a(b) with c(d); */
    DATA **dpp;
    int id;

    if (id1 >= get_n_vars() || id2 >= get_n_vars() || id1 < 0 || id2 < 0)
        ErrMsg(ER_IMPOSVAL, "id values out of range");
    if (id1 < id2) { /* swap id and col */
        id = id1; id1 = id2; id2 = id;
        id = col1; col1 = col2; col2 = id;
    }
    dpp = get_gstat_data();
    if (push_to_merge_table(dpp[id1], id2, col1, col2))
        ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
    return(id1);
}

static DATA_GRIDMAP *gstat_S_fillgrid(std::vector<double> & gridparams) {
    double x_ul, y_ul, cellsizex, cellsizey;
    unsigned int rows, cols;
    
    cellsizex = gridparams[2];
    cellsizey = gridparams[3];
    rows = (unsigned int) gridparams[5];
    cols = (unsigned int) gridparams[4];
    x_ul = gridparams[0] - 0.5 * cellsizex;
    y_ul = gridparams[1] + (rows - 0.5) * cellsizey;
    
    return gsetup_gridmap(x_ul, y_ul, cellsizex, cellsizey, rows, cols);
}

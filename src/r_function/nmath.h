/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2012 Ross Ihaka and the R Core team.
 *  Copyright (C) 2002-3    The R Foundation
 *
 *  Copied FROM https://github.com/atks/Rmath/tree/master
 *  Edited by khaled besrour on 20/10/2022.
 */

#ifndef nmath_h
#define nmath_h
#define M_SQRT_2dPI     0.797884560802865355879892119869        /* sqrt(2/pi) */
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2011  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

/* Required by C99 but might be slow */
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

#include <math.h>
#include <float.h> /* DBL_MIN etc */

/* Mathlib standalone */

#include <stdio.h>
#include <stdlib.h> /* for exit */
#define MATHLIB_ERROR(fmt,x)    { printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)        printf(fmt,x)
#define MATHLIB_WARNING2(fmt,x,x2)    printf(fmt,x,x2)
#define MATHLIB_WARNING3(fmt,x,x2,x3)    printf(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)

#define ISNAN(x) (isnan(x)!=0)
#define R_FINITE(x)    R_finite(x)
int R_finite(double);

#define ML_POSINF    (1.0 / 0.0)
#define ML_NEGINF    ((-1.0) / 0.0)
#define ML_NAN        (0.0 / 0.0)

#define _(String) String

#define ML_VALID(x)    (!ISNAN(x))

#define ME_NONE        0
/*    no error */
#define ME_DOMAIN    1
/*    argument out of domain */
#define ME_RANGE    2
/*    value out of range */
#define ME_NOCONV    4
/*    process did not converge */
#define ME_PRECISION    8
/*    does not have "full" precision */
#define ME_UNDERFLOW    16
/*    and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

/* For a long time prior to R 2.3.0 ML_ERROR did nothing.
   We don't report ME_DOMAIN errors as the callers collect ML_NANs into
   a single warning.
 */
#define ML_ERROR(x, s) { \
   if(x > ME_DOMAIN) { \
       char *msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
       msg = _("argument out of domain in '%s'\n");    \
       break; \
       case ME_RANGE: \
       msg = _("value out of range in '%s'\n");    \
       break; \
       case ME_NOCONV: \
       msg = _("convergence failed in '%s'\n");    \
       break; \
       case ME_PRECISION: \
       msg = _("full precision may not have been achieved in '%s'\n"); \
       break; \
       case ME_UNDERFLOW: \
       msg = _("underflow occurred in '%s'\n");    \
       break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0           Rf_bd0
#define chebyshev_eval    Rf_chebyshev_eval
#define chebyshev_init    Rf_chebyshev_init
#define gammalims    Rf_gammalims
#define lfastchoose    Rf_lfastchoose
#define lgammacor    Rf_lgammacor
#define stirlerr           Rf_stirlerr
/* in Rmath.h
#define gamma_cody      Rf_gamma_cody
*/

#endif /* MATHLIB_PRIVATE_H */

#endif /* nmath_h */

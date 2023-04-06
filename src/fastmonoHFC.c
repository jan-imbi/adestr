/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998--2021  The R Core Team
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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
 *  https://www.R-project.org/Licenses/
 */

/*	Spline Interpolation
 *	--------------------
 *	This modifies the spline evaluation algorithm
 *	from the stats package to work with
 *	monoH.FC splines.
 */

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "modreg.h"

// #include <R_ext/Applic.h>
// #include "Defn.h" // for R_xlen_t, UNIMPLEMENTED_TYPE


SEXP getListElement(SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < Rf_length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
    return elmt;
}

static void
spline_eval(int method, R_xlen_t nu, double *u, double *v,
            R_xlen_t n, double *x, double *y, double *b, double *c, double *d,
            double ml, double mr)
{
  /* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
   * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
   */
  const R_xlen_t n_1 = n - 1;
  R_xlen_t i, l;
  double dx;

  for(l = 0; l < nu; l++) v[l] = u[l];

  for(l = 0, i = 0; l < nu; l++) {
    double ul = v[l];
    if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
      /* reset i  such that  x[i] <= ul <= x[i+1] : */
      i = 0;
      R_xlen_t j = n;
      do {
        R_xlen_t k = (i+j) / 2;
        if(ul < x[k]) j = k; else i = k;
      } while(j > i+1);
    }
    dx = ul - x[i];

    if (method == 1 && ul < x[0]) {
      v[l] = y[i] + dx*ml;
    } else if (method == 1 && ul > x[n_1]) {
      v[l] = y[i] + dx*mr;
    } else {
      v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]));
    }
  }
}

// TODO: move to ../../../main/coerce.c
static R_xlen_t asXlen(SEXP x) {
  if (Rf_isVectorAtomic(x) && XLENGTH(x) >= 1) {
    switch (TYPEOF(x)) {
    case INTSXP:
      return (R_xlen_t) INTEGER(x)[0];
    case REALSXP:
      return (R_xlen_t) REAL(x)[0];
    }
  }
  return NA_INTEGER;
}

SEXP MySplineEval(SEXP xout, SEXP z)
{
  xout = PROTECT(Rf_coerceVector(xout, REALSXP));
  R_xlen_t nu = XLENGTH(xout), nx = asXlen(getListElement(z, "n"));
  SEXP yout = PROTECT(Rf_allocVector(REALSXP, nu));
  int method = Rf_asInteger(getListElement(z, "method"));
  double ml = Rf_asReal(getListElement(z, "ml"));
  double mr = Rf_asReal(getListElement(z, "mr"));
  SEXP x = getListElement(z, "x"), y = getListElement(z, "y"),
    b = getListElement(z, "b"), c = getListElement(z, "c"),
    d = getListElement(z, "d");
  spline_eval(method, nu, REAL(xout), REAL(yout),
              nx, REAL(x), REAL(y), REAL(b), REAL(c), REAL(d),
             ml, mr);
  UNPROTECT(2);
  return yout;
}

/* .Call calls */
extern SEXP MySplineEval(SEXP, SEXP);
extern SEXP MymonoFC_m(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"MySplineEval", (DL_FUNC) &MySplineEval, 2},
  {"MymonoFC_m", (DL_FUNC) &MymonoFC_m, 2},
  {NULL, NULL, 0}
};

void R_init_adestr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}



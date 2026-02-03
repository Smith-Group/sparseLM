#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Call routines */
extern SEXP sparselm(SEXP func, SEXP fjac, SEXP p, SEXP x, SEXP nconvars,
                     SEXP Jnnz, SEXP JtJnnz, SEXP itmax, SEXP opts, SEXP dif,
                     SEXP rho);
extern SEXP add_assign_col_inplace(SEXP x, SEXP i, SEXP j, SEXP value,
                                   SEXP validate);

static const R_CallMethodDef CallEntries[] = {
  {"C_sparselm", (DL_FUNC) &sparselm, 11},
  {"C_add_assign_col_inplace", (DL_FUNC) &add_assign_col_inplace, 5},
  {NULL, NULL, 0}
};

void R_init_sparseLM(DllInfo *dll) attribute_visible;

void R_init_sparseLM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

#include <R.h>
#include <Rdefines.h>
#include "splm.h"

typedef struct SparseLMData {

	SEXP func;
	SEXP fjac;
	SEXP pnames;
	SEXP rho;

} SparseLMData;

/* helper to match sparse matrix class */
static int sparselm_match_dg_class(SEXP obj)
{
    int matclass = -2;
    SEXP klass = PROTECT(Rf_getAttrib(obj, R_ClassSymbol));
    R_len_t n = Rf_length(klass);
    if (n > 0)
    {
        for (R_len_t i = 0; i < n; i++)
        {
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgTMatrix") == 0)
            {
                matclass = 0;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgCMatrix") == 0)
            {
                matclass = 1;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgRMatrix") == 0)
            {
                matclass = 2;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgeMatrix") == 0)
            {
                matclass = -1;
                break;
            }
        }
    }
    UNPROTECT(1);
    return matclass;
}

int sparselm_dgCMatrix_to_splm_ccsm(SEXP dgCMatrix, struct splm_ccsm *ccsm) {

	int nprotect = 0;

	if (!Rf_isS4(dgCMatrix)) { warning("Rf_isS4"); return 1; }
	if (sparselm_match_dg_class(dgCMatrix) != 1) { warning("dgCMatrix"); return 1; }
	if (!R_has_slot(dgCMatrix, Rf_install("Dim"))) { warning("Dim"); return 1; }
	if (!R_has_slot(dgCMatrix, Rf_install("i"))) { warning("i"); return 1; }
	if (!R_has_slot(dgCMatrix, Rf_install("p"))) { warning("p"); return 1; }
	if (!R_has_slot(dgCMatrix, Rf_install("x"))) { warning("x"); return 1; }
	
	SEXP Dim_slot = PROTECT(R_do_slot(dgCMatrix, Rf_install("Dim")));
	nprotect++;
	
	if (!Rf_isInteger(Dim_slot) || Rf_length(Dim_slot) != 2 || INTEGER_ELT(Dim_slot, 0) != ccsm->nr || INTEGER_ELT(Dim_slot, 1) != ccsm->nc) {
		warning("Dim_slot");
		UNPROTECT(nprotect);
		return 1;
	}
	
	SEXP i_slot = PROTECT(R_do_slot(dgCMatrix, Rf_install("i")));
	nprotect++;
	
	if (!Rf_isInteger(i_slot) || Rf_length(i_slot) != ccsm->nnz) {
		warning("i_slot");
		UNPROTECT(nprotect);
		return 1;
	}
	
	SEXP p_slot = PROTECT(R_do_slot(dgCMatrix, Rf_install("p")));
	nprotect++;
	
	if (!Rf_isInteger(p_slot) || Rf_length(p_slot) != ccsm->nc+1) {
		warning("p_slot");
		UNPROTECT(nprotect);
		return 1;
	}
	
	SEXP x_slot = PROTECT(R_do_slot(dgCMatrix, Rf_install("x")));
	nprotect++;
	
	if (!Rf_isReal(x_slot) || Rf_length(x_slot) != ccsm->nnz) {
		warning("x_slot");
		UNPROTECT(nprotect);
		return 1;
	}
	
	int i;
	int *i_ptr = INTEGER(i_slot);
	double *x_ptr = REAL(x_slot);
	int *p_ptr = INTEGER(p_slot);
	
	if (p_ptr[0] != 0 || p_ptr[ccsm->nc] != ccsm->nnz) {
		warning("p_slot");
		UNPROTECT(nprotect);
		return 1;
	}
	for (i = 0; i < ccsm->nc; i++) {
		if (p_ptr[i] < 0 || p_ptr[i] > ccsm->nnz || p_ptr[i] > p_ptr[i + 1]) {
			warning("p_slot");
			UNPROTECT(nprotect);
			return 1;
		}
	}
	
	for (i = 0; i < ccsm->nnz; i++) {
		if (i_ptr[i] < 0 || i_ptr[i] >= ccsm->nr) {
			warning("i_slot");
			UNPROTECT(nprotect);
			return 1;
		}
		ccsm->rowidx[i] = i_ptr[i];
		ccsm->val[i] = x_ptr[i];
	}
	
	for (i = 0; i <= ccsm->nc; i++) {
		ccsm->colptr[i] = p_ptr[i];
	}
	
	UNPROTECT(nprotect);
	
	return 0;
}

void sparselm_func(double *p, double *hx, int nvars, int nobs, void *adata) {

	SparseLMData *sparselm_data;
	
	sparselm_data = (SparseLMData*)adata;
	
	SEXP p_obj = PROTECT(allocVector(REALSXP, nvars));
	R_len_t i;
	for (i = 0; i < nvars; ++i) {
		REAL(p_obj)[i] = p[i];
	}
	SET_NAMES(p_obj, sparselm_data->pnames);
	
	SEXP R_fcall = PROTECT(lang2(sparselm_data->func, p_obj));
	SEXP hx_obj = PROTECT(eval(R_fcall, sparselm_data->rho));
	
	if (!IS_NUMERIC(hx_obj) || length(hx_obj) != nobs) {
        error("evaluation of fn function returns non-sensible value!");
    }
    
    for (i = 0; i < nobs; ++i) {
		hx[i] = REAL(hx_obj)[i];
	}
    
	UNPROTECT(3);
}

void sparselm_fjac(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata) {

	SparseLMData *sparselm_data;
	
	sparselm_data = (SparseLMData*)adata;
	
	SEXP p_obj = PROTECT(allocVector(REALSXP, nvars));
	R_len_t i;
	for (i = 0; i < nvars; ++i) {
		REAL(p_obj)[i] = p[i];
	}
	SET_NAMES(p_obj, sparselm_data->pnames);
	
	SEXP R_fcall = PROTECT(lang2(sparselm_data->fjac, p_obj));
	SEXP jac_obj = PROTECT(eval(R_fcall, sparselm_data->rho));
	
	if (sparselm_dgCMatrix_to_splm_ccsm(jac_obj, jac)) {
		error("fjac did not return expected dgCMatrix object");
	}
    
	UNPROTECT(3);
}

SEXP sparselm(SEXP func, SEXP fjac, SEXP p, SEXP x, SEXP nconvars, SEXP Jnnz, SEXP JtJnnz, SEXP itmax, SEXP opts, SEXP dif, SEXP rho) {

	int nprotect = 0;

	SparseLMData sparselm_data;
 	//double opts_vec[SPLM_OPTS_SZ];
 	double info_vec[SPLM_INFO_SZ];
	
	if (!isEnvironment(rho)) {
		Rf_error("rho is not an environment");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	
	sparselm_data.func = func;
	sparselm_data.fjac = fjac;
	sparselm_data.pnames = GET_NAMES(p);
	sparselm_data.rho = rho;
	
	if (!Rf_isReal(p) || Rf_length(p) < 1) {
		Rf_error("p is not numeric vector with at least one element");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	
	p = PROTECT(duplicate(p));
	nprotect++;
	
	int nvars = Rf_length(p);
	
	if (!Rf_isReal(x) || Rf_length(x) < 1) {
		Rf_error("x is not numeric vector with at least one element");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	
	int nobs = Rf_length(x);
	
	if (!Rf_isInteger(nconvars) || Rf_length(nconvars) != 1 || INTEGER(nconvars)[0] < 0 || INTEGER(nconvars)[0] >= nvars) {
		Rf_error("nconvars not integer of length 1 and less than nvars");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	if (!Rf_isInteger(Jnnz) || Rf_length(Jnnz) != 1) {
		Rf_error("Jnnz not positive integer of length 1");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	if (!Rf_isInteger(JtJnnz) || Rf_length(JtJnnz) != 1) {
		Rf_error("JtJnnz not positive integer of length 1");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	if (!Rf_isInteger(itmax) || Rf_length(itmax) != 1 || INTEGER(itmax)[0] < 1) {
		Rf_error("itmax not positive integer of length 1");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	if (!Rf_isLogical(dif) || Rf_length(dif) != 1) {
		Rf_error("dif not logical of length 1");
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	if (!Rf_isReal(opts) || Rf_length(opts) != SPLM_OPTS_SZ) {
		Rf_error("opts is not numeric vector of length %i", SPLM_OPTS_SZ);
		UNPROTECT(nprotect);
		return R_NilValue;
	}
	
	int niter;

	if (LOGICAL(dif)[0]) {
		niter = sparselm_difccs(&sparselm_func, &sparselm_fjac, REAL(p), REAL(x), nvars, INTEGER(nconvars)[0], nobs, INTEGER(Jnnz)[0], INTEGER(JtJnnz)[0], INTEGER(itmax)[0], REAL(opts), info_vec, (void *)&sparselm_data);
	} else {
		niter = sparselm_derccs(&sparselm_func, &sparselm_fjac, REAL(p), REAL(x), nvars, INTEGER(nconvars)[0], nobs, INTEGER(Jnnz)[0], INTEGER(JtJnnz)[0], INTEGER(itmax)[0], REAL(opts), info_vec, (void *)&sparselm_data);
	}
	
	SEXP niter_integer = PROTECT(NEW_INTEGER(1));
	nprotect++;
	INTEGER(niter_integer)[0] = niter;
	
	SEXP info_numeric = PROTECT(NEW_NUMERIC(SPLM_INFO_SZ));
	nprotect++;
	for (int i = 0; i < SPLM_INFO_SZ; i++) {
		REAL(info_numeric)[i] = info_vec[i];
	}
	
	SEXP out_list = PROTECT(NEW_LIST(3));
	nprotect++;
	
	SET_ELEMENT(out_list, 0, p);
	SET_ELEMENT(out_list, 1, niter_integer);
	SET_ELEMENT(out_list, 2, info_numeric);

	UNPROTECT(nprotect);

	return out_list;
}

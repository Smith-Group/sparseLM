#include <R.h>
#include <Rdefines.h>

SEXP add_assign_col_inplace(SEXP x, SEXP i, SEXP j, SEXP value, SEXP validate) {
	
	if (!Rf_isS4(x)) {
		Rf_error("x is not an S4 object");
	}
	
	SEXP class = PROTECT(Rf_getAttrib(x, R_ClassSymbol));
    R_len_t n = Rf_length(class);
    R_len_t ii;
    for (ii = 0; ii < n; ii++) {
        if (strcmp(CHAR(STRING_ELT(class, ii)), "dgCMatrix") == 0) {
			break;
		}
	}
    UNPROTECT(1);
	if (ii >= n) {
		Rf_error("x is not a dgCMatrix object");
	}
	
	if (!R_has_slot(x, Rf_install("Dim"))) {
		Rf_error("x does not have Dim slot");
	}
	if (!R_has_slot(x, Rf_install("i"))) {
		Rf_error("x does not have i slot");
	}
	if (!R_has_slot(x, Rf_install("p"))) {
		Rf_error("x does not have p slot");
	}
	if (!R_has_slot(x, Rf_install("x"))) {
		Rf_error("x does not have x slot");
	}
	
	int nprotect = 0;
	
	SEXP Dim_slot = PROTECT(R_do_slot(x, Rf_install("Dim")));
	nprotect++;
	
	if (!Rf_isInteger(Dim_slot) || Rf_length(Dim_slot) != 2) {
		UNPROTECT(nprotect);
		Rf_error("Dim slot of x is not integer vector of length 2");
	}
	
	SEXP i_slot = PROTECT(R_do_slot(x, Rf_install("i")));
	nprotect++;
	
	if (!Rf_isInteger(i_slot)) {
		UNPROTECT(nprotect);
		Rf_error("i slot of x is not integer vector");
	}
	
	SEXP p_slot = PROTECT(R_do_slot(x, Rf_install("p")));
	nprotect++;
	
	if (!Rf_isInteger(p_slot) || Rf_length(p_slot) != INTEGER_ELT(Dim_slot, 1)+1) {
		UNPROTECT(nprotect);
		Rf_error("p slot of x is not integer vector of length ncol(x)+1");
	}
	
	SEXP x_slot = PROTECT(R_do_slot(x, Rf_install("x")));
	nprotect++;
	
	if (!Rf_isReal(x_slot)) {
		UNPROTECT(nprotect);
		Rf_error("x slot of x is not numeric");
	}
	
	if (Rf_length(i_slot) != Rf_length(x_slot)) {
		UNPROTECT(nprotect);
		Rf_error("i and x slots have different lengths");
	}
	if (INTEGER_ELT(p_slot, Rf_length(p_slot) - 1) != Rf_length(x_slot)) {
		UNPROTECT(nprotect);
		Rf_error("last entry of p does not match length of x slot");
	}
	
	if (!Rf_isInteger(i) || Rf_length(i) > INTEGER_ELT(Dim_slot, 0)) {
		UNPROTECT(nprotect);
		Rf_error("i is not an integer vector of length less than or equal to dim(x)[1]");
	}

	if (!Rf_isInteger(j) || Rf_length(j) != 1 || INTEGER_ELT(j, 0) < 1 || INTEGER_ELT(j, 0) > Rf_length(p_slot)-1) {
		UNPROTECT(nprotect);
		Rf_error("j is not an integer vector of length 1 referencing a valid column");
	}
	
	if (!Rf_isReal(value) || Rf_length(value) != Rf_length(i)) {
		UNPROTECT(nprotect);
		Rf_error("value is not a numeric vector of length length(i)");
	}
	
	if (!Rf_isLogical(validate) || Rf_length(validate) != 1) {
		UNPROTECT(nprotect);
		Rf_error("validate is not a logical scalar");
	}
	int validate_flag = Rf_asLogical(validate);
	if (validate_flag == NA_LOGICAL) {
		UNPROTECT(nprotect);
		Rf_error("validate is NA");
	}
	
	int value_len = Rf_length(value);
	int *value_i = INTEGER(i);
	if (validate_flag) {
		int nrow = INTEGER_ELT(Dim_slot, 0);
		for (int k = 0; k < value_len; k++) {
			if (value_i[k] == NA_INTEGER) {
				UNPROTECT(nprotect);
				Rf_error("i contains NA");
			}
			if (value_i[k] < 1 || value_i[k] > nrow) {
				UNPROTECT(nprotect);
				Rf_error("i contains out-of-range row index");
			}
			if (k > 0 && value_i[k] <= value_i[k-1]) {
				UNPROTECT(nprotect);
				Rf_error("i must be strictly increasing");
			}
		}
		for (int k = 0; k < value_len; k++) {
			if (ISNA(REAL(value)[k]) || ISNAN(REAL(value)[k])) {
				UNPROTECT(nprotect);
				Rf_error("value contains NA or NaN");
			}
		}
		for (int k = 1; k < Rf_length(p_slot); k++) {
			if (INTEGER_ELT(p_slot, k) < INTEGER_ELT(p_slot, k - 1)) {
				UNPROTECT(nprotect);
				Rf_error("p slot is not nondecreasing");
			}
		}
	}
	
	int col_start_idx = INTEGER_ELT(p_slot, INTEGER_ELT(j, 0)-1);
	int col_end_idx = INTEGER_ELT(p_slot, INTEGER_ELT(j, 0));
	
	if (col_start_idx < 0 || col_start_idx > Rf_length(x_slot) || col_end_idx < 0 || col_end_idx > Rf_length(x_slot) ) {
		Rf_error("column bounds in p slot of x outside sparse data");
	}
	
	int col_len = col_end_idx - col_start_idx;
	int *col_i = INTEGER(i_slot) + col_start_idx;
	double *col_x = REAL(x_slot) + col_start_idx;
	double *value_x = REAL(value);
	
// 	for (int k = 0; k < col_len; k++) {
// 		printf("col_i[%i] = %i  col_x[%i] = %f\n", k, col_i[k], k, col_x[k]);
// 	}
// 	for (int k = 0; k < value_len; k++) {
// 		printf("value_i[%i] = %i  value_x[%i] = %f\n", k, value_i[k], k, value_x[k]);
// 	}
	
	int col_ii = 0;
	for (int value_ii = 0; value_ii < value_len; value_ii++) {
		int const value_iii = value_i[value_ii]-1;
// 		printf("value_iii = %i\n", value_iii);
		while (col_ii < col_len && col_i[col_ii] < value_iii) {
			col_ii++;
		}
		if (col_ii < col_len && col_i[col_ii] == value_iii) {
			col_x[col_ii] += value_x[value_ii];
		} else {
			UNPROTECT(nprotect);
			Rf_error("nonzero not present at x[%i,%i]", value_i[value_ii], INTEGER_ELT(j, 0));
		}
	}
	
	UNPROTECT(nprotect);
	
	return x;
}

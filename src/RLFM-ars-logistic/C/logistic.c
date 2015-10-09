#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>

#include "../../arslogistic/arsspline.h"
#include "util.h"

#define MY_REAL(x) (TYPEOF(x) == NILSXP ? NULL : REAL(x))
#define MY_INTEGER(x) (TYPEOF(x) == NILSXP ? NULL : INTEGER(x))

SEXP logistic_ars_Call (
	SEXP offset,
	SEXP beta_mean,
	SEXP beta_var,
	SEXP X,
	SEXP y,
	SEXP nFactors,
	SEXP nObs,
	SEXP beta_sample,
	SEXP qcent,
	SEXP ncent,
	SEXP ninit,
	SEXP x_lower,
	SEXP x_upper,
	SEXP nSamples
) {
	fitLogistic(MY_REAL(offset),
			    MY_REAL(beta_mean),
			    MY_REAL(beta_var),
			    MY_REAL(X),
			    MY_REAL(y),
			    MY_INTEGER(nFactors),
			    MY_INTEGER(nObs),
			    MY_REAL(beta_sample),
			    MY_REAL(qcent),
			    MY_INTEGER(ncent),
			    MY_INTEGER(ninit),
			    MY_REAL(x_lower),
			    MY_REAL(x_upper),
			    MY_INTEGER(nSamples)
			);
	return R_NilValue;
}

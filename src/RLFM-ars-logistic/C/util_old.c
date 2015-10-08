/*
        Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
        Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Liang Zhang
*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include "util.h"

void sym_eigen(const double* x, const int *nrow, double *eigen_val, double *eigen_vec){
    char jobz = 'V', uplo = 'L';
    double work_size, *work;
    int i,j, info, lwork=-1;
    
    for(i=0; i<(*nrow)*(*nrow); i++) eigen_vec[i] = x[i];

    F77_NAME(dsyev)(&jobz, &uplo, nrow, eigen_vec, nrow, eigen_val, &work_size, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");

    lwork = work_size;
    work  = (double*)Calloc(lwork, double);
    F77_NAME(dsyev)(&jobz, &uplo, nrow, eigen_vec, nrow, eigen_val, work, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");

    Free(work);
}

void sym_inv_byCholesky(
    double *A /* n x n matrix */, const int *n, const int *check_sym
){
    char uplo = 'L';
    int info, i, j;
    
    if(*check_sym > 0) CHK_SYMMETRIC(A, *n, i, j);
    
    F77_NAME(dpotrf)(&uplo, n, A, n, &info);
    if(info != 0) error("error in dpotrf(...): info=%d", info);
    F77_NAME(dpotri)(&uplo, n, A, n, &info);
    if(info != 0) error("error in dpotri(...): info=%d", info);
    for(i=1; i<(*n); i++){
        for(j=0; j<i; j++){
            A[C_MAT(j,i,*n)] = A[C_MAT(i,j,*n)];
        }
    }
}

void sym_inv3DA_byCholesky(
    double *invA /* k x n x n matrix */, const double *A /* k x n x n matrix */,
    const int *k, const int *n, double *temp /* n x n */, const int *check_sym
){
    for(int i=0; i<*k; i++){
        for(int j=0; j<*n; j++)
            for(int m=0; m<*n; m++) temp[C_MAT(j,m,*n)] = A[C_3DA(i,j,m,*k,*n)];
        sym_inv_byCholesky(temp, n, check_sym);
        
        for(int j=0; j<*n; j++)
            for(int m=0; m<*n; m++) invA[C_3DA(i,j,m,*k,*n)] = temp[C_MAT(j,m,*n)];
    }
}

void sum_margin(
    // OUTPUT
    double *ans,
    // INPUT
    const double *A, const int *nrow, const int *ncol, 
    const int *side // side=1: Sum up each row and return a vector with length nrow
                    // side=2: Sum up each column and return a vector with length ncol
){
    int i, j, end;
    if((*side) == 1){
        end=(*nrow)*(*ncol);
        for(i=0; i<*nrow; i++) ans[i] = 0;
        for(j=0; j<end; j+=(*nrow))
            for(i=0; i<*nrow; i++) ans[i] += A[i+j];
    }else if((*side) == 2){
        for(j=0; j<*ncol; j++){
            end = (j+1)*(*nrow);
            ans[j] = 0;
            for(i=j*(*nrow); i<end; i++) ans[j] += A[i];
        }
    }else{
        error("Unknown side=%d (please specify side=1 (for rows) or side=2 (for columns)");
    }
}

void center_array(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n // array length
		  ){
  // double sum = 0.0;
  double mean = x[0];
  // get sum
  //for(int i=0; i<n; i++){sum += x[i];}
  //mean = sum / (double) n;
  //get mean ... more stable ...
  for(int i=1; i<n; i++){ 
    mean = ( (double) i / (double) (i+1) ) * mean + ( 1.0 / (double) (i+1) ) * x[i];
  }

  //  Rprintf("mean = %f\n", mean);
  // center
  for(int i=0; i<n; i++){x[i] -= mean;}
}

void center_array_2d(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim // 1 for subtracting row mean, 2 for column
){

  double sum, mean;
  int i,j;

  if(dim == 1){
    for(i = 0; i<n; i++){
      sum = 0.0;
      for(j = 0; j<m; j++){sum += x[C_MAT(i, j, n)];}
      mean = sum / (double) m;
      for(j = 0; j<m; j++){x[C_MAT(i, j, n)] -= mean;}
    }
  } else if(dim == 2){
    for(j = 0; j<m; j++){
      //sum = 0.0;
      //for(i = 0; i<n; i++){sum += x[C_MAT(i, j, n)];}
      //mean = sum / (double) n;

      mean = x[C_MAT(0, j, n)];
      for(i=1; i<n; i++){
	mean = ((double)i/(double)(i+1)) * mean + (1.0/(double) (i+1))*x[C_MAT(i,j,n)];
      }
      //      Rprintf("fac mean col %d = %f\n", j, mean);
      for(i = 0; i<n; i++){x[C_MAT(i, j, n)] -= mean;}
  }

  }
}

void center_array_online(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n, // array length
		  int * oiNum
		  ){
  // double sum = 0.0;
  double mean = 0.0;
  int nav = 0;
  // get sum
  //for(int i=0; i<n; i++){sum += x[i];}
  //mean = sum / (double) n;
  //get mean ... more stable ...
  for(int i=0; i<n; i++){ 
    if(oiNum[i]>0){
      mean = ( (double) nav / (double) (nav+1) ) * mean +	\
	( 1.0 / (double) (nav+1) ) * x[i];
      nav += 1;
    }
  }

  //  Rprintf("mean = %f\n", mean);
  // center
  for(int i=0; i<n; i++){
    if(oiNum[i]>0){
      x[i] -= mean;
    }
  }
}

void center_array_2d_online(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim, // 1 for subtracting row mean, 2 for column
		     int * oiNum
){

  double sum, mean;
  int i,j;
  int nav;

  // dim == 1 is broken!!!
  if(dim == 1){
    for(i = 0; i<n; i++){
      if(oiNum[i]>0){
	sum = 0.0;
	for(j = 0; j<m; j++){sum += x[C_MAT(i, j, n)];}
	mean = sum / (double) m;
	for(j = 0; j<m; j++){x[C_MAT(i, j, n)] -= mean;}
      }
    }
  } else if(dim == 2){
    for(j = 0; j<m; j++){
      //sum = 0.0;
      //for(i = 0; i<n; i++){sum += x[C_MAT(i, j, n)];}
      //mean = sum / (double) n;
      nav = 0;
      mean = 0.0;
      for(i=0; i<n; i++){
	if(oiNum[i]>0){
	  mean = ((double)nav/(double)(nav+1)) * mean +	\
	    (1.0/(double) (nav+1))*x[C_MAT(i,j,n)];
	  nav += 1;
	}
      }
      //      Rprintf("fac mean col %d = %f\n", j, mean);
      for(i = 0; i<n; i++){
	if(oiNum[i]>0){
	  x[C_MAT(i, j, n)] -= mean;
	}
      }
    }

  }
}


void print_vector(char* prefix, double* vector, int length){
    int i;
    if(prefix != NULL) Rprintf("%s", prefix);
    for(i=0; i<length; i++)
        Rprintf("%f ", vector[i]);
    Rprintf("\n");
}

void print_matrix(const char* prefix, const double* matrix, const int nrow, const int ncol){
    int i,j;
    for(i=0; i<nrow; i++){
        if(prefix != NULL) Rprintf("%s", prefix);
        for(j=0; j<ncol; j++)
            Rprintf("%f\t", matrix[C_MAT(i,j,nrow)]);
        Rprintf("\n");
    }
}

// The intput/output are all R indices (start from 1, NOT 0)
void generateObsIndex(
    //OUTPUT
    int* obsIndex, // E.g., consider user i
    int* start,    //  y[ obsIndex[ start[i]+(0:(num[i]-1)) ] ]
    int* num,      //  are the observations of user i
    //INPUT
    const int* effIndex /* user or item */,
    const int* nObs, const int* nEff,
    //OTHER
    const int* debug
){
    int *ind, i, j, cIndex, eIndex, oIndex;
    ind = (int*)Calloc(*nEff, int);
    
    for(i=0; i<*nObs; i++){
        eIndex = effIndex[i];
        if(*debug > 0){
            CHK_R_INDEX(eIndex, *nEff);
        }
        num[R_VEC(eIndex)]++;
    }
    
    start[0] = 1; ind[0] = 1;
    for(i=1; i<*nEff; i++){
        start[i] = start[i-1]+num[i-1];
        ind[i] = start[i];
    }
    
    for(i=0; i<*nObs; i++){
        cIndex = R_VEC(effIndex[i]);
        obsIndex[R_VEC(ind[cIndex])] = i+1;
        ind[cIndex]++;
    }
    
    if(*debug > 0){
        for(i=0; i<*nEff; i++){
            if(ind[i] != start[i]+num[i]) error("logical error (level 1)");
            if(*debug > 1){
                for(j=0; j<num[i]; j++){
                    oIndex = obsIndex[R_VEC(start[i]+j)];
                    if(effIndex[R_VEC(oIndex)] != i+1)
                        error("logical error (level 2)");
                }
            }
        }
    }
    
    Free(ind);
}

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0
#define _(A)        A

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}

int nParam_debug = 0;
double* gr_debug = NULL;

void my_cgmin(int n, double *Bvec, double *X, double *Fmin,
	   optimfn fminfn, optimgr fmingr, int *fail,
	   double abstol, double intol, void *ex, int type, int trace,
	   int *fncount, int *grcount, int maxit)
{
    Rboolean accpoint;
    double *c, *g, *t;
    int count, cycle, cyclimit;
    double f;
    double G1, G2, G3, gradproj;
    int funcount=0, gradcount=0, i;
    double newstep, oldstep, setstep, steplength=1.0;
    double tol;

    if (maxit <= 0) {
	*Fmin = fminfn(n, Bvec, ex);
	*fncount = *grcount = 0;
	*fail = FALSE;
	return;
    }
    if (trace) {
	Rprintf("  Conjugate gradients function minimizer\n");
	switch (type) {
	case 1:	    Rprintf("Method: Fletcher Reeves\n");	break;
	case 2:	    Rprintf("Method: Polak Ribiere\n");		break;
	case 3:	    Rprintf("Method: Beale Sorenson\n");	break;
	default:
	    error(_("unknown 'type' in CG method of optim"));
	}
    }
    c = vect(n); g = vect(n); t = vect(n);

    // DEBUG
    nParam_debug = n;
    gr_debug = g;
    // END

    print_vector("INITIAL: ", gr_debug, nParam_debug);
    
    setstep = 1.7;
    *fail = 0;
    cyclimit = n;
    tol = intol * n * sqrt(intol);

    if (trace) Rprintf("tolerance used in gradient test=%g\n", tol);
    f = fminfn(n, Bvec, ex);
    print_vector("AFTER fminfn IN cgmin: ", gr_debug, nParam_debug);
    if (!R_FINITE(f)) {
	error(_("Function cannot be evaluated at initial parameters"));
    } else {
	*Fmin = f;
	funcount = 1;
	gradcount = 0;
	do {
	    for (i = 0; i < n; i++) {
		t[i] = 0.0;
		c[i] = 0.0;
	    }
	    cycle = 0;
	    oldstep = 1.0;
	    count = 0;
	    do {
		cycle++;
		if (trace) {
		    Rprintf("gradcount=%d, funcount=%d, Fmin=%f\n", gradcount, funcount, *Fmin);
		    Rprintf("parameters:\n");
		    for (i = 1; i <= n; i++) {
			Rprintf("%10.5f ", Bvec[i - 1]);
			if (i / 7 * 7 == i && i < n)
			    Rprintf("\n");
		    }
		    Rprintf("\n");
		}
		gradcount++;
		if (gradcount > maxit) {
		    *fncount = funcount;
		    *grcount = gradcount;
		    *fail = 1;
		    return;
		}
		if (trace) {
		    Rprintf("BEFORE: gradient:\n");
		    for (i = 1; i <= n; i++) {
    			Rprintf("%10.5f ", g[i - 1]);
    			if (i / 7 * 7 == i && i < n)
    			    Rprintf("\n");
		    }
		    Rprintf("\n");
		}
		fmingr(n, Bvec, g, ex);
		if (trace) {
		    Rprintf("AFTER: gradient:\n");
		    for (i = 1; i <= n; i++) {
    			Rprintf("%10.5f ", g[i - 1]);
    			if (i / 7 * 7 == i && i < n)
    			    Rprintf("\n");
		    }
		    Rprintf("\n");
		}
		G1 = 0.0;
		G2 = 0.0;
		for (i = 0; i < n; i++) {
		    X[i] = Bvec[i];
		    switch (type) {

		    case 1: /* Fletcher-Reeves */
			G1 += g[i] * g[i];
			G2 += c[i] * c[i];
			break;

		    case 2: /* Polak-Ribiere */
			G1 += g[i] * (g[i] - c[i]);
			G2 += c[i] * c[i];
			break;

		    case 3: /* Beale-Sorenson */
			G1 += g[i] * (g[i] - c[i]);
			G2 += t[i] * (g[i] - c[i]);
			break;

		    default:
			error(_("unknown type in CG method of optim"));
		    }
		    c[i] = g[i];
		}
		if (G1 > tol) {
		    if (G2 > 0.0)
			G3 = G1 / G2;
		    else
			G3 = 1.0;
		    gradproj = 0.0;
		    for (i = 0; i < n; i++) {
			t[i] = t[i] * G3 - g[i];
			gradproj += t[i] * g[i];
		    }
		    steplength = oldstep;

		    accpoint = FALSE;
		    do {
			count = 0;
			for (i = 0; i < n; i++) {
			    Bvec[i] = X[i] + steplength * t[i];
			    if (reltest + X[i] == reltest + Bvec[i])
				count++;
			}
			if (count < n) { /* point is different */
			    f = fminfn(n, Bvec, ex);
			    funcount++;
			    accpoint = (R_FINITE(f) &&
					f <= *Fmin + gradproj * steplength * acctol);

			    if (!accpoint) {
                    steplength *= stepredn;
                    if (trace){
                        Rprintf("* %f\n", f);
                    }
			    } else *Fmin = f; /* we improved, so update value */
			}
		    } while (!(count == n || accpoint));
		    if (count < n) {
			newstep = 2 * (f - *Fmin - gradproj * steplength);
			if (newstep > 0) {
			    newstep = -(gradproj * steplength * steplength / newstep);
			    for (i = 0; i < n; i++)
				Bvec[i] = X[i] + newstep * t[i];
			    *Fmin = f;
			    f = fminfn(n, Bvec, ex);
			    funcount++;
			    if (f < *Fmin) {
				*Fmin = f;
				if (trace) Rprintf("i< ");
			    } else { /* reset Bvec to match lowest point */
				if (trace) Rprintf("i> ");
				for (i = 0; i < n; i++)
				    Bvec[i] = X[i] + steplength * t[i];
			    }
			}
		    }
		}
		oldstep = setstep * steplength;
		if (oldstep > 1.0)
		    oldstep = 1.0;
	    } while ((count != n) && (G1 > tol) && (cycle != cyclimit));

	} while ((cycle != 1) ||
		 ((count != n) && (G1 > tol) && *Fmin > abstol));

    }
    if (trace) {
	Rprintf("Exiting from conjugate gradients minimizer\n");
	Rprintf("    %d function evaluations used\n", funcount);
	Rprintf("    %d gradient evaluations used\n", gradcount);
    }
    *fncount = funcount;
    *grcount = gradcount;
}

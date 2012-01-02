/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include "util.h"

void error2(const char *filename, int lineno, const char *fmt, ...){
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	error("ERROR in file %s at line %d: %s", filename, lineno, buf);
}

void do_nothing(const void* x){ }

// indices start from 1 (not 0)
void copy_double_array(double *out, const double *in, const int *size_out, const int *size_in, const int *start_out, const int *start_in, const int *length){
	if(*start_in  < 1  || *start_in  > *size_in)  error("in: start = %d, size = %d", *start_in, *size_in);
	if(*start_out < 1  || *start_out > *size_out) error("out: start = %d, size = %d", *start_out, *size_out);
	if(*length < 0 || (*start_in)+(*length)-1 > *size_in) error("in: start = %d, size = %d, length = %d", *start_in, *size_in, *length);
	if((*start_out)+(*length)-1 > *size_out) error("out: start = %d, size = %d, length = %d", *start_out, *size_out, *length);
	memcpy(&(out[(*start_out)-1]), &(in[(*start_in)-1]), sizeof(double)*(*length));
}
void copy_int_array(int *out, const int *in, const int *size_out, const int *size_in, const int *start_out, const int *start_in, const int *length){
	if(*start_in  < 1  || *start_in  > *size_in)  error("in: start = %d, size = %d", *start_in, *size_in);
	if(*start_out < 1  || *start_out > *size_out) error("out: start = %d, size = %d", *start_out, *size_out);
	if(*length < 0 || (*start_in)+(*length)-1 > *size_in) error("in: start = %d, size = %d, length = %d", *start_in, *size_in, *length);
	if((*start_out)+(*length)-1 > *size_out) error("out: start = %d, size = %d, length = %d", *start_out, *size_out, *length);
	memcpy(&(out[(*start_out)-1]), &(in[(*start_in)-1]), sizeof(int)*(*length));
}

void check_double_matrix_symmetric(const double *A, const int n){
	if(n <= 1) return;
	double max = A[0];
	for(int i=1; i<n*n; i++){
		double A_i = fabs(A[i]);
		if(A_i > max) max = A_i;
	}
    for(int i=1; i<n; i++){
        for(int j=0; j<i; j++){
            double diff = fabs(A[C_MAT(j,i,n)] - A[C_MAT(i,j,n)]);
            if(diff/max > 1e-8) STOP3("Matrix is not symmetric: %.16g vs %.16g (relative diff = %.16g)", A[C_MAT(j,i,n)], A[C_MAT(i,j,n)], diff/max);
        }
    }
}

void sym_eigen(const double* x, const int *nrow, double *eigen_val, double *eigen_vec){
    char jobz = 'V', uplo = 'L';
    double work_size, *work;
    int i, info, lwork=-1;
    
    for(i=0; i<(*nrow)*(*nrow); i++) eigen_vec[i] = x[i];

    F77_NAME(dsyev)(&jobz, &uplo, nrow, eigen_vec, nrow, eigen_val, &work_size, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");

    lwork = (int)work_size;
    work  = (double*)Calloc(lwork, double);
    F77_NAME(dsyev)(&jobz, &uplo, nrow, eigen_vec, nrow, eigen_val, work, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");

    Free(work);
}
void sym_eigen2(const double* x, const int *nrow, double *eigen_val, double *eigen_vec, double *workspace, const int *workspace_size, const int *check_sym){
    char jobz = 'V', uplo = 'L';
    int i, info;
    if(*check_sym > 0) CHK_SYMMETRIC(x, *nrow);
    for(i=0; i<(*nrow)*(*nrow); i++) eigen_vec[i] = x[i];
    F77_NAME(dsyev)(&jobz, &uplo, nrow, eigen_vec, nrow, eigen_val, workspace, workspace_size, &info);
    if(info != 0) error("error in dsyev(...)");
}
int sym_eigen2_workspace_size(const int *nrow){
    char jobz = 'V', uplo = 'L';
    double work_size;
    int info, lwork=-1;
    F77_NAME(dsyev)(&jobz, &uplo, nrow, NULL, nrow, NULL, &work_size, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");
    lwork = (int)work_size;
    return lwork;
}


void sym_inv_byCholesky(
    double *A /* n x n matrix */, const int *n, const int *check_sym
){
    char uplo = 'L';
    int info, i, j;
    
    if(*check_sym > 0) CHK_SYMMETRIC(A, *n);
    
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

void print_vector(const char* prefix, const double* vector, const int length){
    int i;
    if(prefix != NULL) Rprintf("%s", prefix);
    for(i=0; i<length; i++)
        Rprintf(" %11.7f", vector[i]);
    Rprintf("\n");
}

void print_intVector(const char* prefix, const int* vector, const int length){
    int i;
    if(prefix != NULL) Rprintf("%s", prefix);
    for(i=0; i<length; i++)
        Rprintf("%d ", vector[i]);
    Rprintf("\n");
}

void print_matrix(const char* prefix, const double* matrix, const int nrow, const int ncol){
    int i,j;
    for(i=0; i<nrow; i++){
        if(prefix != NULL) Rprintf("%s", prefix);
        for(j=0; j<ncol; j++)
            Rprintf(" %11.7f", matrix[C_MAT(i,j,nrow)]);
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
    
    for(i=0; i<*nEff; i++) num[i] = 0;
    
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

void normalizeToSumUpToOne2(double *output, const double *input, const int length){
    double sum = 0;
    for(int k=0; k<length; k++) sum += input[k];
    for(int k=0; k<length; k++) output[k] = input[k] / sum;
}

void normalizeToSumUpToOne(double *vector, const int length){
    double sum = 0;
    for(int k=0; k<length; k++) sum += vector[k];
    for(int k=0; k<length; k++) vector[k] /= sum;
}

void indexWithQuantities(int *output, int *vector, const int *length){
    int k=0;
    for(int i=0; i<*length; i++){
        for(int j=0; j<vector[i]; j++){ output[k] = i+1; k++;}
    }
}

/**
 * out[k,j] = sum_{x s.t. groupBy[x] = j} matrix[k,select[x]] * weight[x]
 *
 * select and groupBy are R indices (starting from 1)
 */
void selectColumn_agg_sum(
    double *out, const int *nrowOut, const int *ncolOut,
    const double *matrix, const int *nrow, const int *ncol,
    const int *select, const int *groupBy, const int *num,
    const double *weight, const int *nWeights
){
    if(*nrowOut != *nrow) error("error in %s at line %d", __FILE__, __LINE__);
    if((*nWeights) != 0 && (*nWeights) != (*num)) error("error in %s at line %d", __FILE__, __LINE__);
    for(int k=0; k<(*nrowOut)*(*ncolOut); k++) out[k] = 0;
    for(int x=0; x<*num; x++){
        int groupIndex = groupBy[x]-1;
        int colIndex   = select[x]-1;
        CHK_C_INDEX(groupIndex, *ncolOut);
        CHK_C_INDEX(colIndex,   *ncol);
        double w = 1;
        if((*nWeights) == (*num)) w = weight[x];
        for(int k=0; k<*nrow; k++) out[C_MAT(k,groupIndex,*nrow)] += matrix[C_MAT(k,colIndex,*nrow)] * w;
    }
}

/**
 * margin = 1: each row sum up to one
 * margin = 2: each column sum up to one
 */
void normalize_sumToOne2D(
    double *out, const double *matrix, 
    const int *nrow, const int *ncol, const int *margin
){
    if(*margin == 1){
        for(int i=0; i<*nrow; i++){
            double sum = 0;
            for(int j=0; j<*ncol; j++) sum += matrix[C_MAT(i,j,*nrow)];
            if(sum != 0){
                for(int j=0; j<*ncol; j++) out[C_MAT(i,j,*nrow)] = matrix[C_MAT(i,j,*nrow)] / sum;
            }else{
                for(int j=0; j<*ncol; j++) out[C_MAT(i,j,*nrow)] = 0;
            }
        }
    }else if(*margin == 2){
        for(int j=0; j<*ncol; j++){
            double sum = 0;
            for(int i=0; i<*nrow; i++) sum += matrix[C_MAT(i,j,*nrow)];
            if(sum != 0){
                for(int i=0; i<*nrow; i++) out[C_MAT(i,j,*nrow)] = matrix[C_MAT(i,j,*nrow)] / sum;
            }else{
                for(int i=0; i<*nrow; i++) out[C_MAT(i,j,*nrow)] = 0;
            }
        }        
    }else DIE_HERE;
}

/**
 *  output[i] = input[i] / sum_{j : by[j]=by[i]} input[j]
 *  if 0 = sum_{j : by[j]=by[i]} input[j], then output 0
 */
void normalize_sumToOne_groupby(
	double *output, const double *input, const int *by,
	const int *length
){
	int max;
	if((*length) <= 0) return;
	max = by[0];
	for(int i=0; i<(*length); i++){
		if(by[i]<0) error("by[%d] = %d (cannot < 0)", i, by[i]);
		if(input[i]<0) error("input[%d] = %f (should not < 0)", i, input[i]);
		if(max < by[i]) max = by[i];
	}
	double* sum = (double*)Calloc(max+1,double);

	for(int i=0; i<(*length); i++) sum[by[i]] += input[i];
	for(int i=0; i<(*length); i++) output[i] = (sum[by[i]] == 0) ? 0 : input[i]/sum[by[i]];

	Free(sum);
}


void print_doublePointer(double *pointer){
    Rprintf("Address: %p\n",pointer);
}

void print_intPointer(int *pointer){
    Rprintf("Address: %p\n",pointer);
}

void get_doublePointer(double *pointer, int *address){
	(*address) = (long)pointer;
}

void get_intPointer(int *pointer, int *address){
	(*address) = (long)pointer;
}

void compute_uBv_dense(
    // Output
    double *score, // nObs x 1
    // Input
    const int *userIndex, const int *itemIndex, // Both nObs x 1; index starts from 1 (not 0)
    const double *u, // nUsers x nUserFactors
    const double *B, // nUserFeatures x nItemFeatures
    const double *v, // nItems x nItemFactors
    const int *nObs, const int *nUsers, const int *nItems,
    const int *nUserFactors, const int *nItemFactors
){
    for(int k=0; k<(*nObs); k++){
        int i = userIndex[k];
        int j = itemIndex[k];
        if(i < 1 || i > (*nUsers)) error("userIndex out of bound\n");
        if(j < 1 || j > (*nItems)) error("itemIndex out of bound\n");
        score[k] = 0;
        for(int p=1; p<=(*nUserFactors); p++){
            for(int q=1; q<=(*nItemFactors); q++){
                score[k] += B[R_MAT(p,q,*nUserFactors)] * u[R_MAT(i,p,*nUsers)] * v[R_MAT(j,q,*nItems)];
            }
        }
    }
}

inline double compute_szuBv_c_single_dense(
	int i, // user i (start from 1)
	int j, // item j (start from 1)
	const double *s, // nUsers x nUserClusters
	const double *z, // nItems x nItemClusters
    const double *u, // nUsers x nFactorsPerUser
	const double *B, // nUserClusters x nItemClusters x nFeaturesPerUser x nFactorsPerItem
    const double *v, // nItems x nFactorsPerItem
    const double *c, // nUserClusters x nItemClusters
    const int nUsers, const int nItems,
    const int nUserClusters, const int nItemClusters,
    const int nFactorsPerUser, const int nFactorsPerItem,
    const int identity_B, const int has_c
){
	double ans = 0, uzsB, zsB, sB;
	i--; j--;
	int x = 0;
	if(identity_B){
		if(nFactorsPerItem != nFactorsPerUser) error("nFactorsPerItem != nFactorsPerUser (%d vs. %d)", nFactorsPerItem, nFactorsPerUser);
		for(int f=0; f<nFactorsPerUser; f++) ans += u[C_MAT(i,f,nUsers)] * v[C_MAT(j,f,nItems)];
	}else{
		for(int n=0; n<nFactorsPerItem; n++){
			uzsB = 0;
			for(int m=0; m<nFactorsPerUser; m++){
				zsB = 0;
				for(int k=0; k<nItemClusters; k++){
					sB = 0;
					for(int ell=0; ell<nUserClusters; ell++){
						sB += s[C_MAT(i,ell,nUsers)] * B[x];
						x++;
					}
					zsB += z[C_MAT(j,k,nItems)] * sB;
				}
				uzsB += u[C_MAT(i,m,nUsers)] * zsB;
			}
			ans += v[C_MAT(j,n,nItems)] * uzsB;
		}
	}
	if(has_c){
		for(int ell=0; ell<nUserClusters; ell++){
			for(int k=0; k<nItemClusters; k++){
				ans += s[C_MAT(i,ell,nUsers)] * z[C_MAT(j,k,nItems)] * c[C_MAT(ell,k,nUserClusters)];
			}
		}
	}
	return ans;
}

void compute_szuBv_c_dense(
	double *ans, // nObs x 1
	const int* userIndex, // nObs x 1 (start from 1)
	const int* itemIndex, // nObs x 1 (start from 1)
	const double *s, // nUsers x nUserClusters
	const double *z, // nItems x nItemClusters
	const double *u, // nUsers x nFactorsPerUser
	const double *B, // nUserClusters x nItemClusters x nFeaturesPerUser x nFactorsPerItem
	const double *v, // nItems x nFactorsPerItem
	const double *c, // nUserClusters x nItemClusters
	const int *nObs, const int *nUsers, const int *nItems,
	const int *nUserClusters, const int *nItemClusters,
	const int *nFactorsPerUser, const int *nFactorsPerItem,
	const int *identity_B, const int *has_c,
	const int *debug
){
	for(int m=0; m<*nObs; m++){
		int i = userIndex[m];
		int j = itemIndex[m];
		if((*debug) > 0){
			CHK_R_INDEX(i,*nUsers);
			CHK_R_INDEX(j,*nItems);
		}
		ans[m] = compute_szuBv_c_single_dense(
					i,j,s,z,u,B,v,c,*nUsers,*nItems,*nUserClusters,
					*nItemClusters,*nFactorsPerUser,*nFactorsPerItem,
					*identity_B, *has_c
		);
	}
}

// Let (i,j) = (index_u[k], index_v[k])
//   For m in 1:(ncol_u * ncol_v)
//	 	X[k,m] = X[k,(p,q)] = u[i,p] * v[j,q]
//   For m in (ncol_u * ncol_v) + 1:ncol_x
//	 	X[k,m] = x[k, m - (ncol_u * ncol_v)]
//
// W = diag(w) * diag(w)
// XWX = t(X) %*% W %*% X
// XWY = t(X) %*% W %*% y
void compute_XWX_XWy_dyadic(
	double *XWX, // (ncol_u * ncol_v + ncol_x) x (ncol_u * ncol_v + ncol_x)
	double *XWy, // (ncol_u * ncol_v + ncol_x)
	const double *y,    // nObx x 1
	const int *index_u, // nObs x 1 (starting from 1)
	const int *index_v, // nObs x 1 (starting from 1)
	const double *u, // nrow_u x ncol_u
	const double *v, // nrow_v x ncol_v
	const double *x, // nrow_x x ncol_x (nrow_x = nObs)
	const double *w, // nWeights x 1
	const int *nObs,
	const int *nrow_u, const int *ncol_u, const int *nrow_v, const int *ncol_v,
	const int *nrow_x, const int *ncol_x,
	const int *nWeights, const double *tol, const int *debug
){
	if((*nWeights) != 0 && (*nWeights) != (*nObs)) error("(*nWeights) != 0 && (*nWeights) != *nObs");
	if((*nrow_x) != 0 && (*nrow_x) != (*nObs)) error("(*nrow_x) != 0 && (*nrow_x) != (*nObs)");
	if(((*nrow_x) != 0 && (*ncol_x) == 0) || ((*nrow_x) == 0 && (*ncol_x) != 0))
		error("Problem with nrow_x and ncol_x");
	int len = (*ncol_u) * (*ncol_v) + (*ncol_x);
	double *wx = (double*)Calloc(len, double);

	for(int k=0; k<len*len; k++) XWX[k] = 0;
	for(int k=0; k<len; k++) XWy[k] = 0;

	for(int k=0; k<(*nObs); k++){
		double thisWeight = 1;
		if((*nWeights) > 0) thisWeight = w[k];
		if(thisWeight < 0) error("w[%d] = %f < 0", k+1, thisWeight);
		if(thisWeight < (*tol)) continue;

		int i = index_u[k] - 1;
		int j = index_v[k] - 1;
		if((*debug) > 0){
			CHK_C_INDEX(i,*nrow_u);
			CHK_C_INDEX(j,*nrow_v);
		}
		int m = 0;
		for(int q=0; q<*ncol_v; q++){
			for(int p=0; p<*ncol_u; p++){
				wx[m] = thisWeight * u[C_MAT(i,p,*nrow_u)] * v[C_MAT(j,q,*nrow_v)];
				m++;
			}
		}
		if((*ncol_x) > 0){
			for(int n=0; n<(*ncol_x); n++){
				wx[m] = thisWeight * x[C_MAT(k,n,*nrow_x)];
				m++;
			}
		}
		if(m!=len) error("m!=len");
		for(m=0; m<len; m++){
			for(int n=0; n<len; n++){
				XWX[C_MAT(m,n,len)] += wx[m] * wx[n];
			}
		}
		double wy = thisWeight * y[k];
		for(m=0; m<len; m++){
			XWy[m] += wx[m] * wy;
		}
	}
	Free(wx);
}

int rdiscrete(double* probabilities, int nOutcomes, double *rnd_out){
    double rnd = unif_rand();
    if(rnd_out != NULL) (*rnd_out) = rnd;
    for(int i=0; i<nOutcomes; i++){
        if(rnd < probabilities[i]) return i;
        rnd -= probabilities[i];
        if(rnd < -1e-12) error("ERROR: In draw (%s:%d), probabilities do not sum up to one\n", __FILE__, __LINE__);
    }
    if(rnd > 1e-12) error("ERROR: In draw (%s:%d), probabilities do not sum up to one\n", __FILE__, __LINE__);
    return nOutcomes-1;
}

void rdiscrete_many(
	int *output, // nCases x nRep: values in {1, ..., nOutComes}
	const double* probabilities, // nCases x nOutComes
	const int *nCases, const int *nRep, const int *nOutComes
){
	GetRNGstate();
	double *temp = (double*)Calloc(*nOutComes, double);
	for(int i=0; i<*nCases; i++){
		double sum = 0;
		for(int k=0; k<*nOutComes; k++){
			temp[k] = probabilities[C_MAT(i,k,*nCases)];
			sum += temp[k];
		}
		if(abs(sum - 1) > 1e-12) error("ERROR: Probabilities do not sum up to one!");
		for(int k=0; k<*nRep; k++) output[C_MAT(i,k,*nCases)] = rdiscrete(temp,*nOutComes,NULL)+1;
	}
	Free(temp);
	PutRNGstate();
}

// Let (i,j) = (index_u[k], index_v[k])
// X[k,m] = X[k,(p,q)] = u[i,p] * v[j,q]
// W = diag(w) * diag(w)
// XWX = t(X) %*% W %*% X
// XWY = t(X) %*% W %*% y
void compute_XWX_XWy_dyadic_old(
	double *XWX, // (ncol_u * ncol_v) x (ncol_u * ncol_v)
	double *XWy, // (ncol_u * ncol_v)
	const double *y,    // nObx x 1
	const int *index_u, // nObs x 1 (starting from 1)
	const int *index_v, // nObs x 1 (starting from 1)
	const double *u, // nrow_u x ncol_u
	const double *v, // nrow_v x ncol_v
	const double *w, // nWeights x 1
	const int *nObs,
	const int *nrow_u, const int *ncol_u, const int *nrow_v, const int *ncol_v,
	const int *nWeights, const double *tol, const int *debug
){
	if((*nWeights) != 0 && (*nWeights) != (*nObs)) error("(*nWeights) != 0 && (*nWeights) != *nObs");
	int len = (*ncol_u) * (*ncol_v);
	double *wx = (double*)Calloc(len, double);

	for(int k=0; k<len*len; k++) XWX[k] = 0;
	for(int k=0; k<len; k++) XWy[k] = 0;

	for(int k=0; k<(*nObs); k++){
		double thisWeight = 1;
		if((*nWeights) > 0) thisWeight = w[k];
		if(thisWeight < 0) error("w[%d] = %f < 0", k+1, thisWeight);
		if(thisWeight < (*tol)) continue;

		int i = index_u[k] - 1;
		int j = index_v[k] - 1;
		if((*debug) > 0){
			CHK_C_INDEX(i,*nrow_u);
			CHK_C_INDEX(j,*nrow_v);
		}
		int m = 0;
		for(int q=0; q<*ncol_v; q++){
			for(int p=0; p<*ncol_u; p++){
				wx[m] = thisWeight * u[C_MAT(i,p,*nrow_u)] * v[C_MAT(j,q,*nrow_v)];
				m++;
			}
		}
		for(m=0; m<len; m++){
			for(int n=0; n<len; n++){
				XWX[C_MAT(m,n,len)] += wx[m] * wx[n];
			}
		}
		double wy = thisWeight * y[k];
		for(m=0; m<len; m++){
			XWy[m] += wx[m] * wy;
		}
	}
	Free(wx);
}

/**
 * Important note: XWX and XWy will be changed by this function!!
 * Regression weight beta = (X'WX + prior_var^{-1})^{-1} %*% (X'WY + prior_var^{-1} * prior_mean)
 */
void compute_weight_from_XWX_XWy(
	double *beta, // nFeatures x 1 (OUTPUT: Posterior Mean)
	double *XWX,  // nFeatures x nFeatures (OUTPUT: Posterior Var)
	double *XWy,  // nFeatures x 1
	const double *prior_mean, // nFeatures x 1
	const double *prior_var,  // 1x1
	const int *nFeatures, const int *hasPrior,
	const int *debug
){
	if(*hasPrior){
		for(int i=0; i<*nFeatures; i++){
			double lambda = (1.0 / prior_var[0]);
			XWX[C_MAT(i,i,*nFeatures)] += lambda;
			XWy[i] += lambda * prior_mean[i];
		}
	}
	sym_inv_byCholesky(XWX, nFeatures, debug);
	for(int i=0; i<*nFeatures; i++){
		beta[i] = 0;
		for(int j=0; j<*nFeatures; j++) beta[i] += XWX[C_MAT(i,j,*nFeatures)] * XWy[j];
	}
}

/**
 * Bayesian gaussian linear regression
 */
void bayesian_gaussian_regression(
	double *postMean,     // (output)   nFeatures x 1
	double *postVar,      // (output)   nFeatures x nFeatures
	const double *y,      // (response) nObs x 1
	const double *X,      // (features) nObs x nFeatures
	const double *offset, // nObs x 1
	const double *prior_mean, // nFeatures x 1
	const double *prior_var,  // nFeatures x nFeatures
	const double *var_y,      // nObs x 1
	const int *nObs, const int *nFeatures, const int *nOffset,
	const int *debug, const int *verbose
){
	if((*nOffset) != 0 && (*nOffset) != (*nObs)) error("#Offsets should be either 0 or #Observations: %f vs. %f", *nOffset, *nObs);

	double *XWX  = postVar;
	double *XWy  = (double*)Calloc((*nFeatures), double);
	double *x  = (double*)Calloc((*nFeatures), double);
	double *inv_prior_var = (double*)Calloc((*nFeatures)*(*nFeatures), double);
	int nFeatures_sqr = (*nFeatures) * (*nFeatures);

	for(int p=0; p<nFeatures_sqr; p++) inv_prior_var[p] = prior_var[p];
	sym_inv_byCholesky(inv_prior_var, nFeatures, debug);

	for(int p=0; p<nFeatures_sqr; p++) XWX[p] = 0;
	for(int p=0; p<*nFeatures; p++){
		XWy[p] = 0;
		for(int q=0; q<*nFeatures; q++){
			XWy[p] += inv_prior_var[C_MAT(p,q,*nFeatures)] * prior_mean[q];
			XWX[C_MAT(p,q,*nFeatures)] = inv_prior_var[C_MAT(p,q,*nFeatures)];
		}
	}

	if((*verbose) > 10){
		Rprintf("Init XWX:\n"); print_matrix(NULL,XWX,*nFeatures,*nFeatures);
		Rprintf("Init XWy:\n"); print_vector(NULL, XWy, *nFeatures);
	}

	for(int k=0; k<*nObs; k++){
		for(int p=0; p<(*nFeatures); p++) x[p] = X[C_MAT(k,p,*nObs)];
		// Update X'WX & X'Wy
		double y_new = y[k];
		if((*nOffset) == (*nObs)) y_new -= offset[k];
		for(int p=0; p<(*nFeatures); p++){
			for(int q=0; q<(*nFeatures); q++){
				XWX[C_MAT(p,q,*nFeatures)] += (x[p] * x[q]) / var_y[k];
			}
			XWy[p] += (x[p] * y_new) / var_y[k];
		}
	}
	int zero = 0;
	compute_weight_from_XWX_XWy(postMean, XWX, XWy, NULL, NULL, nFeatures, &zero, debug);
	if((*verbose) > 10){
		Rprintf("Data:\n");
		for(int i=0; i<*nObs; i++){
			Rprintf(" %11.7f  %11.7f", y[i], offset[i]);
			for(int j=0; j<*nFeatures; j++) Rprintf(" %11.7f", X[C_MAT(i,j,*nObs)]);
			Rprintf("\n");
		}
		Rprintf("XWX:\n"); print_matrix(NULL,XWX,*nFeatures,*nFeatures);
		Rprintf("XWy:\n"); print_vector(NULL, XWy, *nFeatures);
		Rprintf("Regression Weights:\n"); print_vector(NULL, postMean, *nFeatures);
	}

	Free(XWy);
	Free(x);
	Free(inv_prior_var);
}

void computeMeanSumvar(
    // OUTPUT
    double *mean, double *sumvar,
    // INPUT
    const double *sum, const double *sos /* sum of squares */, const int length, const int nSamples
){
    int k;
    double ratio = ((double)nSamples) / (nSamples - 1.0);
    sumvar[0] = 0;
    for(k=0; k<length; k++){
        mean[k] = sum[k] / nSamples;
        sumvar[0] += (sos[k] / (nSamples - 1.0)) - (ratio * mean[k] * mean[k]);
    }
}

void computeMeanVar(
    // OUTPUT
    double *mean /*OUT*/, double *sumvar /*OUT*/, double *outputVar /*IN:sum-of-products; OUT:variance*/,
    // INPUT
    const double *sum, const int nEffects, const int nFactors, const int nSamples
){
    int k;
    double ratio = ((double)nSamples) / (nSamples - 1.0);
    sumvar[0] = 0;
    for(k=0; k<nEffects*nFactors; k++) mean[k] = sum[k] / nSamples;
    for(k=0; k<nEffects; k++){
        if(nFactors == 1){
            outputVar[k] = (outputVar[k] / (nSamples - 1.0)) - (ratio * mean[k] * mean[k]);
            sumvar[0] += outputVar[k];
        }else{
            for(int f1=0; f1<nFactors; f1++){
                for(int f2=0; f2<=f1; f2++){
                    outputVar[C_3DA(k,f1,f2,nEffects,nFactors)] = (outputVar[C_3DA(k,f1,f2,nEffects,nFactors)] / (nSamples - 1.0)) - (ratio * mean[C_MAT(k,f1,nEffects)] * mean[C_MAT(k,f2,nEffects)]);
                    if(f1 != f2){
                        outputVar[C_3DA(k,f2,f1,nEffects,nFactors)] = outputVar[C_3DA(k,f1,f2,nEffects,nFactors)];
                    }
                }
                sumvar[0] += outputVar[C_3DA(k,f1,f1,nEffects,nFactors)];
            }
        }
    }
}

void computeCov(
    // OUTPUT
    double *outputCov /*IN:sum-of-products; OUT:covariance; size: num1 x num2*/,
    // INPUT
    const double *mean1, const double *mean2,
    const int num1, const int num2, const int nSamples
){
	for(int i=0; i<num1; i++) for(int k=0; k<num2; k++){
		outputCov[C_MAT(i,k,num1)] = outputCov[C_MAT(i,k,num1)]/nSamples - mean1[i]*mean2[k];
	}
}


// X is n x n
inline double min_diag(double *X, int n){
	if(n <= 0) error("Number of elements is <= zero when computing minimum");
	double ans = X[0];
	for(int i=1; i<n; i++) if(X[C_MAT(i,i,n)] < ans) ans = X[C_MAT(i,i,n)];
	return ans;
}

static int ok_to_discount(const double *XWX, const double discount, const double *inv_prior_var, const int nFeatures){
	for(int p=0; p<nFeatures; p++)
		if(XWX[C_MAT(p,p,nFeatures)]*discount < inv_prior_var[C_MAT(p,p,nFeatures)]) return 0;
	return 1;
}

/**
 * Weight of the examples in the past batches will be discounted by this number
 * after each batch
 */
void online_gaussian_batch_predict(
	double *prediction,   // nObs x 1 (output)
	double *regWeight,    // nBatches x nFeatures (output)
	const double *y,      // nObs x 1 (response)
	const double *X,      // nObs x nFeatures
	const double *offset, // nObs x 1
	const int *batch_id,  // nObs x 1 (ID needs to be sorted, start from 1)
	const double *prior_mean, // nFeatures x 1
	const double *prior_var,  // nFeatures x nFeatures
	const double *discount,   // 1x1 (0 <= discount <= 1)
	const int *nObs, const int *nFeatures, const int *nOffset, const int *nBatches,
	const int *output_regWeight, const int *debug, const int *verbose
){
	int prev_id = 1, zero = 0;
	if((*discount) < 0 || (*discount) > 1) error("The discounting factor should be between 0 and 1: %f", *discount);
	if((*nOffset) != 0 && (*nOffset) != (*nObs)) error("#Offsets should be either 0 or #Observations: %f vs. %f", *nOffset, *nObs);

	double *XWX  = (double*)Calloc((*nFeatures)*(*nFeatures), double);
	double *XWy  = (double*)Calloc((*nFeatures), double);
	double *x  = (double*)Calloc((*nFeatures), double);
	double *beta = (double*)Calloc((*nFeatures), double);
	double *XWX_temp = (double*)Calloc((*nFeatures)*(*nFeatures), double);
	double *XWy_temp = (double*)Calloc((*nFeatures), double);
	double *inv_prior_var = (double*)Calloc((*nFeatures)*(*nFeatures), double);
	int nFeatures_sqr = (*nFeatures) * (*nFeatures);

	for(int p=0; p<nFeatures_sqr; p++) inv_prior_var[p] = prior_var[p];
	sym_inv_byCholesky(inv_prior_var, nFeatures, debug);

	for(int p=0; p<nFeatures_sqr; p++) XWX[p] = 0;
	for(int p=0; p<*nFeatures; p++){
		beta[p] = prior_mean[p];
		XWy[p]  = 0;
		for(int q=0; q<*nFeatures; q++){
			XWy[p] += inv_prior_var[C_MAT(p,q,*nFeatures)] * prior_mean[q];
			XWX[C_MAT(p,q,*nFeatures)] = inv_prior_var[C_MAT(p,q,*nFeatures)];
		}
	}

	if((*verbose) > 10){
		Rprintf("Init XWX:\n"); print_matrix(NULL,XWX,*nFeatures,*nFeatures);
		Rprintf("Init XWy:\n"); print_vector(NULL, XWy, *nFeatures);
	}

	for(int k=0; k<*nObs; k++){
		if(batch_id[k] <= 0) error("Batch ID starts from 1: %d", batch_id[k]);
		if(batch_id[k] > *nBatches) error("Batch ID cannot be > %d: %d", *nBatches, batch_id[k]);
		if(prev_id > batch_id[k]) error("Data is not ordered by batch ID");
		if(prev_id != batch_id[k]){
			if(k > 0){
				// Update the model
				for(int p=0; p<nFeatures_sqr; p++) XWX_temp[p] = XWX[p];
				for(int p=0; p<(*nFeatures); p++) XWy_temp[p] = XWy[p];
				compute_weight_from_XWX_XWy(beta, XWX_temp, XWy_temp, NULL, NULL, nFeatures, &zero, debug);
				// Apply discount
				if(ok_to_discount(XWX, *discount, inv_prior_var, *nFeatures)){
					for(int p=0; p<nFeatures_sqr; p++) XWX[p] *= (*discount);
					for(int p=0; p<(*nFeatures); p++) XWy[p] *= (*discount);
				}
			}
			if(*output_regWeight){
				for(int b = prev_id-1; b < batch_id[k]-1; b++){
					for(int p=0; p<*nFeatures; p++) regWeight[C_MAT(b,p,*nBatches)] = beta[p];
				}
			}
			if((*verbose) > 10){
				Rprintf("Data for batch %d:\n", prev_id);
				for(int i=0; i<k; i++){
					Rprintf(" %11.7f  %11.7f", y[i], offset[i]);
					for(int j=0; j<*nFeatures; j++) Rprintf(" %11.7f", X[C_MAT(i,j,*nObs)]);
					Rprintf("\n");
				}
				Rprintf("Regression Weights:\n"); print_vector(NULL, beta, *nFeatures);
			}
		}
		for(int p=0; p<(*nFeatures); p++) x[p] = X[C_MAT(k,p,*nObs)];

		// Make prediction
		double pred = 0;
		for(int p=0; p<(*nFeatures); p++) pred += x[p] * beta[p];
		if((*nOffset) == (*nObs)) pred += offset[k];
		prediction[k] = pred;

		// Update X'WX & X'Wy
		double y_new = y[k];
		if((*nOffset) == (*nObs)) y_new -= offset[k];
		for(int p=0; p<(*nFeatures); p++){
			for(int q=0; q<(*nFeatures); q++){
				XWX[C_MAT(p,q,*nFeatures)] += x[p] * x[q];
			}
			XWy[p] += x[p] * y_new;
		}

		prev_id = batch_id[k];
	}
	if(*output_regWeight){
		if((*nObs) > 0){
			for(int p=0; p<nFeatures_sqr; p++) XWX_temp[p] = XWX[p];
			for(int p=0; p<(*nFeatures); p++) XWy_temp[p] = XWy[p];
			compute_weight_from_XWX_XWy(beta, XWX_temp, XWy_temp, NULL, NULL, nFeatures, &zero, debug);
		}
		for(int b = prev_id-1; b < *nBatches; b++){
			for(int p=0; p<*nFeatures; p++) regWeight[C_MAT(b,p,*nBatches)] = beta[p];
		}
		if((*verbose) > 10){
			Rprintf("------------------------------------------------------\n");
			Rprintf("    Batch %d\n", prev_id);
			Rprintf("------------------------------------------------------\n");
			Rprintf("Data:\n");
			for(int i=0; i<*nObs; i++){
				Rprintf(" %11.7f  %11.7f", y[i], offset[i]);
				for(int j=0; j<*nFeatures; j++) Rprintf(" %11.7f", X[C_MAT(i,j,*nObs)]);
				Rprintf("\n");
			}
			Rprintf("XWX:\n"); print_matrix(NULL,XWX,*nFeatures,*nFeatures);
			Rprintf("XWy:\n"); print_vector(NULL, XWy, *nFeatures);
			Rprintf("Regression Weights:\n"); print_vector(NULL, beta, *nFeatures);
		}
	}

	Free(XWX);
	Free(XWy);
	Free(XWX_temp);
	Free(XWy_temp);
	Free(x);
	Free(beta);
	Free(inv_prior_var);
}

/**
 * - Weight of the examples in the past batches will be discounted by this number
 *   after each batch
 * - Observations must be ordered by itemIndex
 * - Observations for an item after the (nBatches*nObsPerBatch)th are ignored.
 *   They will have batch_id = -1
 */
void perItem_online_factor_batch_predict(
	double *prediction, // nObs x 1 (output)
	int *batch_id,      // nObs x 1 (output, ID starts from 1)
	double *beta,       // nItems x nBatches (output)
	double *v,          // nItems x nBatches x nFactors (output)
	const double *y,    // nObs x 1 (response)
	const int *userIndex,  // nObs x 1 (index start from 1)
	const int *itemIndex,  // nObs x 1 (index start from 1)
	const double *u,       // nUsers x nFactors
	const double *offset,  // nObs x 1
	const double *beta_prior_mean, // nItems x 1
	const double *beta_prior_var,  // nItems x 1
	const double *v_prior_mean, // nItems x nFactors
	const double *v_prior_var,  // nItems x 1 or nItems x nFactors x nFactors
	const double *discount,     // 1x1 (0 <= discount <= 1)
	const int *nObs, const int *nUsers, const int *nItems,
	const int *nFactors, const int *nBatches, const int *nObsPerBatch, const int *v_prior_var_length,
	const int *output_Factors, const int *debug, const int *verbose
){
	int *num = (int*)Calloc(*nItems,int);
	for(int k=0; k<*nObs; k++){
		int item_id = itemIndex[k]-1;
		if(item_id < 0 || item_id >= *nItems) error("Item index out of bound: %d", itemIndex[k]);
		if(num[item_id] < (*nBatches)*(*nObsPerBatch)) num[item_id]++;
		prediction[k] = 0;
		batch_id[k] = -1;
	}
	int nFeatures = (*nFactors)+1;
	double *X = (double*)Calloc((*nBatches)*(*nObsPerBatch)*nFeatures, double);
	double *regWeight = (double*)Calloc((*nBatches)*nFeatures, double);
	double *prior_mean = (double*)Calloc(nFeatures, double);
	double *prior_var  = (double*)Calloc(nFeatures*nFeatures, double);
	int nObs_thisItem = 0;
	int startIndex_thisItem = 0;
	int id_thisItem = -1; // start from 0

	int nProcessed = 0;

	for(int k=0; k<*nObs; k++){
		int item_id = itemIndex[k]-1;
		int user_id = userIndex[k]-1;
		if(item_id < 0 || item_id >= *nItems) error("Item index out of bound: %d", itemIndex[k]);
		if(user_id < 0 || user_id >= *nUsers) error("User index out of bound: %d", userIndex[k]);
		if(id_thisItem > item_id) error("Item index is not sorted!");
		if(id_thisItem < item_id){
			if(nObs_thisItem > 0){
				// Make predictions
				if(nObs_thisItem != num[id_thisItem]) error("internal error!");
				online_gaussian_batch_predict(
					&(prediction[startIndex_thisItem]),
					regWeight,
					&(y[startIndex_thisItem]),
					X,
					&(offset[startIndex_thisItem]),
					&(batch_id[startIndex_thisItem]),
					prior_mean,
					prior_var,
					discount,
					&nObs_thisItem, &nFeatures, &nObs_thisItem, nBatches,
					output_Factors, debug, verbose
				);
				nProcessed++;
				if(*output_Factors){
					for(int b=0; b<*nBatches; b++){
						beta[C_MAT(id_thisItem,b,*nItems)] = regWeight[C_MAT(b,0,*nBatches)];
						for(int f=1; f<nFeatures; f++)
							v[C_3DA(id_thisItem,b,f-1,*nItems,*nBatches)] = regWeight[C_MAT(b,f,*nBatches)];
					}
				}
				if((*verbose) > 0 && nProcessed % 1000 == 0){
					Rprintf("Finished: %d tasks\n", nProcessed);
				}
			}
			// Setup the regression problem for the next item
			id_thisItem = item_id;
			startIndex_thisItem = k;
			nObs_thisItem = 0;
			// Setup prior mean
			prior_mean[0] = beta_prior_mean[id_thisItem];
			for(int p=1; p<nFeatures; p++)
				prior_mean[p] = v_prior_mean[C_MAT(id_thisItem,p-1,*nItems)];
			// Setup prior var
			for(int p=0; p<nFeatures*nFeatures; p++) prior_var[p]=0;
			prior_var[0] = beta_prior_var[id_thisItem];
			if((*v_prior_var_length) == (*nItems)){
				for(int p=1; p<nFeatures; p++) prior_var[C_MAT(p,p,nFeatures)] = v_prior_var[id_thisItem];
			}else if((*v_prior_var_length) == (*nItems)*(*nFactors)*(*nFactors)){
				for(int p=1; p<nFeatures; p++) for(int q=1; q<nFeatures; q++)
					prior_var[C_MAT(p,q,nFeatures)] = v_prior_var[C_3DA(id_thisItem,p-1,q-1,*nItems,*nFactors)];
			}else error("length(v_prior_var) = %d: Should be either %d or %d", *v_prior_var_length, *nItems, (*nItems)*(*nFactors)*(*nFactors));
		}
		if(nObs_thisItem >= (*nBatches)*(*nObsPerBatch)) continue;
		// Setup the regression problem
		X[C_MAT(nObs_thisItem,0,num[item_id])] = 1;
		for(int p=1; p<nFeatures; p++)
			X[C_MAT(nObs_thisItem,p,num[item_id])] = u[C_MAT(user_id,p-1,*nUsers)];
		batch_id[k] = (nObs_thisItem / (*nObsPerBatch)) + 1;
		nObs_thisItem++;
	}
	if(nObs_thisItem > 0){
		// Make predictions
		if(nObs_thisItem != num[id_thisItem]) error("internal error!");
		online_gaussian_batch_predict(
			&(prediction[startIndex_thisItem]),
			regWeight,
			&(y[startIndex_thisItem]),
			X,
			&(offset[startIndex_thisItem]),
			&(batch_id[startIndex_thisItem]),
			prior_mean,
			prior_var,
			discount,
			&nObs_thisItem, &nFeatures, &nObs_thisItem, nBatches,
			output_Factors, debug, verbose
		);
		if(*output_Factors){
			for(int b=0; b<*nBatches; b++){
				beta[C_MAT(id_thisItem,b,*nItems)] = regWeight[C_MAT(b,0,*nBatches)];
				for(int f=1; f<nFeatures; f++)
					v[C_3DA(id_thisItem,b,f-1,*nItems,*nBatches)] = regWeight[C_MAT(b,f,*nBatches)];
			}
		}
	}
	if((*verbose) > 0){
		Rprintf("Finished: %d tasks\n", nProcessed);
	}

	Free(X);
	Free(regWeight);
	Free(prior_mean);
	Free(prior_var);
	Free(num);
}

/**
 * out[k] = t(b) %*% A[k,,] %*% b
 */
void compute_bAb_3DA(
    // Output
    double *out, // nCases x 1
    // Input
    const double *b, // nDim
    const double *A, // nCases x nDim x nDim
    const int *nCases, const int *nDim
){
    for(int k=0; k<(*nCases); k++){
    	out[k] = 0;
    	for(int i=0; i<*nDim; i++){
    		for(int j=0; j<*nDim; j++){
    			out[k] += b[i]*b[j]*A[C_3DA(k,i,j,*nCases,*nDim)];
    		}
    	}
    }
}

/**
 * out[k] = A[k,,] %*% b
 */
void compute_Ab_3DA(
    // Output
    double *out, // nCases x nDim
    // Input
    const double *b, // nDim
    const double *A, // nCases x nDim x nDim
    const int *nCases, const int *nDim
){
    for(int k=0; k<(*nCases); k++){
    	for(int i=0; i<*nDim; i++){
        	double ans = 0;
    		for(int j=0; j<*nDim; j++){
    			ans += A[C_3DA(k,i,j,*nCases,*nDim)]*b[j];
    		}
    		out[C_MAT(k,i,*nCases)] = ans;
    	}
    }
}

/**
 *  output[m] = u[src_id[m], , src_ctx[m]]' v[dst_id[m, , dst_ctx[m]]
 *  All ids/indices start from 1
 */
void computeMultiResponseUV(
	// OUTPUT
	double *output,
	// INPUT
	const double *u, // nSrcNodes x nFactors x nSrcContexts
	const double *v, // nDstNodes x nFactors x nDstContexts
	const int *src_id,  const int *dst_id, // nObs x 1
	const int *src_ctx, const int *dst_ctx, // nObs x 1
	const int *nObs_, const int *nFactors_,
	const int *nSrcNodes_, const int *nSrcContexts_,
	const int *nDstNodes_, const int *nDstContexts_,
	const int *debug_
){
	int nObs = (*nObs_), nFactors = (*nFactors_), debug = (*debug_),
		nSrcNodes = (*nSrcNodes_), nSrcContexts = (*nSrcContexts_),
		nDstNodes = (*nDstNodes_), nDstContexts = (*nDstContexts_);
	for(int m=0; m<nObs; m++){
		int i = src_id[m]-1;
		int j = dst_id[m]-1;
		int k_s = (nSrcContexts==1 ? 0 : src_ctx[m]-1);
		int k_d = (nDstContexts==1 ? 0 : dst_ctx[m]-1);
		if(debug > 0){
			CHK_C_INDEX(i,nSrcNodes); CHK_C_INDEX(k_s,nSrcContexts);
			CHK_C_INDEX(j,nDstNodes); CHK_C_INDEX(k_d,nDstContexts);
		}
		double ans = 0;
		for(int n=0; n<nFactors; n++){
			ans += u[C_3DA(i,n,k_s,nSrcNodes,nFactors)] * v[C_3DA(j,n,k_d,nDstNodes,nFactors)];
		}
		output[m] = ans;
	}
}

/**
 * Draw a multivariate Gaussian vector
 * option == 0 : A    is the variance-covariance matrix
 * option == 1 : A^-1 is the variance-covariance matrix
 * option == 2 : eigen_val and eigen_vec are the eigen value and vector of
 *               the variance-covariance matrix (A is not used)
 * In any case, the output eigen_val and eigen_vec are the eigen value and vector
 * of the variance-covariance matrix.
 *
 * GetRNGstate() and PutRNGstate() are NOT called inside draw_multivar_gaussian().
 */
void draw_multivar_gaussian(
	// OUTPUT
	double* out /* nDim x 1 */,
	// INPUT/OUTPUT
	double* eigen_val /* nDim x 1 */,
	double* eigen_vec /* nDim x nDim */,
	// INPUT
	const double* mean /* nDim x 1*/,
	const double* A    /* nDim x nDim  or  NULL */,
	const int* nDim,
	const int* option, const int* debug,
	// WORK/TEMP SPACE (see workspace_size_for_draw_multivar_gaussian(nDim))
	double *workspace /* workspace_size x 1 */,
	const int* workspace_size /* 0 or at least 3*nDim or call workspace_size_for_draw_multivar_gaussian(nDim) */,
	double *temp1 /* nDim x 1 */,
	double *temp2 /* nDim x nDim */
){
	if((*option) == 0 || (*option) == 1){
		// eigen decomposition of matrix A
	    char jobz = 'V', uplo = 'L';
	    double work_size, *work=NULL;
	    int i, info, lwork=-1;

	    for(i=0; i<(*nDim)*(*nDim); i++) eigen_vec[i] = A[i];

	    if((*workspace_size) == 0){
		    F77_NAME(dsyev)(&jobz, &uplo, nDim, eigen_vec, nDim, eigen_val, &work_size, &lwork, &info);
		    if(info != 0) error("error in dsyev(...)");
		    lwork = (int)work_size;
		    work  = (double*)Calloc(lwork, double);
	    }else{
	    	work = workspace;
	    	lwork = (*workspace_size);
	    }

	    F77_NAME(dsyev)(&jobz, &uplo, nDim, eigen_vec, nDim, eigen_val, work, &lwork, &info);
	    if(info != 0) error("error in dsyev(...)");

	    if((*workspace_size) == 0) Free(work);

	    if((*option) == 1){
	    	for(i=0; i<*nDim; i++) eigen_val[i] = 1/eigen_val[i];
	    }
	}
	// Now, eigen_val and eigen_vec are the eigen value and vector of
	//      the variance-covariance matrix

	// Draw a sample
	int j,k;
	for(j=0; j<*nDim; j++) out[j] = sqrt(eigen_val[j]);
	for(j=0; j<*nDim; j++)
		for(k=0; k<*nDim; k++) temp2[C_MAT(j,k,*nDim)] = eigen_vec[C_MAT(j,k,*nDim)] * out[k];

	for(j=0; j<*nDim; j++) temp1[j] = norm_rand();

	for(j=0; j<*nDim; j++){
		out[j] = mean[j];
		for(k=0; k<*nDim; k++) out[j] += temp2[C_MAT(j,k,*nDim)] * temp1[k];
	}
}
int workspace_size_for_draw_multivar_gaussian(int nDim){
    char jobz = 'V', uplo = 'L';
    double work_size;
    int info, lwork=-1;
    F77_NAME(dsyev)(&jobz, &uplo, &nDim, NULL, &nDim, NULL, &work_size, &lwork, &info);
    if(info != 0) error("error in dsyev(...)");
    lwork = (int)work_size;
    return lwork;
}

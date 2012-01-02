/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include <R_ext/Applic.h>

#define R_VEC(i) (i-1)
#define R_MAT(i,j,nrow) (((j-1)*(nrow))+(i-1))
#define R_3DA(i,j,k,nrow,ncol) ((((k-1)*(ncol))+(j-1))*(nrow) + (i-1))
#define R_4DA(i,j,k,m,dim1,dim2,dim3) (((((m-1)*(dim3)+(k-1))*(dim2))+(j-1))*(dim1) + (i-1))

#define C_MAT(i,j,nrow) (((j)*(nrow))+(i))
#define C_3DA(i,j,k,nrow,ncol) ((((k)*(ncol))+(j))*(nrow) + (i))
#define C_4DA(i,j,k,m,dim1,dim2,dim3) (((((m)*(dim3)+(k))*(dim2))+(j))*(dim1) + (i))

#define CHK_C_INDEX(index, num) if(index < 0 || index >= num) error("index out of bound: index=%d, bound=%d (file: %s, line: %d)", index, num, __FILE__, __LINE__)
#define CHK_C_INDEX_2D(row_i, col_j, nrow, ncol) if(row_i < 0 || row_i >= nrow || col_j < 0 || col_j >= ncol) error("index out of bound: i=%d, j=%d nrow=%d ncol=%d (file: %s, line: %d)", row_i, col_j, nrow, ncol, __FILE__, __LINE__)

#define CHK_R_INDEX(index, num) if(index < 1 || index > num) error("index out of bound: index=%d, bound=%d (file: %s, line: %d)", index, num, __FILE__, __LINE__)
#define CHK_SYMMETRIC(A, nrow) for(int iii=0; iii<nrow; iii++) for(int jjj=0; jjj<iii; jjj++) if(fabs(A[C_MAT(iii,jjj,nrow)] - A[C_MAT(jjj,iii,nrow)]) > 1e-10) error("A symmetric matrix is not symmetric: %f vs %f, diff=%e (file: %s, line: %d)", A[C_MAT(iii,jjj,nrow)], A[C_MAT(jjj,iii,nrow)], A[C_MAT(iii,jjj,nrow)] - A[C_MAT(jjj,iii,nrow)], __FILE__, __LINE__)

#define CHK_SAME_NUMBER(msg, x, y) if((x != 0 || y != 0) && (fabs(x-y) / fmax(fabs(x), fabs(y)) > 1e-8)) error("Error: %s. The two number should be the same: %f vs %f (file: %s, line: %d)", msg, x, y, __FILE__, __LINE__)



#define CHK_MAT_SYM(msg, A, nrow) for(int CHK_MAT_SYM_i=0; CHK_MAT_SYM_i<nrow; CHK_MAT_SYM_i++) for(int CHK_MAT_SYM_j=0; CHK_MAT_SYM_j<CHK_MAT_SYM_i; CHK_MAT_SYM_j++) CHK_SAME_NUMBER(msg, A[C_MAT(CHK_MAT_SYM_i,CHK_MAT_SYM_j,nrow)], A[C_MAT(CHK_MAT_SYM_j,CHK_MAT_SYM_i,nrow)])

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define SQR(x) ((x) * (x))

#define STOP(msg) error("Error at %s on line %d: %s\n", __FILE__, __LINE__, msg)
#define STOP1(msg,x) error2(__FILE__, __LINE__, msg, x)
#define STOP2(msg,x1,x2) error2(__FILE__, __LINE__, msg, x1, x2)
#define STOP3(msg,x1,x2,x3) error2(__FILE__, __LINE__, msg, x1, x2, x3)
#define STOP4(msg,x1,x2,x3,x4) error2(__FILE__, __LINE__, msg, x1, x2, x3, x4)
#define DIE_HERE error("Error in file: %s, at line: %d", __FILE__, __LINE__)

void error2(const char *filename, int lineno, const char *fmt, ...);

// Matrix inversion (only for a symmetric matrix)
void sym_inv_byCholesky(double *A /* n x n matrix */, const int *n, const int *check_sym);

// Matrix inversion (only for a symmetric matrix) 3-dim array
void sym_inv3DA_byCholesky(
    double *invA /* k x n x n matrix */, const double *A /* k x n x n matrix */,
    const int *k, const int *n, double *temp /* n x n */, const int *check_sym
);

// Compute the eigenvalues and vectors of a symmetric matrix x
void sym_eigen(const double* x, const int *nrow, double *eigen_val, double *eigen_vec);

// Compute the eigenvalues and vectors of a symmetric matrix x
void sym_eigen2(const double* x, const int *nrow, double *eigen_val, double *eigen_vec, double *workspace, const int *workspace_size, const int *check_sym);
int sym_eigen2_workspace_size(const int *nrow);

void sum_margin(
    // OUTPUT
    double *ans,
    // INPUT
    const double *A, const int *nrow, const int *ncol, 
    const int *side // side=1: Sum up each row and return a vector with length nrow
                    // side=2: Sum up each column and return a vector with length ncol
);

// The intput/output are all R indices (start from 1, NOT 0)
void generateObsIndex(
    //OUTPUT
    int* obsIndex, // E.g., consider user i (i starts from 0)
    int* start,    //  y[ obsIndex[ start[i]+(0:(num[i]-1)) ] ]
    int* num,      //  are the observations of user i
    //INPUT
    const int* effIndex /* user or item */,
    const int* nObs, const int* nEff,
    //OTHER
    const int* debug
);


void normalizeToSumUpToOne2(double *output, const double *input, const int length);

void normalizeToSumUpToOne(double *vector, const int length);

void indexWithQuantities(int *output, int *vector, const int *length);

void print_vector(const char* prefix, const double* vector, const int length);
void print_intVector(const char* prefix, const int* vector, const int length);
void print_matrix(const char* prefix, const double* matrix, const int nrow, const int ncol);

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
);

double compute_szuBv_single_dense(
	int i, // user i (start from 1)
	int j, // item j (start from 1)
	const double *s, // nUsers x nUserClusters
	const double *z, // nItems x nItemClusters
    const double *u, // nUsers x nFactorsPerUser
	const double *B, // nUserClusters x nItemClusters x nFeaturesPerUser x nFactorsPerItem
    const double *v, // nItems x nFactorsPerItem
    const int nUsers, const int nItems,
    const int nUserClusters, const int nItemClusters,
    const int nFactorsPerUser, const int nFactorsPerItem
);

void compute_szuBv_dense(
	double *ans, // nObs x 1
	const int* userIndex, // nObs x 1 (start from 1)
	const int* itemIndex, // nObs x 1 (start from 1)
	const double *s, // nUsers x nUserClusters
	const double *z, // nItems x nItemClusters
	const double *u, // nUsers x nFactorsPerUser
	const double *B, // nUserClusters x nItemClusters x nFeaturesPerUser x nFactorsPerItem
	const double *v, // nItems x nFactorsPerItem
	const int *nObs, const int *nUsers, const int *nItems,
	const int *nUserClusters, const int *nItemClusters,
	const int *nFactorsPerUser, const int *nFactorsPerItem,
	const int *debug
);

int rdiscrete(double* probabilities, int nOutcomes, double *rnd_out);

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
);

void computeMeanSumvar(
    // OUTPUT
    double *mean   /* length x 1 */,
    double *sumvar /* 1x1 */,
    // INPUT
    const double *sum /* length x 1 */,
    const double *sos /* sum of squares: length x 1 */,
    const int length, const int nSamples
);

void computeMeanVar(
    // OUTPUT
    double *mean   /* nLevels x nFactors */,
    double *sumvar /* 1x1 */,
    double *outputVar /*IN:sum-of-products; OUT:variance; nLevels x nFactors x nFactors*/,
    // INPUT
    const double *sum /* nLevels x nFactors */,
    const int nLevels, const int nFactors, const int nSamples
);

void computeCov(
    // OUTPUT
    double *outputCov /*IN:sum-of-products; OUT:covariance; size: num1 x num2*/,
    // INPUT
    const double *mean1, const double *mean2,
    const int num1, const int num2, const int nSamples
);

/**
 *  output[m] = u[src_id[m], , src_ctx[m]]' v[dst_id[m, , dst_ctx[m]]
 */
void computeMultiResponseUV(
	// OUTPUT
	double *output,
	// INPUT
	const double *u, // nSrcNodes x nFactors x nSrcContexts
	const double *v, // nDstNodes x nFactors x nDstContexts
	const int *src_id,  const int *dst_id, // nObs x 1
	const int *src_ctx, const int *dst_ctx, // nObs x 1
	const int *nObs, const int *nFactors,
	const int *nSrcNodes, const int *nSrcContexts,
	const int *nDstNodes, const int *nDstContexts,
	const int *debug
);


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
);
int workspace_size_for_draw_multivar_gaussian(int nDim);

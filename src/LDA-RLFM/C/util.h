/*
	Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <R_ext/Applic.h>

#define R_VEC(i) (i-1)
#define R_MAT(i,j,nrow) (((j-1)*(nrow))+(i-1))
#define R_3DA(i,j,k,nrow,ncol) ((((k-1)*(ncol))+(j-1))*(nrow) + (i-1))

#define C_MAT(i,j,nrow) (((j)*(nrow))+(i))
#define C_3DA(i,j,k,nrow,ncol) ((((k)*(ncol))+(j))*(nrow) + (i))

#define CHK_C_INDEX(index, num) if(index < 0 || index >= num) error("index out of bound: index=%d, bound=%d (file: %s, line: %d)", index, num, __FILE__, __LINE__)
#define CHK_C_INDEX_2D(row_i, col_j, nrow, ncol) if(row_i < 0 || row_i >= nrow || col_j < 0 || col_j >= ncol) error("index out of bound: i=%d, j=%d nrow=%d ncol=%d (file: %s, line: %d)", row_i, col_j, nrow, ncol, __FILE__, __LINE__)

#define CHK_R_INDEX(index, num) if(index < 1 || index > num) error("index out of bound: index=%d, bound=%d (file: %s, line: %d)", index, num, __FILE__, __LINE__)
#define CHK_SYMMETRIC(A, nrow, i, j) for(i=0; i<nrow; i++) for(j=0; j<i; j++) if(abs(A[C_MAT(i,j,nrow)] - A[C_MAT(j,i,nrow)]) > 1e-10) error("A symmetric matrix is not symmetric: %f vs %f, diff=%e (file: %s, line: %d)", A[C_MAT(i,j,nrow)], A[C_MAT(j,i,nrow)], A[C_MAT(i,j,nrow)] - A[C_MAT(j,i,nrow)], __FILE__, __LINE__)

#define CHK_SAME_NUMBER(msg, x, y) if((x != 0 || y != 0) && (abs(x-y) / fmax(abs(x), abs(y)) > 1e-8)) error("Error: %s. The two number should be the same: %f vs %f (file: %s, line: %d)", msg, x, y, __FILE__, __LINE__)

#define CHK_MAT_SYM(msg, A, nrow) for(int CHK_MAT_SYM_i=0; CHK_MAT_SYM_i<nrow; CHK_MAT_SYM_i++) for(int CHK_MAT_SYM_j=0; CHK_MAT_SYM_j<CHK_MAT_SYM_i; CHK_MAT_SYM_j++) CHK_SAME_NUMBER(msg, A[C_MAT(CHK_MAT_SYM_i,CHK_MAT_SYM_j,nrow)], A[C_MAT(CHK_MAT_SYM_j,CHK_MAT_SYM_i,nrow)])

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define SQR(x) ((x) * (x))

#define DIE_HERE error("Error in file: %s, at line: %d", __FILE__, __LINE__)

typedef struct {
    int* user/*nObs x 1*/;           int* item/*nObs x 1*/;
    double* y/*nObs x 1*/;           double* xb/*nObs x 1*/; 
    double* g0w/*nUsers x 1*/;       double* d0z/*nItems x 1*/;
    double* Gw/*nUsers x nFactors*/; double* Dz/*nItems x nFactors*/; 
    double* var_y; double* var_alpha; double* var_beta; double* var_u; double* var_v;
    int* nObs; int* nUsers; int* nItems; int* nFactors;
    int* User_obsIndex; int* User_oiStart; int* User_oiNum;
    int* Item_obsIndex; int* Item_oiStart; int* Item_oiNum;
    int* debug;
} RegParams;

typedef struct {
    int* user/*nObs x 1*/;           int* item/*nObs x 1*/;
    double* y/*nObs x 1*/;           double* xb/*nObs x 1*/; 
    double* g0w/*nUsers x 1*/;       double* d0z/*nItems x 1*/;
    double* Gw/*nUsers x nFactors*/; double* Dz/*nItems x nFactors*/; 
    double* var_y; double* var_alpha; double* var_beta; double* var_u; double* var_v;
    int* nObs; int* nUsers; int* nItems; int* nFactors;
    int* nVar_y; int* nVar_alpha; int* nVar_beta; int* nVar_u; int* nVar_v;
    double* inv_var_u; double* inv_var_v;
    int* User_obsIndex; int* User_oiStart; int* User_oiNum;
    int* Item_obsIndex; int* Item_oiStart; int* Item_oiNum;
    int* debug;
} RegParams2;

typedef struct {
    int* user/*nObs x 1*/;               int* item/*nObs x 1*/;
    double* y/*nObs x 1*/;               double* x/*nObs x nJointFeatures*/; 
    double* w/*nUsers x nUserFeatures*/; double* z/*nItems x nItemFeatures*/;
    double* lambda_b; double* lambda_g0; double* lambda_d0; double* lambda_G; double* lambda_D;
    int* nObs; int* nUsers; int* nItems; int* nFactors;
    int* nJointFeatures; int* nUserFeatures; int* nItemFeatures;
    int* debug;
} RegInputFeatureOnly;

/*
typedef struct{
    const int* nObs;     int* corpusSize; int* nUsers;         int* nItems;     int* nTerms;
    int* nFactors; int* nTopics;    int* nCorpusWeights;
    int* nVar_y;   int* nVar_alpha; int* nVar_beta;      int* nVar_gamma;
    int* nVar_u;   int* nVar_v;     int* nVar_s;
    int* nEtaCandidates; int* nLambdaCandidates;
} LDA_RLFM_Dim;
// setup LDA_RLFM_DIM
LDA_RLFM_Dim get_LDA_RLFM_Dim(const int *dim, int num);
*/

// Matrix inversion (only for a symmetric matrix)
void sym_inv_byCholesky(double *A /* n x n matrix */, const int *n, const int *check_sym);

// Matrix inversion (only for a symmetric matrix) 3-dim array
void sym_inv3DA_byCholesky(
    double *invA /* k x n x n matrix */, const double *A /* k x n x n matrix */,
    const int *k, const int *n, double *temp /* n x n */, const int *check_sym
);

// Compute the eigenvalues and vectors of a symmetric matrix x
void sym_eigen(const double* x, const int *nrow, double *eigen_val, double *eigen_vec);

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
    int* obsIndex, // E.g., consider user i
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

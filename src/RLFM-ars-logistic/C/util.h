/*
        Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
        Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Liang Zhang
*/

#include <R_ext/Applic.h>

#define R_VEC(i) (i-1)
#define R_MAT(i,j,nrow) (((j-1)*(nrow))+(i-1))
#define R_3DA(i,j,k,nrow,ncol) ((((k-1)*(ncol))+(j-1))*(nrow) + (i-1))

#define C_MAT(i,j,nrow) (((j)*(nrow))+(i))
#define C_3DA(i,j,k,nrow,ncol) ((((k)*(ncol))+(j))*(nrow) + (i))

#define CHK_C_INDEX(index, num) if(index < 0 || index >= num) error("index out of bound: index=%d, bound=%d", index, num)
#define CHK_C_INDEX_2D(row_i, col_j, nrow, ncol) if(row_i < 0 || row_i >= nrow || col_j < 0 || col_j >= ncol) error("index out of bound: i=%d, j=%d nrow=%d ncol=%d", row_i, col_j, nrow, ncol)

#define CHK_R_INDEX(index, num) if(index < 1 || index > num) error("index out of bound: index=%d, bound=%d", index, num)
#define CHK_SYMMETRIC(A, nrow, i, j) for(i=0; i<nrow; i++) for(j=0; j<i; j++) if(abs(A[C_MAT(i,j,nrow)] - A[C_MAT(j,i,nrow)]) > 1e-10) error("A symmetric matrix is not symmetric: %f vs %f, diff=%e", A[C_MAT(i,j,nrow)], A[C_MAT(j,i,nrow)], A[C_MAT(i,j,nrow)] - A[C_MAT(j,i,nrow)])

#define MAX(x,y) (x) > (y) ? (x) : (y)

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

void center_array(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n // array length
);

void center_array_2d(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim // 1 for subtracting row mean, 2 for column
);

void center_array_online(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n, // array length
		  int * oiNum
);

void center_array_2d_online(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim, // 1 for subtracting row mean, 2 for column
		     int * oiNum
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

void print_vector(char* prefix, double* vector, int length);
void print_matrix(const char* prefix, const double* matrix, const int nrow, const int ncol);

void my_cgmin(int n, double *Bvec, double *X, double *Fmin,
	   optimfn fminfn, optimgr fmingr, int *fail,
	   double abstol, double intol, void *ex, int type, int trace,
	   int *fncount, int *grcount, int maxit);

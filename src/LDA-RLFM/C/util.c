/*
	Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include "util.h"

void do_nothing(const void* x){ }

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

void print_vector(const char* prefix, const double* vector, const int length){
    int i;
    if(prefix != NULL) Rprintf("%s", prefix);
    for(i=0; i<length; i++)
        Rprintf("%f ", vector[i]);
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


// setup LDA_RLFM_DIM
/*
LDA_RLFM_Dim get_LDA_RLFM_Dim(const int *dim, int num){
    if(num != 17) error("get_LDA_RLFM_Dim has an error: (num=%d)",num);
    LDA_RLFM_Dim out;
    out.nObs      = dim + 0;       out.corpusSize       = dim + 1;
    out.nUsers    = &(dim[2]);       out.nItems         = &(dim[3]);
    out.nTerms    = &(dim[4]);       out.nFactors       = &(dim[5]);
    out.nTopics   = &(dim[6]);       out.nCorpusWeights = &(dim[7]);
    out.nVar_y    = &(dim[8]);       out.nVar_alpha     = &(dim[9]);
    out.nVar_beta = &(dim[10]);      out.nVar_gamma     = &(dim[11]);
    out.nVar_u    = &(dim[12]);      out.nVar_v         = &(dim[13]);
    out.nVar_s    = &(dim[14]);      out.nEtaCandidates = &(dim[15]);
    out.nLambdaCandidates = &(dim[16]);
    return out;
}
*/

void printPointer(double *pointer){
    Rprintf("Address: %p\n",pointer);
}


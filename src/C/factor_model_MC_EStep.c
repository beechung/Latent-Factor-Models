/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#include <time.h>
#include "util.h"
#include "factor_model_util.h"

// ----------------------------------------------------------------------------
//                              MCEM_EStep
// ----------------------------------------------------------------------------
//  Notation: fErr_{ij} = y_{ij} - (alpha_i + beta_j + v_i' v_j)
//
//  {alpha,beta,v,fErr}_mean       are the Monte-Carlo means of alpha, beta, ...
//  {alpha,beta,v,fErr}_sumvar are the sums of the Monte-Carlo variances over all alpha's, ...
//
//  if outputFactorVar == 1, then
//      {alpha,beta,v}_outputVar will contain the Monte-Carlo variance (cov matrix) for each individual user
//  otherwise (outputFactorVar == 0), {alpha,beta,v}_outputVar will not be changed
//
//  nVar_{y,alpha,...} specifies the length of input var_{y,alpha,...}
//
//  SET nVar_{alpha,beta,...} = 0 to fix the factor values (i.e., prior variance = 0)
// ----------------------------------------------------------------------------
void MCEM_EStep(
    // INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
    double* alpha_mean/*nUsers x 1*/,     double* beta_mean/*nUsers x 1*/,
    double* v_mean/*nUsers x nFactors*/,
    // OUTPUT
    double* alpha_sumvar/*1x1*/,   double* alpha_outputVar/*nUsers x 1*/,
    double* beta_sumvar/*1x1*/,    double* beta_outputVar/*nUsers x 1*/,
    double* v_sumvar/*1x1*/,       double* v_outputVar/*nUsers x nFactors x nFactors*/,
    double* fErr_mean/*nObs x 1*/, double* fErr_sumvar/*1x1*/,
    double* y_pred_square/*nObs x 1: E[(predicted y)^2] for logistic*/,
    // INPUT
    const int* nSamples,                    const int* nBurnIn,
    const int* fromIndex/*nObs x 1*/,       const int* toIndex/*nObs x 1*/,
    const double* y/*nObs x 1*/,            const double* xb/*nObs x 1*/,
    const double* g0x_user/*nUsers x 1*/,   const double* d0x_user/*nUsers x 1*/,
    const double* Gx_user/*nUsers x nFactors*/,
    const double* var_y, const double* var_alpha, const double* var_beta, const double* var_v,
    const int* dim /*7 x 1*/,      const int* nDim /*must be 7*/,
    //  dim = {nObs, nUsers, nFactors,
    //         nVar_y, nVar_alpha, nVar_beta, nVar_v}
    const int* outputFactorVar,
    // OTHER
    const int* debug,  const int* verbose
){
    int user_i, user_j;
    const int one = 1;
    int verbose_nextLevel = (*verbose) - 1;
    clock_t t_begin_in=0;

    if(*verbose > 0) Rprintf("START MCEM_EStep.C\n");

    if(*nDim != 7) error("nDim should be 7: nDim=%d)",*nDim);
    const int* nObs = dim+0;  const int* nUsers = dim+1; const int* nFactors = dim+2;
    const int* nVar_y = dim+3;    const int* nVar_alpha = dim+4; const int* nVar_beta = dim+5;
    const int* nVar_v = dim+6;

    if(*verbose > 5){
    	Rprintf("  nObs=%d, nUsers=%d, nFactors=%d\n",*nObs,*nUsers,*nFactors);
    	Rprintf("  nVar_y=%d, nVar_alpha=%d, nVar_beta=%d, nVar_v=%d\n",*nVar_y,*nVar_alpha,*nVar_beta,*nVar_v);
    }

    // Allocate space for sum and sum-of-squares (or sum of products of a pair of factors)
    double *alpha_sum = (double*)Calloc(*nUsers, double);
    double *beta_sum  = (double*)Calloc(*nUsers, double);
    double *v_sum     = (double*)Calloc((*nUsers)*(*nFactors), double);
    double *fErr_sum  = (double*)Calloc(*nObs, double);
    double *fErr_sos  = (double*)Calloc(*nObs, double);
    double *rest      = (double*)Calloc(*nObs, double);

    for(int k=0; k<*nObs; k++) y_pred_square[k] = 0;

    double *alpha_sos=NULL, *beta_sos=NULL, *v_sos=NULL;
    if((*outputFactorVar) == 0){
        alpha_sos = (double*)Calloc(*nUsers, double);
        beta_sos  = (double*)Calloc(*nUsers, double);
        v_sos     = (double*)Calloc((*nUsers)*(*nFactors), double);
    }else if((*outputFactorVar) == 1){
        // use alpha_outputVar, beta_outputVar and v_outputVar
        for(int k=0; k<*nUsers; k++) alpha_outputVar[k] = 0;
        for(int k=0; k<*nUsers; k++) beta_outputVar[k] = 0;
        for(int k=0; k<(*nUsers)*(*nFactors)*(*nFactors); k++) v_outputVar[k] = 0;
    }else error("outputFactorVar = %d should be only 0 or 1", *outputFactorVar);

    // Allocate space for the observation indices
    int *author_obsIndex = (int*)Calloc(*nObs, int);  int *author_oiStart = (int*)Calloc(*nUsers,int);  int *author_oiNum = (int*)Calloc(*nUsers,int);
    int  *voter_obsIndex = (int*)Calloc(*nObs, int);  int  *voter_oiStart = (int*)Calloc(*nUsers,int);  int  *voter_oiNum = (int*)Calloc(*nUsers,int);


    // Create Observation indices for authors and voters
    generateObsIndex(author_obsIndex, author_oiStart, author_oiNum,   toIndex, nObs, nUsers, debug);
    generateObsIndex( voter_obsIndex,  voter_oiStart,  voter_oiNum, fromIndex, nObs, nUsers, debug);

    // Initialize factors
    // use the memory space of the output to store the current alpha, beta and v
    double *alpha = alpha_mean; double *beta = beta_mean; double *v = v_mean;

    for(int sampleNo=0; sampleNo<(*nSamples)+(*nBurnIn); sampleNo++){

        if(*verbose > 0){
            t_begin_in = clock();
            Rprintf("SAMPLE %3d:",sampleNo);
            if(sampleNo < *nBurnIn) Rprintf(" Burn-in");
            else                    Rprintf(" Ready  ");
        }
        //----------------------------------------
        // Draw samples
        //----------------------------------------
        if(*nVar_alpha > 0){
            // Compute y - xb - beta - vv
            for(int k=0; k<*nObs; k++){
                user_i = toIndex[k] - 1; user_j = fromIndex[k] - 1;
                if(*debug > 0){ CHK_C_INDEX(user_i, *nUsers); CHK_C_INDEX(user_j, *nUsers);}
                double vv = 0;  for(int f=0; f<*nFactors; f++) vv += v[C_MAT(user_i,f,*nUsers)] * v[C_MAT(user_j,f,*nUsers)];
                rest[k] = y[k] - xb[k] - beta[user_j] - vv;
            }
            // print_vector("rest: ", rest, *nObs);

            // Sample alpha
            gaussianPosterior_mainEffect(alpha, NULL, NULL, &one, toIndex, rest, g0x_user, NULL, var_y, var_alpha, nObs, nUsers, nVar_y, nVar_alpha, debug);
        }
        if(*nVar_beta > 0){
            // Compute y - xb - alpha - vv
            for(int k=0; k<*nObs; k++){
                user_i = toIndex[k] - 1; user_j = fromIndex[k] - 1;
                if(*debug > 0){ CHK_C_INDEX(user_i, *nUsers); CHK_C_INDEX(user_j, *nUsers);}
                double vv = 0;  for(int f=0; f<*nFactors; f++) vv += v[C_MAT(user_i,f,*nUsers)] * v[C_MAT(user_j,f,*nUsers)];
                rest[k] = y[k] - xb[k] - alpha[user_i] - vv;
            }
            // Sample beta
            gaussianPosterior_mainEffect(beta, NULL, NULL, &one, fromIndex, rest, d0x_user, NULL, var_y, var_beta, nObs, nUsers, nVar_y, nVar_beta, debug);
        }

        if(*verbose > 0){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf("  draw main: %.1f sec", secUsed);
            t_begin_in = clock();
        }

		// Sample v
		if(*nVar_v > 0){
			// Compute rest = y - xb - alpha - beta
			for(int k=0; k<*nObs; k++){
	            user_i = toIndex[k] - 1; user_j = fromIndex[k] - 1;
	            if(*debug > 0){ CHK_C_INDEX(user_i, *nUsers); CHK_C_INDEX(user_j, *nUsers);}
				rest[k] = y[k] - xb[k] - alpha[user_i] - beta[user_j];
			}

			gaussianPosterior_SelfInteraction(
				v, NULL, NULL, &one, fromIndex, toIndex, rest, Gx_user,
				var_y, var_v, nObs, nUsers, nFactors, nVar_y, nVar_v,
				author_obsIndex, author_oiStart, author_oiNum, voter_obsIndex, voter_oiStart, voter_oiNum, debug
			);
		}

		if(*verbose > 0){
			double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
			Rprintf(" + factor: %.1f sec", secUsed);
			t_begin_in = clock();
		}

        // DEBUG CODE:
        // print_matrix("  s = ", s, *nUsers, *nTopics);

        if(sampleNo < *nBurnIn){
            if(*verbose > 0) Rprintf("\n");
            continue;
        }

        //----------------------------------------
        // Update statistics & output
        //----------------------------------------
        // update {alpha, beta, v}_sum
        for(int k=0; k<(*nUsers); k++) alpha_sum[k] += alpha[k];
        for(int k=0; k<(*nUsers); k++) beta_sum[k]  += beta[k];
        for(int k=0; k<(*nUsers)*(*nFactors); k++) v_sum[k] += v[k];

        // update fErr_sum, fErr_sos
        for(int k=0; k<*nObs; k++){
            user_i = toIndex[k] - 1; user_j = fromIndex[k] - 1;
            if(*debug > 0){ CHK_C_INDEX(user_i, *nUsers); CHK_C_INDEX(user_j, *nUsers);}
            double vv = 0;  for(int f=0; f<*nFactors; f++) vv += v[C_MAT(user_i,f,*nUsers)] * v[C_MAT(user_j,f,*nUsers)];
            double o = y[k] - alpha[user_i] - beta[user_j] - vv;
            double pred_y = alpha[user_i] + beta[user_j] + vv + xb[k];
            y_pred_square[k] += pred_y*pred_y;
            fErr_sum[k] += o;
            fErr_sos[k] += o*o;
        }

        // update sample variances of alpha, beta and v
        if((*outputFactorVar) == 0){
            for(int k=0; k<(*nUsers); k++)          alpha_sos[k] += alpha[k]*alpha[k];
            for(int k=0; k<(*nUsers); k++)           beta_sos[k] +=  beta[k]* beta[k];
            for(int k=0; k<(*nUsers)*(*nFactors); k++)  v_sos[k] +=     v[k]*    v[k];
        }else if((*outputFactorVar) == 1){
            for(int k=0; k<*nUsers; k++) alpha_outputVar[k] += alpha[k]*alpha[k];
            for(int k=0; k<*nUsers; k++)  beta_outputVar[k] +=  beta[k]* beta[k];
            for(int k=0; k<(*nUsers); k++)
                // ONLY store the lower triangle
                for(int f1=0; f1<*nFactors; f1++) for(int f2=0; f2<=f1; f2++)
                    v_outputVar[C_3DA(k,f1,f2,*nUsers,*nFactors)]
                                += v[C_MAT(k,f1,*nUsers)] * v[C_MAT(k,f2,*nUsers)];
        }else error("outputFactorVar = %d should be only 0 or 1", *outputFactorVar);

        if(*verbose > 0){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf(" + update: %.1f sec\n", secUsed);
        }
    }

    if((*outputFactorVar) == 0){
        computeMeanSumvar(alpha_mean, alpha_sumvar, alpha_sum, alpha_sos, *nUsers, *nSamples);
        computeMeanSumvar( beta_mean,  beta_sumvar,  beta_sum,  beta_sos, *nUsers, *nSamples);
        computeMeanSumvar(v_mean, v_sumvar, v_sum, v_sos, (*nUsers)*(*nFactors), *nSamples);
    }else if((*outputFactorVar) == 1){
        computeMeanVar(alpha_mean, alpha_sumvar, alpha_outputVar, alpha_sum, *nUsers, 1, *nSamples);
        computeMeanVar( beta_mean,  beta_sumvar,  beta_outputVar,  beta_sum, *nUsers, 1, *nSamples);
        computeMeanVar(v_mean, v_sumvar, v_outputVar, v_sum, (*nUsers), (*nFactors), *nSamples);
    }else error("outputFactorVar = %d should be only 0 or 1", *outputFactorVar);

    computeMeanSumvar(fErr_mean, fErr_sumvar,  fErr_sum,  fErr_sos, *nObs, *nSamples);
    for(int k=0; k<*nObs; k++) y_pred_square[k] /= (*nSamples);

    // Free the allocated space
    //   The R Free() would only free a pointer if it is NOT NULL.
    if((*outputFactorVar) == 0){ Free(alpha_sos); Free(beta_sos); Free(v_sos); }

    Free(alpha_sum);        Free(beta_sum);         Free(v_sum);
    Free(fErr_sum);          Free(fErr_sos);    	    Free(rest);
    Free(author_obsIndex);  Free(author_oiStart);   Free(author_oiNum);
    Free( voter_obsIndex);  Free( voter_oiStart);   Free( voter_oiNum);

    if(*verbose > 0) Rprintf("END   MCEM_EStep.C\n");
}

/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

/*
To Compile:
R CMD SHLIB C/util.c C/factor_model_util.c C/factor_model_multicontext.c -o C/c_funcs.so
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#include <time.h>
#include "util.h"
#include "factor_model_util.h"

inline double compute_uvw(
	const int m, const double *u, const double *v, const double *w, const int *edgeContext,
	const int from_i, const int to_j, const int nrowU, const int nrowV, const int nEdgeContexts, const int nFactors
){
    double uvw = 0;
    if(nEdgeContexts == 0){
    	if(nrowU == 0)	for(int f=0; f<nFactors; f++) uvw += v[C_MAT(from_i,f,nrowV)] * v[C_MAT(to_j,f,nrowV)];
    	else            for(int f=0; f<nFactors; f++) uvw += u[C_MAT(from_i,f,nrowU)] * v[C_MAT(to_j,f,nrowV)];
    }else{
    	int edge_k = edgeContext[m]-1;
    	CHK_C_INDEX(edge_k, nEdgeContexts);
    	if(nrowU == 0)	for(int f=0; f<nFactors; f++) uvw += v[C_MAT(from_i,f,nrowV)] * v[C_MAT(to_j,f,nrowV)] * w[C_MAT(edge_k,f,nEdgeContexts)];
    	else            for(int f=0; f<nFactors; f++) uvw += u[C_MAT(from_i,f,nrowU)] * v[C_MAT(to_j,f,nrowV)] * w[C_MAT(edge_k,f,nEdgeContexts)];
    }
    return uvw;
}

/*************************************************************************************
 *                              MCEM_EStep
 *************************************************************************************
 *  Notation: fScore[ijk] = alpha[i,k] + beta[j,k] + gamma[k] + sum(v[i,] * v[j,] * w[k,]) if nrowU == 0
 *                                                        ... + sum(u[i,] * v[j,] * w[k,]) if nrowU != 0
 *			  from i to j
 *
 *        obs[ijk] ~ N(fScore[ijk], var_obs[ijk])
 *      alpha[i,k] ~ N(alpha_prior[i,k] + q[k] * alpha_global[i],  var_alpha[i,k])
 * alpha_global[i] ~ N(0,  var_alpha_global[i])
 *       beta[j,k] ~ N(beta_prior[j,k] + r[k] * beta_global[j],  var_beta[j,k])
 *  beta_global[j] ~ N(0,  var_beta_global[j])
 *        gamma[k] ~ N(gamma_prior[k], var_gamma)
 *           u[i,] ~ N(u_prior[i,],  var_u[i,,])
 *           v[j,] ~ N(v_prior[j,],  var_v[j,,])
 *           w[k,] ~ N(w_prior[k,],  var_w[k,,])
 *
 *  {alpha,beta,gamma,u,v,w,fScore}_mean      are the Monte-Carlo means of alpha, ...
 *  {alpha,beta,gamma,u,v,w,fScore}_outputVar are the Monte-Carlo variances of alpha, ...
 *
 *  nVar_{y,alpha,...} specifies the length of input var_{y,alpha,...}
 *
 *  OPTIONS:
 *    * Set nVar_{alpha,beta,...} = 0 to fix the factor values (i.e., prior variance = 0)
 *    * Set nEdgeContexts = 0 to disable w_k
 *    * Set nrowU  = 0 to use <w_k, v_i, v_j>
 *                != 0 to use <w_k, u_i, v_j>
 *    * Set nAlphaContexts = 1 to ignore the hierarchical model.  In this case,
 *      q, alpha_global and var_alpha_global will not be used.
 *      Same for beta.
 *    * Set nGamma = 0 to disable gamma[k]; Otherwise, nGamma must = nEdgeContexts
 *    * Set nTestObs > 0 to enable computation of fScores for test cases.
 *      The i,j,k of all test cases must not be out of bound.
 *
 *  SPECIAL CASE 1:
 *    fScore[ijk] = alpha[i,k] + beta[j,k] + sum(v[i,,k] * v[j,,k])
 *        v[j,,k] ~ N(v_prior[j,,k],  var_v0[k])
 *        v: nrowV x nLocalFactors x nEdgeContexts
 *      * Set nFactors = nLocalFactors * nEdgeContexts
 *      * w = array(0.0, dim=c(nEdgeContexts,nFactors));
 *        for(k in 1:nEdgeContexts) w[k, (k-1)*nLocalFactors + (1:nLocalFactors) ] = 1;
 *      * var_v = rep(0.0, nFactors);
 *        for(k in 1:nEdgeContexts) var_v[ (k-1)*nLocalFactors + (1:nLocalFactors) ] = var_v0[k];
 *      * nVar_w = 0; nVar_v = nFactors;
 */
void MCEM_EStep_multicontext(
    // INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
    double* alpha_mean/*nAlpha x nAlphaContexts*/, double* beta_mean/*nBeta x nBetaContexts*/,
    double* gamma_mean/*nGamma x 1*/,
    double* u_mean/*nrowU x nFactors*/,            double* v_mean/*nrowV x nFactors*/,
    double* w_mean/*nEdgeContexts x nFactors*/,
    // OUTPUT
    double* alpha_global_mean/*nAlpha x 1*/,            double* beta_global_mean/*nBeta x 1*/,
    double* fScore_mean/*nObs x 1*/,                    double* fScore_outputVar/*nObs x 1*/,
    double* alpha_outputVar/*nAlpha x nAlphaContexts*/, double* alpha_global_outputVar/*nAlpha x 1*/,
    double* alpha_outputCov/*nAlpha x nAlphaContexts*/,
    double* beta_outputVar/*nBeta x nBetaContexts*/,    double* beta_global_outputVar/*nBeta x 1*/,
    double* beta_outputCov/*nBeta x nBetaContexts*/,
    double* gamma_outputVar/*nGamma x 1*/,
    double* u_outputVar/*nrowU x nFactors*/,            double* v_outputVar/*nrowV x nFactors*/,
    double* w_outputVar/*nEdgeContexts x nFactors*/,
    double* test_fScore_mean/*nTestObs x 1*/,           double* test_fScore_var/*nTestObs x 1*/,
    // INPUT
    const int* numSamples,                  const int* numBurnIn,
    const int* edgeFrom/*nObs x 1*/,        const int* edgeTo/*nObs x 1*/,
    const int* alphaContext/*nObs x 1*/,    const int* betaContext/*nObs x 1*/,
    const int* edgeContext/*nObs x 1*/,
    const double* q/*nAlphaContexts x 1*/, const double* r/*nBetaContexts x 1*/,
    const double* obs/*nObs x 1*/,
    const double* alpha_prior/*0 or nAlpha or nAlpha x nAlphaContexts*/,
    const double* beta_prior /*0 or nBeta  or  nBeta x nBetaContexts*/,
    const double* gamma_prior/*nGamma x 1*/,
    const double* u_prior/*nrowU x nFactors*/, const double* v_prior/*nrowV x nFactors*/,
    const double* w_prior/*nEdgeContexts x nFactors*/,
    const double* var_obs,  const double* var_alpha,       const double* var_alpha_global,
    const double* var_beta, const double* var_beta_global, const double* var_gamma,
    const double* var_u, const double* var_v,     const double* var_w,
    const int* test_edgeFrom/*nTestObs x 1*/,     const int* test_edgeTo/*nTestObs x 1*/,
    const int* test_alphaContext/*nTestObs x 1*/, const int* test_betaContext/*nTestObs x 1*/,
    const int* test_edgeContext/*nTestObs x 1*/,
    const int* dim /*22 x 1*/,      const int* nDim /*must be 22*/,
    //  dim = {0:nObs, 1:nAlpha, 2:nBeta, 3:nrowU, 4:nrowV, 5:nAlphaContexts, 6:nBetaContexts,
    //         7:nEdgeContexts, 8:nFactors,
    //         9:nVar_y, 10:nVar_alpha, 11:nVar_beta, 12:nVar_u, 13:nVar_v, 14:nVar_w
    //        15:nVar_alpha_global, 16:nVar_beta_global, 17:nAlpha_prior, 18:nBeta_prior,
    //        19:nGamma, 20:nVar_gamma, 21:nTestObs}
    // OTHER
    const int* debug,  const int* verbose
){
    int from_i, to_j, node_k, edge_k, alpha_k, beta_k;
    const int one = 1;
    int verbose_nextLevel = (*verbose) - 1;
    clock_t t_begin_in=0;

    if(*verbose > 0) Rprintf("START MCEM_EStep_multicontext.C\n");

    // Dimensionality
    if(*nDim != 22) error("nDim should be 22: nDim=%d)",*nDim);
    const int nObs = dim[0];  const int nAlpha = dim[1];  const int nBeta = dim[2];  const int nrowU = dim[3];  const int nrowV = dim[4];
    const int nAlphaContexts = dim[5]; const int nBetaContexts = dim[6]; const int nEdgeContexts = dim[7];  const int nFactors = dim[8];
    const int nVar_y = dim[9];  const int nVar_alpha = dim[10]; const int nVar_beta = dim[11];
    const int nVar_u = dim[12]; const int nVar_v = dim[13];  const int nVar_w = dim[14];
    const int nVar_alpha_global = dim[15]; const int nVar_beta_global = dim[16];
    const int nAlpha_prior = dim[17];      const int nBeta_prior = dim[18];
    const int nGamma = dim[19];            const int nVar_gamma = dim[20];
    const int nTestObs = dim[21];
    const int nSamples = numSamples[0];    const int nBurnIn = numBurnIn[0];

    if(nrowU != 0){
    	if(nrowU != nAlpha) error("nrowU != nAlpha");
    	if(nrowV != nBeta)  error("nrowV != nBeta");
    }else if(nFactors != 0){
    	if(nrowV != MAX(nAlpha,nBeta)) error("nrowV != MAX(nAlpha,nBeta)");
    	if(nVar_u != 0) error("nrowU == 0 && nVar_u != 0");
    }
    if(nGamma == 0){
    	if(nVar_gamma != 0) error("nGamma == 0 but nVar_gamma != 0");
    }else if(nGamma != nEdgeContexts) error("nGamma != nEdgeContexts");

    if(*verbose > 5){
    	Rprintf("  nObs=%d, nAlpha=%d, nBeta=%d, nGamma=%d, nrowU=%d, nrowV=%d\n",nObs,nAlpha,nBeta,nGamma,nrowU,nrowV);
    	Rprintf("  nAlphaContexts=%d, nBetaContexts=%d, nEdgeContexts=%d, nFactors=%d\n",nAlphaContexts,nBetaContexts,nEdgeContexts,nFactors);
    	Rprintf("  nVar_y=%d, nVar_alpha=%d, nVar_beta=%d, nVar_gamma=%d, nVar_u=%d, nVar_v=%d, nVar_w=%d\n",nVar_y,nVar_alpha,nVar_beta,nVar_gamma,nVar_u,nVar_v,nVar_w);
    }

    // Allocate space for sum and sum-of-squares (or sum of products of a pair of factors)
    double *alpha_sum  = (double*)Calloc(nAlpha * nAlphaContexts, double);
    double *beta_sum   = (double*)Calloc(nBeta  * nBetaContexts,  double);
    double *gamma_sum  = (double*)Calloc(nGamma,  double);
    double *u_sum      = (double*)Calloc(nrowU  * nFactors, double);
    double *v_sum      = (double*)Calloc(nrowV  * nFactors, double);
    double *w_sum      = (double*)Calloc(nEdgeContexts * nFactors, double);
    double *rest       = (double*)Calloc(nObs, double);

    double *fScore_sum = fScore_mean;       for(int k=0; k<nObs; k++) fScore_sum[k] = 0;
    double *fScore_sos = fScore_outputVar;  for(int k=0; k<nObs; k++) fScore_sos[k] = 0;
    double *alpha_sos = alpha_outputVar;    for(int k=0; k<nAlpha*nAlphaContexts; k++) alpha_sos[k] = 0;
    double *beta_sos  = beta_outputVar;     for(int k=0; k<nBeta *nBetaContexts;  k++) beta_sos[k]  = 0;
    double *gamma_sos = gamma_outputVar;    for(int k=0; k<nGamma;  k++) gamma_sos[k]  = 0;
    double *u_sos = u_outputVar;            for(int k=0; k<nrowU*nFactors; k++) u_sos[k] = 0;
    double *v_sos = v_outputVar;            for(int k=0; k<nrowV*nFactors; k++) v_sos[k] = 0;
    double *w_sos = w_outputVar;            for(int k=0; k<nEdgeContexts*nFactors; k++) w_sos[k] = 0;
    double *test_fScore_sum = test_fScore_mean; for(int k=0; k<nTestObs; k++) test_fScore_sum[k] = 0;
    double *test_fScore_sos = test_fScore_var;  for(int k=0; k<nTestObs; k++) test_fScore_sos[k] = 0;

    double *tempspace_context = (double*)Calloc(MAX(nAlpha*nAlphaContexts, nBeta*nBetaContexts), double);
    double *tempspace_global  = (double*)Calloc(MAX(nAlpha, nBeta), double);
    double *alpha_global_sum=NULL, *alpha_global_sos=NULL, *beta_global_sum=NULL, *beta_global_sos=NULL;
    if(nAlphaContexts > 1){
    	alpha_global_sum = (double*)Calloc(nAlpha, double);
    	alpha_global_sos = alpha_global_outputVar;  for(int k=0; k<nAlpha; k++) alpha_global_sos[k] = 0;
    	for(int k=0; k<nAlpha*nAlphaContexts; k++) alpha_outputCov[k] = 0;
    }
    if(nBetaContexts > 1){
    	beta_global_sum = (double*)Calloc(nBeta, double);
    	beta_global_sos = beta_global_outputVar;  for(int k=0; k<nBeta; k++) beta_global_sos[k] = 0;
    	for(int k=0; k<nBeta *nBetaContexts;  k++) beta_outputCov[k]  = 0;
    }

    // Allocate space for the observation indices
    int *from_obsIndex = (int*)Calloc(nObs, int);  int *from_oiStart = (int*)Calloc(nAlpha,int);  int *from_oiNum = (int*)Calloc(nAlpha,int);
    int *to_obsIndex   = (int*)Calloc(nObs, int);  int *to_oiStart   = (int*)Calloc(nBeta,int);   int *to_oiNum   = (int*)Calloc(nBeta,int);
    int *context_obsIndex = NULL, *context_oiStart = NULL, *context_oiNum = NULL;

    // Create Observation indices for authors and voters
    generateObsIndex(from_obsIndex, from_oiStart, from_oiNum, edgeFrom, &nObs, &nAlpha, debug);
    generateObsIndex(  to_obsIndex,   to_oiStart,   to_oiNum, edgeTo,   &nObs, &nBeta,  debug);
    if(nEdgeContexts > 0){
        context_obsIndex = (int*)Calloc(nObs, int);  context_oiStart = (int*)Calloc(nEdgeContexts,int);  context_oiNum = (int*)Calloc(nEdgeContexts,int);
        generateObsIndex(context_obsIndex, context_oiStart, context_oiNum, edgeContext, &nObs, &nEdgeContexts,  debug);
    }

    // Initialize factors
    // use the memory space of the output to store the current alpha, beta, u, v, w
    double *alpha = alpha_mean, *beta = beta_mean, *gamma = gamma_mean, *u = u_mean, *v = v_mean, *w = w_mean,
    	   *alpha_global = alpha_global_mean, *beta_global = beta_global_mean;

    for(int sampleNo=0; sampleNo<nSamples+nBurnIn; sampleNo++){

        if(*verbose > 1){
            t_begin_in = clock();
            Rprintf("SAMPLE %3d:",sampleNo);
            if(sampleNo < nBurnIn) Rprintf(" Burn-in");
            else                   Rprintf(" Ready  ");
        }

        //----------------------------------------
        // Draw samples
        //----------------------------------------

        // Sample alpha
        if(nVar_alpha > 0){
            // Compute y - beta - gamma - uvw
            for(int m=0; m<nObs; m++){
                from_i = edgeFrom[m]-1;  to_j = edgeTo[m]-1;  node_k = nBetaContexts>1 ? betaContext[m]-1 : 0;
                if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(node_k, nBetaContexts);}
                double uvw = compute_uvw(m,u,v,w,edgeContext,from_i,to_j,nrowU,nrowV,nEdgeContexts,nFactors);
                rest[m] = obs[m] - beta[C_MAT(to_j,node_k,nBeta)] - uvw;
                if(nGamma > 0) rest[m] -= gamma[edgeContext[m]-1];
            }
            if(nAlphaContexts == 1){
				gaussianPosterior_mainEffect(alpha, NULL, NULL, &one, edgeFrom, rest, alpha_prior, NULL, var_obs, var_alpha, &nObs, &nAlpha, &nVar_y, &nVar_alpha, debug);
            }else if(nAlphaContexts > 1){
            	gaussianPosterior_mainEffect_2Levels(
            		alpha, alpha_global, alpha, tempspace_context, alpha_global, tempspace_global,
            		edgeFrom, alphaContext, rest, q, alpha_prior, var_obs, var_alpha, var_alpha_global,
            		&nObs, &nAlpha, &nAlphaContexts, &nVar_y, &nVar_alpha, &nVar_alpha_global, &nAlpha_prior,
            		debug, &verbose_nextLevel);
            }else error("nAlphaContexts = %d", nAlphaContexts);
        }

        // Sample beta
        if(nVar_beta > 0){
            // Compute y - alpha - gamma - uvw
            for(int m=0; m<nObs; m++){
                from_i = edgeFrom[m]-1;  to_j = edgeTo[m]-1;  node_k = nAlphaContexts>1 ? alphaContext[m]-1 : 0;
                if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(node_k, nAlphaContexts);}
                double uvw = compute_uvw(m,u,v,w,edgeContext,from_i,to_j,nrowU,nrowV,nEdgeContexts,nFactors);
                rest[m] = obs[m] - alpha[C_MAT(from_i,node_k,nAlpha)] - uvw;
                if(nGamma > 0) rest[m] -= gamma[edgeContext[m]-1];
            }
            if(nBetaContexts == 1){
            	gaussianPosterior_mainEffect(beta, NULL, NULL, &one, edgeTo, rest, beta_prior, NULL, var_obs, var_beta, &nObs, &nBeta, &nVar_y, &nVar_beta, debug);
            }else if(nBetaContexts > 1){
            	gaussianPosterior_mainEffect_2Levels(
            		beta, beta_global, beta, tempspace_context, beta_global, tempspace_global,
            		edgeTo, betaContext, rest, r, beta_prior, var_obs, var_beta, var_beta_global,
            		&nObs, &nBeta, &nBetaContexts, &nVar_y, &nVar_beta, &nVar_beta_global, &nBeta_prior,
            		debug, &verbose_nextLevel);
            }else error("nBetaContexts = %d", nBetaContexts);
        }

        // Sample gamma
        if(nVar_gamma > 0 && nGamma > 0){
            // Compute y - alpha - beta - uvw
            for(int m=0; m<nObs; m++){
                from_i = edgeFrom[m]-1;  to_j = edgeTo[m]-1;
                int node_ka = nAlphaContexts>1 ? alphaContext[m]-1 : 0;
                int node_kb = nBetaContexts>1  ? betaContext[m]-1  : 0;
                if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(node_ka, nAlphaContexts); CHK_C_INDEX(node_kb, nBetaContexts);}
                double uvw = compute_uvw(m,u,v,w,edgeContext,from_i,to_j,nrowU,nrowV,nEdgeContexts,nFactors);
                rest[m] = obs[m] - alpha[C_MAT(from_i,node_ka,nAlpha)] - beta[C_MAT(to_j,node_kb,nBeta)] - uvw;
            }
			gaussianPosterior_mainEffect(gamma, NULL, NULL, &one, edgeContext, rest, gamma_prior, NULL, var_obs, var_gamma, &nObs, &nGamma, &nVar_y, &nVar_gamma, debug);
        }

        if(*verbose > 1){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf("  draw main: %.1f sec", secUsed);
            t_begin_in = clock();
        }

        if(nVar_u > 0 || nVar_v > 0 || nVar_w > 0){
			// Compute rest = y - alpha - beta - gamma
			for(int m=0; m<nObs; m++){
                from_i = edgeFrom[m]-1;  to_j = edgeTo[m]-1;
                alpha_k = nAlphaContexts>1 ? alphaContext[m]-1 : 0;  beta_k = nBetaContexts>1 ? betaContext[m]-1 : 0;
                if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(alpha_k, nAlphaContexts); CHK_C_INDEX(beta_k, nBetaContexts);}
				rest[m] = obs[m] - alpha[C_MAT(from_i,alpha_k,nAlpha)] - beta[C_MAT(to_j,beta_k,nBeta)];
                if(nGamma > 0) rest[m] -= gamma[edgeContext[m]-1];
			}
        }

        // Sample u
        if(nrowU > 0 && nVar_u > 0 && nFactors > 0){
        	gaussianPosterior_3WayInteraction(
        		u,NULL,NULL,&one,edgeFrom,edgeTo,edgeContext,rest,u_prior,
        		v,w,var_obs,var_u,&nObs,&nrowU,&nrowV,&nEdgeContexts,&nFactors,&nVar_y,&nVar_u,
        		from_obsIndex,from_oiStart,from_oiNum,debug);
        }

		// Sample v
		if(nrowU > 0 && nVar_v > 0 && nFactors > 0){
        	gaussianPosterior_3WayInteraction(
        		v,NULL,NULL,&one,edgeTo,edgeFrom,edgeContext,rest,v_prior,
        		u,w,var_obs,var_v,&nObs,&nrowV,&nrowU,&nEdgeContexts,&nFactors,&nVar_y,&nVar_v,
        		to_obsIndex,to_oiStart,to_oiNum,debug);
		}else if(nVar_v > 0 && nFactors > 0){
			gaussianPosterior_SelfPlusOneInteraction(
				v,NULL,NULL,&one,edgeFrom,edgeTo,edgeContext,rest,v_prior,
				var_obs,var_v,w,&nObs,&nrowV,&nEdgeContexts,&nFactors,&nVar_y,&nVar_v,
				from_obsIndex,from_oiStart,from_oiNum,to_obsIndex,to_oiStart,to_oiNum,debug);
		}

		// Sample w
		if(nEdgeContexts > 0 && nVar_w > 0 && nFactors > 0){
			int num_u = nrowU;
			if(nrowU == 0){ u = v; num_u = nrowV;}
        	gaussianPosterior_3WayInteraction(
        		w,NULL,NULL,&one,edgeContext,edgeFrom,edgeTo,rest,w_prior,
        		u,v,var_obs,var_w,&nObs,&nEdgeContexts,&num_u,&nrowV,&nFactors,&nVar_y,&nVar_w,
        		context_obsIndex,context_oiStart,context_oiNum,debug);
		}

		if(*verbose > 1){
			double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
			Rprintf(" + factor: %.1f sec", secUsed);
			t_begin_in = clock();
		}


        if(sampleNo < nBurnIn){
            if(*verbose > 1) Rprintf("\n");
            continue;
        }

        //----------------------------------------
        // Update statistics & output
        //----------------------------------------
        // update {alpha, beta, u, v, w}_sum
        for(int k=0; k<nAlpha*nAlphaContexts; k++) alpha_sum[k] += alpha[k];
        for(int k=0; k<nBeta *nBetaContexts;  k++) beta_sum[k]  += beta[k];
        for(int k=0; k<nGamma; k++) gamma_sum[k] += gamma[k];
        for(int k=0; k<nrowU*nFactors; k++) u_sum[k] += u[k];
        for(int k=0; k<nrowV*nFactors; k++) v_sum[k] += v[k];
        for(int k=0; k<nEdgeContexts*nFactors; k++) w_sum[k] += w[k];

        // update fScore_sum, fScore_sos
        for(int m=0; m<nObs; m++){
            from_i = edgeFrom[m]-1;  to_j = edgeTo[m]-1;
            alpha_k = nAlphaContexts>1 ? alphaContext[m]-1 : 0;  beta_k = nBetaContexts>1 ? betaContext[m]-1 : 0;
            if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(alpha_k, nAlphaContexts); CHK_C_INDEX(beta_k, nBetaContexts);}
            double uvw = compute_uvw(m,u,v,w,edgeContext,from_i,to_j,nrowU,nrowV,nEdgeContexts,nFactors);
            double o = alpha[C_MAT(from_i,alpha_k,nAlpha)] + beta[C_MAT(to_j,beta_k,nBeta)] + uvw;
            if(nGamma > 0) o += gamma[edgeContext[m]-1];
            fScore_sum[m] += o;
            fScore_sos[m] += o*o;
        }
        // update test_fScore_sum, test_fScore_sos
        for(int m=0; m<nTestObs; m++){
            from_i = test_edgeFrom[m]-1;  to_j = test_edgeTo[m]-1;
            alpha_k = nAlphaContexts>1 ? test_alphaContext[m]-1 : 0;  beta_k = nBetaContexts>1 ? test_betaContext[m]-1 : 0;
            if(*debug > 0){ CHK_C_INDEX(from_i, nAlpha); CHK_C_INDEX(to_j, nBeta); CHK_C_INDEX(alpha_k, nAlphaContexts); CHK_C_INDEX(beta_k, nBetaContexts);}
            double uvw = compute_uvw(m,u,v,w,test_edgeContext,from_i,to_j,nrowU,nrowV,nEdgeContexts,nFactors);
            double o = alpha[C_MAT(from_i,alpha_k,nAlpha)] + beta[C_MAT(to_j,beta_k,nBeta)] + uvw;
            if(nGamma > 0) o += gamma[test_edgeContext[m]-1];
            test_fScore_sum[m] += o;
            test_fScore_sos[m] += o*o;
        }

        // update sum-of-squares of alpha, beta u, v, w
        for(int k=0; k<nAlpha*nAlphaContexts; k++) alpha_sos[k] += SQR(alpha[k]);
        for(int k=0; k<nBeta *nBetaContexts;  k++) beta_sos[k]  += SQR(beta[k]);
        for(int k=0; k<nGamma; k++) gamma_sos[k] += SQR(gamma[k]);
        for(int k=0; k<nrowU*nFactors; k++) u_sos[k] += SQR(u[k]);
        for(int k=0; k<nrowV*nFactors; k++) v_sos[k] += SQR(v[k]);
        for(int k=0; k<nEdgeContexts*nFactors; k++) w_sos[k] += SQR(w[k]);

        if(nAlphaContexts > 1){
        	for(int k=0; k<nAlpha; k++) alpha_global_sum[k] += alpha_global[k];
        	for(int k=0; k<nAlpha; k++) alpha_global_sos[k] += SQR(alpha_global[k]);
        	for(int i=0; i<nAlpha; i++) for(int k=0; k<nAlphaContexts; k++){
        		alpha_outputCov[C_MAT(i,k,nAlpha)] += alpha_global[i] * alpha[C_MAT(i,k,nAlpha)];
        	}
        }
        if(nBetaContexts > 1){
        	for(int k=0; k<nBeta; k++) beta_global_sum[k] += beta_global[k];
        	for(int k=0; k<nBeta; k++) beta_global_sos[k] += SQR(beta_global[k]);
        	for(int i=0; i<nBeta; i++) for(int k=0; k<nBetaContexts; k++){
        		beta_outputCov[C_MAT(i,k,nBeta)] += beta_global[i] * beta[C_MAT(i,k,nBeta)];
        	}
        }

        if(*verbose > 1){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf(" + update: %.1f sec\n", secUsed);
        }
    }

    double temp = 0;
	computeMeanVar( alpha_mean, &temp, alpha_outputVar, alpha_sum, nAlpha*nAlphaContexts, 1, nSamples);
	computeMeanVar(  beta_mean, &temp,  beta_outputVar,  beta_sum,  nBeta*nBetaContexts,  1, nSamples);
	computeMeanVar(     u_mean, &temp,     u_outputVar,     u_sum, nrowU*nFactors, 1, nSamples);
	computeMeanVar(     v_mean, &temp,     v_outputVar,     v_sum, nrowV*nFactors, 1, nSamples);
	computeMeanVar(     w_mean, &temp,     w_outputVar,     w_sum, nEdgeContexts*nFactors, 1, nSamples);
	computeMeanVar(fScore_mean, &temp,fScore_outputVar,fScore_sum, nObs, 1, nSamples);
	if(nTestObs>0) computeMeanVar(test_fScore_mean, &temp, test_fScore_var, test_fScore_sum, nTestObs, 1, nSamples);
	if(nGamma > 0) computeMeanVar( gamma_mean, &temp, gamma_outputVar, gamma_sum, nGamma, 1, nSamples);

	if(nAlphaContexts > 1){
		computeMeanVar(alpha_global_mean, &temp, alpha_global_outputVar, alpha_global_sum, nAlpha, 1, nSamples);
		for(int i=0; i<nAlpha; i++) for(int k=0; k<nAlphaContexts; k++){
			int ik = C_MAT(i,k,nAlpha);
			alpha_outputCov[ik] = alpha_outputCov[ik]/nSamples - alpha_global_mean[i]*alpha_mean[ik];
		}
	}
	if(nBetaContexts > 1){
		computeMeanVar(beta_global_mean, &temp, beta_global_outputVar, beta_global_sum, nBeta, 1, nSamples);
		for(int i=0; i<nBeta; i++) for(int k=0; k<nBetaContexts; k++){
			int ik = C_MAT(i,k,nBeta);
			beta_outputCov[ik] = beta_outputCov[ik]/nSamples - beta_global_mean[i]*beta_mean[ik];
		}
	}

    // Free the allocated space
    //   The R Free() would only free a pointer if it is NOT NULL.

    Free(alpha_sum);  Free(beta_sum);  Free(gamma_sum);  Free(u_sum);  Free(v_sum);  Free(w_sum);
    Free(rest);
    Free(from_obsIndex);     Free(from_oiStart);    Free(from_oiNum);
    Free(to_obsIndex);       Free(to_oiStart);      Free(to_oiNum);
    Free(context_obsIndex);  Free(context_oiStart); Free(context_oiNum);
    Free(tempspace_context); Free(tempspace_global);
    if(nAlphaContexts > 1) Free(alpha_global_sum);
    if(nBetaContexts  > 1) Free( beta_global_sum);

    if(*verbose > 0) Rprintf("END   MCEM_EStep_multicontext.C\n");
}

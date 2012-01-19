/*
	Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

/**
 * Prediction function using Monte-Carlo mean
 *   This implements the transduction leaning: When building the LDA topic model,
 *   assume all the test items are available. Specifically, the test items change the
 *   topic-term counters.
 */
 
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#include <time.h>
#include "util.h"

void condMeanVarSample_singleDim(
    double* outSample, double* outMean, double* outVar,
    const int* option /*1:Sample, 2:Mean&Var, 3:Sample&Mean&Var*/,
    const int* thisEffIndex /*user or item*/, const double* rest /*o in the paper*/,
    const double* fittedEff /*g0w or d0z*/,   const double* multiplier /*NULL or x_dyad*b*/, 
    const double* var_y, const double* var_eff /*var_alpha or var_beta*/,
    const int* nObs, const int* nThisEff, const int* nVar_y, const int* nVar_eff,
    const int* debug
);
void condMeanVarSample_multiDim(
    double* outSample, double* outMean, double* outVar,
    const int* option /*1:Sample, 2:Mean&Var, 3:Sample&Mean&Var*/,
    const int* thisEffIndex /*user or item*/, const int* otherEffIndex /*item or user*/, const double* rest /*o in the paper*/,
    const double* fittedEff /*Gw or Dz*/, const double* otherEff /*v or u*/, 
    const double* var_y, const double* var_eff /*var_u or var_v*/,
    const int* nObs, const int* nThisEff, const int* nOtherEff, const int* nFactors, const int* nVar_y, const int *nVar_eff,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    const int* debug
);
void fillInTopicCounts(
    double *cnt_item_topic, double *cnt_topic_term, double *cnt_topic, double *z_avg,
    const int *corpus_topic, const int *corpus_item, const int *corpus_term, const double *corpus_weight,
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nCorpusWeights,
    const int *debug
);
void fillInTopicCounts2(
    double *cnt_item_topic, double *cnt_topic_term, double *cnt_topic, double *z_avg,
    const double *corpus_topic, const int *corpus_item, const int *corpus_term, const double *corpus_weight,
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nCorpusWeights,
    const int *debug
);
void compute_z_avg(double *z_avg, const double *cnt_item_topic, const int *nItems, const int *nTopics);
void add_to_phi(double *phi, const double *cnt_topic_term, const double *cnt_topic, const double *eta, const int *nTopics, const int *nTerms, const int *debug);
void finalize_phi(double *phi, const int *nTopics, const int *nTerms, const int *nSamples, const int *debug);

void condProbSample_topic(
    int *corpus_topic /*corpusSize x 1*/,
    double *cnt_item_topic /*nItems x nTopics*/,
    double *cnt_topic_term /*nTopics x nTerms*/,
    double *cnt_topic /*nTopics x 1*/,
    double *probDist /*option=1 then NULL; option=2,3 then corpusSize x nTopics; option=4 then nItems x nTopics*/,
    const int* option /*1:Sample, 2:Probabilities, 3:Sample&Probabilities, 4:Sample&itemTopicProb*/,
    const double *rest /*nObs x 1: y - x*b*gamma - alpha - beta - uv */,
    const double *s /*nUsers x nTopics*/, const double *var_y,
    const double *eta, const double *lambda,
    const int *user /*nUsers x 1*/, const int *item /*nItems x 1*/,
    const int *corpus_item /*corpusSize x 1*/, const int *corpus_term /*corpusSize x 1*/, const double *corpus_weight /*corpusSize x 1*/,
    const int *nUsers, const int *nItems, const int *nObs, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nVar_y, const int *nLambda,
    const int *obsIndex, const int *oiStart, const int *oiNum, // observation index for items
    const int *cpsIndex, const int *ciStart, const int *ciNum, // corpus index for items
    const double* nRatingExponent,
    const int *debug,
    const int *verbose,
    const int *checkItemTerm
);
void condProbSample_topic2(
    double *corpus_topic /*corpusSize x nTopics*/,
    double *cnt_item_topic /*nItems x nTopics*/,
    double *cnt_topic_term /*nTopics x nTerms*/,
    double *cnt_topic /*nTopics x 1*/,
    double *probDist /*corpusSize x nTopics or NULL*/,
    const int* option /*1:Sample, 2:Probabilities, 3:Sample&Probabilities*/,
    const double *rest /*nObs x 1: y - x*b*gamma - alpha - beta - uv */,
    const double *s /*nUsers x nTopics*/, const double *var_y,
    const double *eta, const double *lambda,
    const int *user /*nUsers x 1*/, const int *item /*nItems x 1*/,
    const int *corpus_item /*corpusSize x 1*/, const int *corpus_term /*corpusSize x 1*/, const double *corpus_weight /*corpusSize x 1*/,
    const int *nUsers, const int *nItems, const int *nObs, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nVar_y, const int *nLambda,
    const int *obsIndex, const int *oiStart, const int *oiNum, // observation index for items
    const int *cpsIndex, const int *ciStart, const int *ciNum, // corpus index for items
    const double* nRatingExponent,
    const int *debug,
    const int *verbose,
    const int *checkItemTerm
);


// ----------------------------------------------------------------------------
//                              MC_predict
// ----------------------------------------------------------------------------
//
//  {alpha,beta,gamma,u,v,s,z_avg}_mean are the Monte-Carlo means.
//
//  nVar_{y,alpha,...} specifies the length of input var_{y,alpha,...}
//
//  SET nFactors = 0 to disable the u'v part
//  SET nTopics  = 0 to disable the s'z part
//  SET nCorpusWeights = 0 to give each term the same weight
//  SET nVar_{alpha,beta,...} = 0 to fix the factor values (i.e., prior variance = 0)
//
//  SET corpus_topic_type = 0 to use corpus_topic_vector
//                          1 to use corpus_topic_matrix
//
// ----------------------------------------------------------------------------
void MC_predict(
    // INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
    // Required:
    double* alpha_mean/*nUsers x 1*/,    double* beta_mean/*nItems x 1*/,     double* gamma_mean/*nUsers x 1*/,
    // Optional:
    double* u_mean/*nUsers x nFactors*/, double* v_mean/*nItems x nFactors*/, double* s_mean/*nUsers x nTopics*/,
    double* corpus_topic_matrix/*corpusSize x nTopics*/,
    int*    corpus_topic_vector/*corpusSize x 1*/,
    // OUTPUT
    double* z_avg_mean/*nItems x nTopics*/,
    double* prediction/*nTestCases*/,
    double* phi/*nTopics x nTerms: phi[k,] is the term distribution of topic k*/,
    // INPUT
    const int* nSamples,                        const int* nBurnIn,
    const int* user/*nObs x 1*/,                const int* item/*nObs x 1*/, 
    const int* corpus_item/*corpusSize x 1*/,   const int* corpus_term/*corpusSize x 1*/,   const double* input_corpus_weight/*corpusSize x 1 or NULL*/, 
    const double* y/*nObs x 1*/,                const double* xb/*nObs x 1*/, 
    const double* g0x_user/*nUsers x 1*/,       const double* d0x_item/*nItems x 1*/,       const double* c0x_user/*nUsers x 1*/,
    const double* Gx_user/*nUsers x nFactors*/, const double* Dx_item/*nItems x nFactors*/, const double* Hx_user/*nUsers x nTopics*/,
    const double* var_y, const double* var_alpha, const double* var_beta, const double* var_gamma,
    const double* var_u, const double* var_v,     const double* var_s,
    const double* eta,   const double* lambda,
    const int* test_user/*nTestCases x 1*/,  const int* test_item/*nTestCases x 1*/, 
    const double* test_xb/*nTestCases x 1*/,
    const double* nRatingExponent,
    const int* dim /*16 x 1*/,      const int* nDim /*must be 16*/,
    //  dim = {nObs, corpusSize, nUsers, nItems, nTerms, nFactors, nTopics, nCorpusWeights, nVar_y, 
    //         nVar_alpha, nVar_beta, nVar_gamma, nVar_u, nVar_v, nVar_s, nTestCases}
    const int* corpus_topic_type/*0 or 1*/,
    // OTHER
    const int* debug,  const int* verbose
){
    int *obsIndex_user, *oiStart_user, *oiNum_user, 
        *obsIndex_item, *oiStart_item, *oiNum_item,
        *cpsIndex,      *ciStart,      *ciNum,
        option=1, user_i, item_j;
    double *alpha, *beta, *gamma, *u, *v, *s,
           *alpha_sum, *beta_sum, *gamma_sum, 
           *u_sum, *v_sum, *s_sum, *z_avg_sum,
           *rest, *cnt_item_topic, *cnt_topic_term, *cnt_topic, *z_avg;
    const int one = 1;
    int verbose_nextLevel = (*verbose) - 1;
    clock_t t_begin, t_begin_in;
    
    if(*verbose > 0) Rprintf("START MC_predict.C\n");
    
    if(*nDim != 16) error("nDim should be 16: nDim=%d)",nDim);
    const int* nObs = dim+0;  const int* corpusSize = dim+1; const int* nUsers = dim+2; const int* nItems = dim+3; const int* nTerms = dim+4;
    const int* nFactors = dim+5; const int* nTopics = dim+6; const int* nCorpusWeights = dim+7;
    const int* nVar_y = dim+8; const int* nVar_alpha = dim+9; const int* nVar_beta = dim+10; const int* nVar_gamma = dim+11;
    const int* nVar_u = dim+12; const int* nVar_v = dim+13; const int* nVar_s = dim+14;
    const int* nTestCases = dim+15;
    
    // Allocate space for sum and sum-of-squares (or sum of products of a pair of factors)
    alpha_sum      = (double*)Calloc(*nUsers, double);
    beta_sum       = (double*)Calloc(*nItems, double);
    u_sum          = (double*)Calloc((*nUsers)*(*nFactors), double);
    v_sum          = (double*)Calloc((*nItems)*(*nFactors), double);
    s_sum          = (double*)Calloc((*nUsers)*(*nTopics),  double);
    z_avg_sum      = (double*)Calloc((*nItems)*(*nTopics),  double);
    rest           = (double*)Calloc(*nObs, double);
    cnt_item_topic = (double*)Calloc((*nItems)*(*nTopics), double);
    cnt_topic_term = (double*)Calloc((*nTopics)*(*nTerms), double);
    cnt_topic      = (double*)Calloc(*nTopics, double);
    z_avg          = (double*)Calloc((*nItems)*(*nTopics), double);
    
    if(*nVar_gamma > 0){
        gamma_sum = (double*)Calloc(*nUsers, double);
    }
    
    for(int t=0; t<*nTestCases; t++) prediction[t] = 0;
    
    // Allocate space for the observation indices
    obsIndex_user = (int*)Calloc(*nObs, int);  oiStart_user = (int*)Calloc(*nUsers,int);  oiNum_user = (int*)Calloc(*nUsers,int);
    obsIndex_item = (int*)Calloc(*nObs, int);  oiStart_item = (int*)Calloc(*nItems,int);  oiNum_item = (int*)Calloc(*nItems,int);
    cpsIndex = (int*)Calloc(*corpusSize,int);       ciStart = (int*)Calloc(*nItems,int);       ciNum = (int*)Calloc(*nItems,int);

    // Use the memory space of the output to store the current alpha, beta, u and v
    alpha = alpha_mean; beta = beta_mean; gamma = gamma_mean; u = u_mean; v = v_mean; s = s_mean;
    
    // Create Observation indices for users and items
    generateObsIndex(obsIndex_user, oiStart_user, oiNum_user, user,        nObs,       nUsers, debug);
    generateObsIndex(obsIndex_item, oiStart_item, oiNum_item, item,        nObs,       nItems, debug);
    if(*nTopics > 0) generateObsIndex(cpsIndex, ciStart, ciNum, corpus_item, corpusSize, nItems, debug);

    const double *corpus_weight = (*nCorpusWeights == 0 ? NULL : input_corpus_weight);
    if(*nCorpusWeights != 0 && *nCorpusWeights != *corpusSize) 
        error("nCorpusWeights = %d, but corpuseSize = %d", *nCorpusWeights, *corpusSize);
    
    // Initialize topic counts & z_avg
    if(*nTopics > 0){
        if(*corpus_topic_type == 0){
            fillInTopicCounts(cnt_item_topic, cnt_topic_term, cnt_topic, z_avg,
                corpus_topic_vector, corpus_item, corpus_term, corpus_weight, 
                nItems, corpusSize, nTopics, nTerms, nCorpusWeights, debug);
        }else if(*corpus_topic_type == 1){
            fillInTopicCounts2(cnt_item_topic, cnt_topic_term, cnt_topic, z_avg,
                corpus_topic_matrix, corpus_item, corpus_term, corpus_weight, 
                nItems, corpusSize, nTopics, nTerms, nCorpusWeights, debug);
        }else error("Unknown corpus_topic_type: %d", *corpus_topic_type);

        for(int k=0; k<(*nTopics)*(*nTerms); k++) phi[k] = 0;
    }
    
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
            // Compute y - (xb)gamma - beta - uv - sz
            for(int k=0; k<*nObs; k++){
                user_i = user[k]; item_j = item[k]; if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
                double uv = 0;  for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
                double sz = 0;  for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
                rest[k] = y[k] - (xb[k])*gamma[R_VEC(user_i)] - beta[R_VEC(item_j)] - uv - sz;
            }
            // print_matrix("z_avg: ", z_avg, *nItems, *nTopics);
            // print_vector("rest: ", rest, *nObs);
            
            // Sample alpha
            condMeanVarSample_singleDim(alpha, NULL, NULL, &option, user, rest, g0x_user, NULL, var_y, var_alpha, nObs, nUsers, nVar_y, nVar_alpha, debug);
        }
        if(*nVar_beta > 0){
            // Compute y - (xb)gamma - alpha - uv - sz
            for(int k=0; k<*nObs; k++){
                user_i = user[k]; item_j = item[k]; if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
                double uv = 0;  for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
                double sz = 0;  for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
                rest[k] = y[k] - (xb[k])*gamma[R_VEC(user_i)] - alpha[R_VEC(user_i)] - uv - sz;
            }
            // Sample beta
            condMeanVarSample_singleDim(beta, NULL, NULL, &option, item, rest, d0x_item, NULL, var_y, var_beta, nObs, nItems, nVar_y, nVar_beta, debug);
        }
        if(*nVar_gamma > 0){
            // Compute y - alpha - beta - uv - sz
            for(int k=0; k<*nObs; k++){
                user_i = user[k]; item_j = item[k]; if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
                double uv = 0;  for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
                double sz = 0;  for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
                rest[k] = y[k] - alpha[R_VEC(user_i)] - beta[R_VEC(item_j)] - uv - sz;
            }
            // Sample gamma
            condMeanVarSample_singleDim(gamma, NULL, NULL, &option, user, rest, c0x_user, xb, var_y, var_gamma, nObs, nUsers, nVar_y, nVar_gamma, debug);
        }

        if(*verbose > 0){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf("  draw main: %.1f sec", secUsed);
            t_begin_in = clock();
        }
        
        if(*nFactors > 0){
            // Compute y - (xb)gamma - alpha - beta - sz
            for(int k=0; k<*nObs; k++){
                user_i = user[k]; item_j = item[k]; if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
                double sz = 0;  for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
                rest[k] = y[k] - (xb[k])*gamma[R_VEC(user_i)] - alpha[R_VEC(user_i)] - beta[R_VEC(item_j)] - sz;
            }
            // Sample u
            if(*nVar_u > 0)
                condMeanVarSample_multiDim(u, NULL, NULL, &option, user, item, rest, Gx_user, v, var_y, var_u, nObs, nUsers, nItems, nFactors, nVar_y, nVar_u, obsIndex_user, oiStart_user, oiNum_user, debug);
            // Sample v
            if(*nVar_v > 0)
                condMeanVarSample_multiDim(v, NULL, NULL, &option, item, user, rest, Dx_item, u, var_y, var_v, nObs, nItems, nUsers, nFactors, nVar_y, nVar_v, obsIndex_item, oiStart_item, oiNum_item, debug);

            if(*verbose > 0){
                double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
                Rprintf(" + factor: %.1f sec", secUsed);
                t_begin_in = clock();
            }
        }
        
        if(*nTopics > 0){
            // Compute y - (xb)gamma - alpha - beta - uv
            for(int k=0; k<*nObs; k++){
                user_i = user[k]; item_j = item[k]; if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
                double uv = 0;  for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
                rest[k] = y[k] - (xb[k])*gamma[R_VEC(user_i)] - alpha[R_VEC(user_i)] - beta[R_VEC(item_j)] - uv;
            }

            // Sample s
            if(*nVar_s > 0)
                condMeanVarSample_multiDim(s, NULL, NULL, &option, user, item, rest, Hx_user, z_avg, var_y, var_s, nObs, nUsers, nItems, nTopics, nVar_y, nVar_s, obsIndex_user, oiStart_user, oiNum_user, debug);

            if(*corpus_topic_type == 0){
                // Sample topics
                condProbSample_topic(corpus_topic_vector, cnt_item_topic, cnt_topic_term, cnt_topic, z_avg,
                    &option, rest, s, var_y, eta, lambda, user, item, corpus_item, corpus_term, corpus_weight,
                    nUsers, nItems, nObs, corpusSize, nTopics, nTerms, nVar_y, &one, 
                    obsIndex_item, oiStart_item, oiNum_item, cpsIndex, ciStart, ciNum, nRatingExponent, debug, &verbose_nextLevel, &one
                );
            }else if(*corpus_topic_type == 1){
                // Sample topics
                condProbSample_topic2(corpus_topic_matrix, cnt_item_topic, cnt_topic_term, cnt_topic, NULL,
                    &option, rest, s, var_y, eta, lambda, user, item, corpus_item, corpus_term, corpus_weight,
                    nUsers, nItems, nObs, corpusSize, nTopics, nTerms, nVar_y, &one, 
                    obsIndex_item, oiStart_item, oiNum_item, cpsIndex, ciStart, ciNum, nRatingExponent, debug, &verbose_nextLevel, &one
                );
            }else error("error");
            // Update z_avg
            compute_z_avg(z_avg, cnt_item_topic, nItems, nTopics);
            
            // print_matrix(" cnt_item_topic: ", cnt_item_topic, *nItems, *nTopics);
            
            if(*verbose > 0){
                double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
                Rprintf(" + topic: %.1f sec", secUsed);
                t_begin_in = clock();
            }
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
        // update {alpha, beta, u, v, s, z_avg}_sum
        for(int k=0; k<*nUsers; k++)               alpha_sum[k] += alpha[k];
        for(int k=0; k<*nItems; k++)               beta_sum[k]  += beta[k];
        for(int k=0; k<(*nUsers)*(*nFactors); k++) u_sum[k]     += u[k];
        for(int k=0; k<(*nItems)*(*nFactors); k++) v_sum[k]     += v[k];
        for(int k=0; k<(*nUsers)*(*nTopics);  k++) s_sum[k]     += s[k];
        for(int k=0; k<(*nItems)*(*nTopics);  k++) z_avg_sum[k] += z_avg[k];
        // update gamma_sum
        if(*nVar_gamma > 0){
            for(int k=0; k<*nUsers; k++) gamma_sum[k] += gamma[k];
        }
        // update prediction
        for(int t=0; t<*nTestCases; t++){
            user_i = test_user[t]; item_j = test_item[t]; 
            if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems);}
            double uv = 0;  for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
            double sz = 0;  for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
            prediction[t] += test_xb[t] * gamma[R_VEC(user_i)] + alpha[R_VEC(user_i)] + beta[R_VEC(item_j)] + uv + sz;
        }
        // update phi
        if(*nTopics > 0)
            add_to_phi(phi, cnt_topic_term, cnt_topic, eta, nTopics, nTerms, debug);

        if(*verbose > 0){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf(" + update: %.1f sec\n", secUsed);
        }
    }

    for(int i=0; i<*nUsers; i++) alpha_mean[i] = alpha_sum[i] / (double)(*nSamples);
    for(int j=0; j<*nItems; j++)  beta_mean[j] =  beta_sum[j] / (double)(*nSamples);
    if(*nVar_gamma > 0){
        for(int i=0; i<*nUsers; i++) gamma_mean[i] = gamma_sum[i] / (double)(*nSamples);
    }
    if(*nFactors > 0){
        for(int k=0; k<(*nUsers)*(*nFactors); k++) u_mean[k] = u_sum[k] / (double)(*nSamples);
        for(int k=0; k<(*nItems)*(*nFactors); k++) v_mean[k] = v_sum[k] / (double)(*nSamples);
    }
    if(*nTopics > 0){
        for(int k=0; k<(*nUsers)*(*nTopics); k++) s_mean[k]     = s_sum[k]     / (double)(*nSamples);
        for(int k=0; k<(*nItems)*(*nTopics); k++) z_avg_mean[k] = z_avg_sum[k] / (double)(*nSamples);
    }
    for(int t=0; t<*nTestCases; t++) prediction[t] /= (double)(*nSamples);

    if(*nTopics > 0) finalize_phi(phi, nTopics, nTerms, nSamples, debug);
    
    // Free the allocated space
    //   The R Free() would only free a point if it is NOT NULL.
    if((*nVar_gamma) > 0){ Free(gamma_sum); }
    Free(alpha_sum);      Free(beta_sum);       Free(u_sum);
    Free(v_sum);          Free(s_sum);          Free(z_avg_sum);
    Free(rest);           Free(z_avg);
    Free(obsIndex_user);  Free(oiStart_user);   Free(oiNum_user);
    Free(obsIndex_item);  Free(oiStart_item);   Free(oiNum_item);
    Free(cpsIndex);       Free(ciStart);        Free(ciNum);
    Free(cnt_item_topic); Free(cnt_topic_term); Free(cnt_topic);
    
    if(*verbose > 0) Rprintf("END   MC_predict.C\n");
}

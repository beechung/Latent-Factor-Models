/*
	Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

/*
  Do not draw multinomial sample; instead, use the probabilities directly

To Compile:
R CMD SHLIB C/util.c C/MCEM_EStep.c -o C/MCEM_EStep.so

*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#include <time.h>
#include "util.h"


void validate_corpus_counts2(
    const double *corpus_topic, const double *cnt_item_topic, const double *cnt_topic_term, const double *cnt_topic,
    const int *corpus_item, const int *corpus_term, const double *corpus_weight, 
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms
);
void computeMeanSumvar(
    double *mean, double *sumvar,
    const double *sum, const double *sos, const int length, const int nSamples
);
void computeMeanVar(
    double *mean, double *sumvar, double *outputVar,
    const double *sum, const int nEffects, const int nFactors, const int nSamples
);
void compute_z_avg(double *z_avg, const double *cnt_item_topic, const int *nItems, const int *nTopics);
void add_LDA_prior_objval_part2(
    double *objval, const double *candidates, const int *nCandidates,
    const double *cnt_dim1, const double *cnt_2dim, const int *dim1, const int *dim2
);
void fillInTopicCounts2(
    double *cnt_item_topic, double *cnt_topic_term, double *cnt_topic, double *z_avg,
    const double *corpus_topic, const int *corpus_item, const int *corpus_term, const double *corpus_weight,
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nCorpusWeights,
    const int *debug
);
void add_to_phi(double *phi, const double *cnt_topic_term, const double *cnt_topic, const double *eta, const int *nTopics, const int *nTerms, const int *debug);
void finalize_phi(double *phi, const int *nTopics, const int *nTerms, const int *nSamples, const int *debug);
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
void add_LDA_lambda_objval_part2(
    double *objval,
    const double *candidates, const int *nCandidates,
    const double *cnt_item_topic, const int *nItems, const int *nTopics,
    const double *nRatingExponent, const int *oiNum
);
double compute_LDA_lambda_objval_part1(
    double candidate,
    const int *nItems, const int *nTopics,
    const double *nRatingExponent, const int *oiNum
);


// FUNCTION: condProbSample_topic2
//  
// Observation index (consider, say, item j, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
//  obsIndex[ oiStart[j]+0 ], ..., obsIndex[ oiStart[j]+oiNum[j]-1 ] are the indices of item j's observations
//
// Corpus index (consider item j); they are R indices (starting from 1, NOT 0)
//  cpsIndex[ ciStart[j]+0 ], ..., cpsIndex[ ciStart[j]+ciNum[j]-1 ] are the indices of item j's terms in the corpus
//  Sanity check: corpus_item[cpsIndex[ ciStart[j]+k ]] = j
//
void condProbSample_topic2(
    // INPUT & OUTPUT
    double *corpus_topic /*corpusSize x nTopics*/,
    double *cnt_item_topic /*nItems x nTopics*/,
    double *cnt_topic_term /*nTopics x nTerms*/,
    double *cnt_topic /*nTopics x 1*/,
    // OUTPUT
    double *probDist /*corpusSize x nTopics or NULL*/,
    // INPUT
    const int* option /*1:Sample, 2:Probabilities, 3:Sample&Probabilities*/,
    const double *rest /*nObs x 1: y - x*b*gamma - alpha - beta - uv */,
    const double *s /*nUsers x nTopics*/, const double *var_y,
    const double *eta, const double *lambda,
    const int *user /*nUsers x 1*/, const int *item /*nItems x 1*/,
    const int *corpus_item /*corpusSize x 1*/, const int *corpus_term /*corpusSize x 1*/, const double *corpus_weight /*corpusSize x 1*/,
    const int *nUsers, const int *nItems, const int *nObs, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nVar_y, const int *nLambda,
    const int *obsIndex, const int *oiStart, const int *oiNum, // observation index for items
    const int *cpsIndex, const int *ciStart, const int *ciNum, // corpus index for items
    const double *nRatingExponent,
    // OTHER
    const int *debug,
    const int *verbose,
    const int *checkItemTerm
){
    double *prob, *z_avg_j, *B_j, *C_j, *s_i, o, *newLambda;
        // B_j: nTopics x 1:       sum_{i in I_j} (o_{ij} s_i^') / sigma_{ij}^2
        // C_j: nTopics x nTopics: sum_{i in I_j} (s_i    s_i^') / sigma_{ij}^2
    int outputSample, outputProb, *draw;
    
    if(*verbose > 0) Rprintf("condProbSample_topic2: begin\n");
    
    if(*debug >= 2) validate_corpus_counts2(corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic, corpus_item, corpus_term, corpus_weight, nItems, corpusSize, nTopics, nTerms);
    
    if(*option == 1){
        outputSample = 1; outputProb = 0;
    }else if(*option == 2){
        outputSample = 0; outputProb = 1;
    }else if(*option == 3){
        outputSample = 1; outputProb = 1;
    }else error("Unknown option: %d", *option);

    
    prob    = (double*)Calloc(*nTopics, double);
    draw    = (int*)   Calloc(*nTopics, int);
    z_avg_j = (double*)Calloc(*nTopics, double);
    s_i     = (double*)Calloc(*nTopics, double);
    B_j     = (double*)Calloc(*nTopics, double);
    C_j     = (double*)Calloc((*nTopics)*(*nTopics), double);
    newLambda = (double*)Calloc(*nLambda, double);

    if(outputSample) GetRNGstate();

    for(int j=0; j<*nItems; j++){

        for(int k=0; k<*nLambda; k++){
            newLambda[k] = lambda[k];
            if(*nRatingExponent != 0) newLambda[k] *= R_pow(1 + oiNum[j], *nRatingExponent);
        }
        
        //---------------------------------------
        // Initialize B_j, C_j
        //---------------------------------------
        int itemID = j+1;
        for(int k=0; k<*nTopics; k++) B_j[k] = 0;
        for(int k=0; k<(*nTopics)*(*nTopics); k++) C_j[k] = 0;
        
        for(int i=0; i<oiNum[j]; i++){
            
            int oIndex = obsIndex[R_VEC(oiStart[j]+i)];
            if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
            
            int userID = user[R_VEC(oIndex)];
            if(*debug > 0) CHK_R_INDEX(userID, *nUsers);
            
            if(*debug > 1) if(item[R_VEC(oIndex)] != itemID) error("error in obsIndex, oiStart, oiNum\n");
            
            o = rest[R_VEC(oIndex)];
            for(int k=1; k<=*nTopics; k++) s_i[R_VEC(k)] = s[R_MAT(userID,k,*nUsers)];
            
            double var_y_thisObs = 0;
            if((*nVar_y) == 1)            var_y_thisObs = var_y[0];
            else if((*nVar_y) == (*nObs)) var_y_thisObs = var_y[R_VEC(oIndex)];
            else error("nVar_y = %d, nObs = %d", *nVar_y, *nObs);
            
            for(int k=0; k<*nTopics; k++) B_j[k] += (o * s_i[k]) / var_y_thisObs;
            for(int k=0; k<*nTopics; k++)
                for(int m=0; m<*nTopics; m++) C_j[C_MAT(k,m,*nTopics)] += (s_i[k] * s_i[m]) / var_y_thisObs;
        }
        // NOTE: Now,
        //   B_j = sum_{i in I_j} (o_{ij} s_i) / sigma_{ij}^2    (size: nTopics x 1)
        //   C_j = sum_{i in I_j} (s_i  s_i^') / sigma_{ij}^2    (size: nTopics x nTopics)

        double W_j = 0; // sum of weights for item j
        for(int k=0; k<*nTopics; k++){
            z_avg_j[k] = cnt_item_topic[C_MAT(j,k,*nItems)];
            W_j += z_avg_j[k];
        }
        if(*checkItemTerm && W_j == 0){
            warning("Item %d has no terms", itemID);
            continue;
        }
        normalizeToSumUpToOne(z_avg_j, *nTopics);

        if(*debug > 1) CHK_MAT_SYM("C_j should be symmetric",C_j,*nTopics);
        if(*verbose > 1){
            Rprintf("  Process item: %d\n", j+1);
            if(*verbose > 3){
                Rprintf("    B_j = "); print_vector("", B_j, *nTopics);
                Rprintf("    C_j = \n");
                print_matrix("          ", C_j, *nTopics, *nTopics);
                Rprintf("    z_avg_j = "); print_vector("", z_avg_j, *nTopics);
           }
        }
                
        //---------------------------------------
        // Process each term in this item
        //---------------------------------------
        for(int n=0; n<ciNum[j]; n++){
            
            // Setup local variables
            int cIndex = cpsIndex[R_VEC(ciStart[j]+n)];
            if(*debug > 0) CHK_R_INDEX(cIndex, *corpusSize);
            
            int thisItem   = corpus_item[ R_VEC(cIndex)];
            int thisTerm   = corpus_term[ R_VEC(cIndex)];
            double absWeight = (corpus_weight != NULL ? corpus_weight[R_VEC(cIndex)] : 1);
            double relWeight = absWeight / W_j;

            if(*debug > 0){ CHK_R_INDEX(thisItem, *nItems); CHK_R_INDEX(thisTerm, *nTerms); }
            if(*debug > 1) if(thisItem != itemID) error("error in cpsIndex, ciStart, ciNum\n");
            
            // remove the topic-weight of the current term
            for(int k=1; k<=*nTopics; k++) z_avg_j[R_VEC(k)] -= relWeight * corpus_topic[R_MAT(cIndex,k,*corpusSize)];
            
            double Bjzj_rm = 0;
            double zjCjzj_rm = 0;
            for(int k=0; k<*nTopics; k++) Bjzj_rm += B_j[k] * z_avg_j[k];
            for(int k=0; k<*nTopics; k++) for(int m=0; m<*nTopics; m++) zjCjzj_rm += C_j[C_MAT(k,m,*nTopics)] * z_avg_j[k] * z_avg_j[m];
            
            // Rprintf("      Bjzj_rm   = %f\n", Bjzj_rm);
            // Rprintf("      zjCjzj_rm = %f\n", zjCjzj_rm);
            
            for(int k=1; k<=*nTopics; k++){
                
                // Compute rating likelihood for this item
                double Bjzj   = Bjzj_rm + (B_j[R_VEC(k)] * relWeight);
                double zjCjzj = zjCjzj_rm + 
                                (C_j[R_MAT(k,k,*nTopics)] * relWeight * (2 * z_avg_j[R_VEC(k)] + relWeight));
                for(int m=1; m<=*nTopics; m++) if(m!=k) zjCjzj += 2 * C_j[R_MAT(k,m,*nTopics)] * z_avg_j[R_VEC(m)] * relWeight;
                // Rprintf("      zjCjzj = %f\n", zjCjzj);
                
                // Compute LDA probability
                double thisLambda = newLambda[0];
                if(*nLambda == *nTopics) thisLambda = newLambda[R_VEC(k)];
                else if(*nLambda != 1) error("nLambda != 1 and != nTopics!");
                                
                double lda_prob = 0;
                double temp = absWeight * corpus_topic[R_MAT(cIndex,k,*corpusSize)];
                
                lda_prob = ( (cnt_topic_term[R_MAT(k,thisTerm,*nTopics)] - temp + eta[0]) /
                                (cnt_topic[R_VEC(k)] - temp + (*nTerms) * eta[0]) ) *
                              ( cnt_item_topic[R_MAT(thisItem,k,*nItems)] - temp + thisLambda);
                
                if(oiNum[j] == 0 && (zjCjzj != 0 || Bjzj != 0)) error("Error in %s at line %d", __FILE__, __LINE__);
                if(*debug > 0 && lda_prob <= 0) error("LDA probability < 0");
                
                prob[R_VEC(k)] = log(lda_prob) + (-0.5 * zjCjzj) + Bjzj;
                
                // Rprintf("      %3d: lda_prob=%f\tLL=%f\n", k, lda_prob, (-0.5 * zjCjzj) + Bjzj);
            }
            
            // Now, prob[k] = log prob[k]
            double max_prob = prob[0];
            for(int k=1; k<*nTopics; k++) if(prob[k] > max_prob) max_prob = prob[k];
            for(int k=0; k<*nTopics; k++) prob[k] = exp(prob[k] - max_prob);
            
            normalizeToSumUpToOne(prob, *nTopics);
            if(outputProb){
                for(int k=1; k<=*nTopics; k++) probDist[R_MAT(cIndex,k,*corpusSize)] = prob[R_VEC(k)];
            }
            
            if(*verbose > 2) Rprintf("    Process the %5dth term: %5d (absWeight: %f, relWeight: %f)\n", n+1, thisTerm, absWeight, relWeight);

            // Update output: corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic
            // Update z_avg_j
            for(int k=1; k<=*nTopics; k++){
                double diff = absWeight * (-corpus_topic[R_MAT(cIndex,k,*corpusSize)] + prob[R_VEC(k)]);
                // update corpus_topic
                corpus_topic[R_MAT(cIndex,k,*corpusSize)] = prob[R_VEC(k)];
                // update cnt_item_topic
                cnt_item_topic[R_MAT(itemID, k, *nItems)] += diff;
                // update cnt_topic_term
                cnt_topic_term[R_MAT(k, thisTerm, *nTopics)] += diff;
                // update cnt_topic
                cnt_topic[R_VEC(k)] += diff;
                // update z_avg_j
                z_avg_j[R_VEC(k)] += relWeight * prob[R_VEC(k)];
            }
            
            // Do some sanity checks
            if(*debug >= 5){
                double temp = 0;
                for(int k=0; k<*nTopics; k++) temp += z_avg_j[k];
                CHK_SAME_NUMBER("sum(z_avg_j) != 1", temp, 1);
                
                // Do a very, very expensive sanity check
                if(*debug > 100) validate_corpus_counts2(corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic, corpus_item, corpus_term, corpus_weight, nItems, corpusSize, nTopics, nTerms);
            }
        }
    }
    
    if(outputSample) PutRNGstate();

    if(*debug >= 2) validate_corpus_counts2(corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic, corpus_item, corpus_term, corpus_weight, nItems, corpusSize, nTopics, nTerms);
    
    Free(prob);
    Free(draw);
    Free(z_avg_j);
    Free(s_i);
    Free(B_j);
    Free(C_j);
    Free(newLambda);
    
    if(*verbose > 0) Rprintf("condProbSample_topic2: end\n");
}


// ----------------------------------------------------------------------------
//                              MCEM_EStep2
//
//  Do not draw multinomial sample; instead, use the probabilities directly
//
// ----------------------------------------------------------------------------
//  Notation: o_{ij} = y_{ij} - (alpha_i + beta_j + u_i' v_j + s_i' z_j)
//            gamma2:  gamma_i^2
//            o_gamma: o_{ij} * gamma_i
//
//  {o_gamma,alpha,beta,gamma,gamma2,u,v,s,z_avg}_mean are the Monte-Carlo means of o*gamma, alpha, ...
//  {alpha,beta,u,v,s}_sumvar are the sums of the Monte-Carlo variances over all o's, alpha's, ...
//  o_sum_adjvar: sum_{ij} ( E[o_{ij}^2] - (E[o_{ij}*gamma_i])^2 / E[gamma_i^2] )
//
//  eta_obj:    A vector of the objective values for different possible eta    values
//  lambda_obj: A vector of the objective values for different possible lambda values
//
//  if outputUserFactorVar == 1, then
//      {alpha,gamma,u,s}_outputVar will contain the Monte-Carlo variance (cov matrix) for each individual user
//  otherwise (outputUserFactorVar == 0), {alpha,gamma,u,s}_outputVar will not be changed
//
//  if outputItemFactorVar == 1, then
//      {beta,v}_outputVar will contain the Monte-Carlo variance (cov matrix) for each individual item
//  otherwise (outputItemFactorVar == 0), {beta,v}_outputVar will not be changed
//
//  nVar_{y,alpha,...} specifies the length of input var_{y,alpha,...}
//
//  SET nFactors = 0 to disable the u'v part
//  SET nTopics  = 0 to disable the s'z part
//  SET nCorpusWeights = 0 to give each term the same weight
//  SET nVar_{alpha,beta,...} = 0 to fix the factor values (i.e., prior variance = 0)
//  SET drawTopicSample = 0 to use the input z_avg_mean without drawing topic samples
//
// ----------------------------------------------------------------------------
void MCEM_EStep2(
    // INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
    double* alpha_mean/*nUsers x 1*/,    double* beta_mean/*nItems x 1*/,     double* gamma_mean/*nUsers x 1*/,
    double* u_mean/*nUsers x nFactors*/, double* v_mean/*nItems x nFactors*/, double* s_mean/*nUsers x nTopics*/,
    double* corpus_topic/*corpusSize x nTopics: the output consists of the topic for each term at the last draw*/,
    double* z_avg_mean/*nItems x nTopics: Use this input when drawTopicSample == 0*/,  
    // OUTPUT
    double* alpha_sumvar/*1x1*/,             double* alpha_outputVar/*nUsers x 1*/,
    double* beta_sumvar/*1x1*/,              double* beta_outputVar/*nItems x 1*/,
    double* gamma_sumvar/*1x1*/,             double* gamma_outputVar/*nUsers x 1*/,        double* gamma2_mean/*nUsers x 1*/,
    double* u_sumvar/*1x1*/,                 double* u_outputVar/*nUsers x nFactors x nFactors*/,
    double* v_sumvar/*1x1*/,                 double* v_outputVar/*nItems x nFactors x nFactors*/,
    double* s_sumvar/*1x1*/,                 double* s_outputVar/*nUsers x nTopics x nTopics*/,
    double* z_avg_outputVar/*nItems x nTopics x nTopics*/,
    double* eta_objval/*nEtaCanidates x 1*/, double* lambda_objval/*nLambdaCandidates x 1*/,
    double* o_gamma_mean/*nObs x 1*/,        double* o_sum_adjvar/*1x1*/,
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
    const double* eta,   const double* lambda,    const double* eta_candidates, const double* lambda_candidates,
    const double* nRatingExponent,
    const int* dim /*17 x 1*/,      const int* nDim /*must be 17*/,
    //  dim = {nObs, corpusSize, nUsers, nItems, nTerms, nFactors, nTopics, nCorpusWeights, nVar_y, 
    //         nVar_alpha, nVar_beta, nVar_gamma, nVar_u, nVar_v, nVar_s, nEtaCandidates, nLambdaCandidates}
    const int* outputUserFactorVar, const int* outputItemFactorVar,
    const int* outputUserTopicVar,  const int* outputItemTopicVar,
    const int* drawTopicSample/*MUST NOT ZERO if you want to do LDA*/,
    // OTHER
    const int* debug,  const int* verbose
){
    int *obsIndex_user, *oiStart_user, *oiNum_user, 
        *obsIndex_item, *oiStart_item, *oiNum_item,
        *cpsIndex,      *ciStart,      *ciNum,
        option=1, user_i, item_j;
    double *alpha, *beta, *gamma, *u, *v, *s,
           *o_gamma_sum, *o_sos, *alpha_sum, *alpha_sos, *beta_sum, *beta_sos, *gamma_sum, *gamma_sos, 
           *u_sum, *u_sos, *v_sum, *v_sos, *s_sum, *s_sos, *z_avg_sum,
           *rest, *cnt_item_topic, *cnt_topic_term, *cnt_topic, *z_avg;
    const int one = 1;
    int verbose_nextLevel = (*verbose) - 1;
    clock_t t_begin, t_begin_in;
    
    if(*verbose > 0) Rprintf("START MCEM_EStep.C\n");
    
    if(*nDim != 17) error("nDim should be 17: nDim=%d)",nDim);
    const int* nObs = dim+0;  const int* corpusSize = dim+1; const int* nUsers = dim+2; const int* nItems = dim+3; const int* nTerms = dim+4;
    const int* nFactors = dim+5; const int* nTopics = dim+6; const int* nCorpusWeights = dim+7;
    const int* nVar_y = dim+8; const int* nVar_alpha = dim+9; const int* nVar_beta = dim+10; const int* nVar_gamma = dim+11;
    const int* nVar_u = dim+12; const int* nVar_v = dim+13; const int* nVar_s = dim+14;
    const int* nEtaCandidates = dim+15; const int* nLambdaCandidates = dim+16;
    
    // Allocate space for sum and sum-of-squares (or sum of products of a pair of factors)
    alpha_sum      = (double*)Calloc(*nUsers, double);
    beta_sum       = (double*)Calloc(*nItems, double);
    u_sum          = (double*)Calloc((*nUsers)*(*nFactors), double);
    v_sum          = (double*)Calloc((*nItems)*(*nFactors), double);
    s_sum          = (double*)Calloc((*nUsers)*(*nTopics),  double);
    z_avg_sum      = (double*)Calloc((*nItems)*(*nTopics),  double);
    o_gamma_sum    = (double*)Calloc(*nObs, double);
    o_sos          = (double*)Calloc(*nObs, double);
    rest           = (double*)Calloc(*nObs, double);
    cnt_item_topic = (double*)Calloc((*nItems)*(*nTopics), double);
    cnt_topic_term = (double*)Calloc((*nTopics)*(*nTerms), double);
    cnt_topic      = (double*)Calloc(*nTopics, double);
    z_avg          = (double*)Calloc((*nItems)*(*nTopics), double);
    
    if(*nVar_gamma > 0){
        gamma_sum = (double*)Calloc(*nUsers, double);
        gamma_sos = (double*)Calloc(*nUsers, double);
    }

    if((*outputUserFactorVar) == 0){
        alpha_sos = (double*)Calloc(*nUsers, double);
        u_sos     = (double*)Calloc((*nUsers)*(*nFactors), double);
    }else if((*outputUserFactorVar) == 1){
        alpha_sos = NULL; u_sos = NULL;
        // use alpha_outputVar and u_outputVar
        for(int k=0; k<*nUsers; k++) alpha_outputVar[k] = 0;
        for(int k=0; k<(*nUsers)*(*nFactors)*(*nFactors); k++) u_outputVar[k] = 0;
        if(*nVar_gamma > 0){ for(int k=0; k<*nUsers; k++) gamma_outputVar[k] = 0;}
    }else error("outputUserFactorVar = %d should be only 0 or 1", *outputUserFactorVar);

    if((*outputItemFactorVar) == 0){
        beta_sos  = (double*)Calloc(*nItems, double);
        v_sos     = (double*)Calloc((*nItems)*(*nFactors), double);
    }else if((*outputItemFactorVar) == 1){
        beta_sos = NULL; v_sos = NULL;
        // use beta_outputVar and v_outputVar
        for(int k=0; k<*nItems; k++) beta_outputVar[k] = 0;
        for(int k=0; k<(*nItems)*(*nFactors)*(*nFactors); k++) v_outputVar[k] = 0;
    }else error("outputItemFactorVar = %d should be only 0 or 1", *outputItemFactorVar);

    if((*outputUserTopicVar) == 0){
        s_sos = (double*)Calloc((*nUsers)*(*nTopics), double);
    }else if((*outputUserTopicVar) == 1){
        s_sos = NULL; 
        // use s_outputVar
        for(int k=0; k<(*nUsers)*(*nTopics)*(*nTopics); k++) s_outputVar[k] = 0;
    }else error("outputUserTopicVar = %d should be only 0 or 1", *outputUserTopicVar);

    if((*outputItemTopicVar) == 1){
        for(int k=0; k<(*nItems)*(*nTopics)*(*nTopics); k++) z_avg_outputVar[k] = 0;
    }else if((*outputItemTopicVar) != 0) error("outputItemTopicVar = %d should be only 0 or 1", *outputItemTopicVar);
    
    for(int k=0; k<*nEtaCandidates;    k++) eta_objval[k]    = 0;
    for(int k=0; k<*nLambdaCandidates; k++) lambda_objval[k] = 0;
    
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
    if(*nTopics > 0) fillInTopicCounts2(cnt_item_topic, cnt_topic_term, cnt_topic, z_avg,
        corpus_topic, corpus_item, corpus_term, corpus_weight, nItems, corpusSize, nTopics, nTerms, nCorpusWeights, debug);
    if(*drawTopicSample == 0){
        for(int k=0; k<(*nItems)*(*nTopics); k++) z_avg[k] = z_avg_mean[k];
    }

    if(*nTopics > 0 && *drawTopicSample != 0)
        for(int k=0; k<(*nTopics)*(*nTerms); k++) phi[k] = 0;
    
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
            
            if(*drawTopicSample != 0){
                // Sample topics
                condProbSample_topic2(corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic, NULL,
                    &option, rest, s, var_y, eta, lambda, user, item, corpus_item, corpus_term, corpus_weight,
                    nUsers, nItems, nObs, corpusSize, nTopics, nTerms, nVar_y, &one, 
                    obsIndex_item, oiStart_item, oiNum_item, cpsIndex, ciStart, ciNum, nRatingExponent, debug, &verbose_nextLevel, &one
                );
                // Update z_avg
                compute_z_avg(z_avg, cnt_item_topic, nItems, nTopics);
            }
            
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
        // update gamma_sum and gamma_sos
        if(*nVar_gamma > 0){
            for(int k=0; k<*nUsers; k++) gamma_sum[k] += gamma[k];
            for(int k=0; k<*nUsers; k++) gamma_sos[k] += gamma[k]*gamma[k];
        }
        // update o_gamma_sum, o_sos
        for(int k=0; k<*nObs; k++){
            user_i = user[k]; item_j = item[k];
            if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems); }
            double uv = 0; for(int f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
            double sz = 0; for(int f=1; f<=*nTopics;  f++) sz += s[R_MAT(user_i,f,*nUsers)] * z_avg[R_MAT(item_j,f,*nItems)];
            double o = y[k] - alpha[R_VEC(user_i)] - beta[R_VEC(item_j)] - uv - sz;
            o_gamma_sum[k] += o * gamma[R_VEC(user_i)];
            o_sos[k] += o*o;
        }
        // update {eta,lambda}_objval
        if(*nTopics > 0){
            add_LDA_prior_objval_part2(eta_objval,    eta_candidates,    nEtaCandidates,    cnt_topic, cnt_topic_term, nTopics, nTerms);
            add_LDA_lambda_objval_part2(lambda_objval, lambda_candidates, nLambdaCandidates, cnt_item_topic, nItems, nTopics, nRatingExponent, oiNum_item);
        }
        
        // update sample variances of alpha and u
        if((*outputUserFactorVar) == 0){
            for(int k=0; k<*nUsers; k++)               alpha_sos[k] += alpha[k]*alpha[k];
            for(int k=0; k<(*nUsers)*(*nFactors); k++) u_sos[k]     += u[k]*u[k];
        }else if((*outputUserFactorVar) == 1){
            // update alpha
            for(int k=0; k<*nUsers; k++) alpha_outputVar[k] += alpha[k]*alpha[k];
            // update gamma
            if(*nVar_gamma > 0){ for(int k=0; k<*nUsers; k++) gamma_outputVar[k] += gamma[k]*gamma[k]; }
            // update u
            for(int k=0; k<(*nUsers); k++)
                // ONLY store the lower triangle
                for(int f1=0; f1<*nFactors; f1++) for(int f2=0; f2<=f1; f2++)
                    u_outputVar[C_3DA(k,f1,f2,*nUsers,*nFactors)] += u[C_MAT(k,f1,*nUsers)] * u[C_MAT(k,f2,*nUsers)];
        }else error("outputUserFactorVar = %d should be only 0 or 1", *outputUserFactorVar);

        // update sample variances of beta and v
        if((*outputItemFactorVar) == 0){
            for(int k=0; k<*nItems; k++)               beta_sos[k] += beta[k]*beta[k];
            for(int k=0; k<(*nItems)*(*nFactors); k++) v_sos[k]    += v[k]*v[k];
        }else if((*outputItemFactorVar) == 1){
            // update beta
            for(int k=0; k<*nItems; k++) beta_outputVar[k] += beta[k]*beta[k];
            // update v
            for(int k=0; k<(*nItems); k++)
                // ONLY store the lower triangle
                for(int f1=0; f1<*nFactors; f1++) for(int f2=0; f2<=f1; f2++)
                    v_outputVar[C_3DA(k,f1,f2,*nItems,*nFactors)] += v[C_MAT(k,f1,*nItems)] * v[C_MAT(k,f2,*nItems)];
        }else error("outputItemFactorVar = %d should be only 0 or 1", *outputItemFactorVar);
        
        // update sample variances of s
        if((*outputUserTopicVar) == 0){
            for(int k=0; k<(*nUsers)*(*nTopics); k++) s_sos[k] += s[k]*s[k];
        }else if((*outputUserTopicVar) == 1){
            for(int k=0; k<(*nUsers); k++)
                // ONLY store the lower triangle
                for(int f1=0; f1<*nTopics; f1++) for(int f2=0; f2<=f1; f2++)
                    s_outputVar[C_3DA(k,f1,f2,*nUsers,*nTopics)] += s[C_MAT(k,f1,*nUsers)] * s[C_MAT(k,f2,*nUsers)];
        }else error("outputUserTopicVar = %d should be only 0 or 1", *outputUserTopicVar);

        // update sample variances of z_avg
        if((*outputItemTopicVar) == 0){
            // do nothing
        }else if((*outputItemTopicVar) == 1){
            for(int k=0; k<(*nItems); k++)
                // ONLY store the lower triangle
                for(int f1=0; f1<*nTopics; f1++) for(int f2=0; f2<=f1; f2++)
                    z_avg_outputVar[C_3DA(k,f1,f2,*nItems,*nTopics)] += z_avg[C_MAT(k,f1,*nItems)] * z_avg[C_MAT(k,f2,*nItems)];
        }else error("outputItemTopicVar = %d should be only 0 or 1", *outputItemTopicVar);
        
        // update phi
        if(*nTopics > 0 && *drawTopicSample != 0)
            add_to_phi(phi, cnt_topic_term, cnt_topic, eta, nTopics, nTerms, debug);

        if(*verbose > 0){
            double secUsed = ((double)(clock() - t_begin_in)) / CLOCKS_PER_SEC;
            Rprintf(" + update: %.1f sec\n", secUsed);
        }
    }

    if((*outputUserFactorVar) == 0){
        computeMeanSumvar(alpha_mean, alpha_sumvar, alpha_sum, alpha_sos, *nUsers, *nSamples);
        computeMeanSumvar(u_mean, u_sumvar, u_sum, u_sos, (*nUsers)*(*nFactors), *nSamples);
        if(*nVar_gamma > 0) computeMeanSumvar(gamma_mean, gamma_sumvar, gamma_sum, gamma_sos, *nUsers, *nSamples);
    }else if((*outputUserFactorVar) == 1){
        computeMeanVar(alpha_mean, alpha_sumvar, alpha_outputVar, alpha_sum, *nUsers, 1, *nSamples);
        computeMeanVar(u_mean, u_sumvar, u_outputVar, u_sum, (*nUsers), (*nFactors), *nSamples);
        if(*nVar_gamma > 0) computeMeanVar(gamma_mean, gamma_sumvar, gamma_outputVar, gamma_sum, *nUsers, 1, *nSamples);
    }else error("outputUserFactorVar = %d should be only 0 or 1", *outputUserFactorVar);
    
    if(*nVar_gamma > 0){
        for(int i=0; i<*nUsers; i++) gamma2_mean[i] = gamma_sos[i] / (double)(*nSamples);
    }else{
        for(int i=0; i<*nUsers; i++) gamma2_mean[i] = gamma[i] * gamma[i];
    }

    if((*outputItemFactorVar) == 0){
        computeMeanSumvar(beta_mean, beta_sumvar,  beta_sum,  beta_sos,  *nItems, *nSamples);
        computeMeanSumvar(v_mean, v_sumvar, v_sum, v_sos, (*nItems)*(*nFactors), *nSamples);
    }else if((*outputItemFactorVar) == 1){
        computeMeanVar(beta_mean, beta_sumvar, beta_outputVar,  beta_sum, *nItems, 1, *nSamples);
        computeMeanVar(v_mean, v_sumvar, v_outputVar, v_sum, (*nItems), (*nFactors), *nSamples);
    }else error("outputItemFactorVar = %d should be only 0 or 1", *outputItemFactorVar);

    if((*outputUserTopicVar) == 0){
        computeMeanSumvar(s_mean, s_sumvar, s_sum, s_sos, (*nUsers)*(*nTopics), *nSamples);
    }else if((*outputUserTopicVar) == 1){
        computeMeanVar(s_mean, s_sumvar, s_outputVar, s_sum, (*nUsers), (*nTopics), *nSamples);
    }else error("outputUserTopicVar = %d should be only 0 or 1", *outputUserTopicVar);

    if((*outputItemTopicVar) == 0){
        for(int k=0; k<(*nItems)*(*nTopics); k++) z_avg_mean[k] = z_avg_sum[k] / (*nSamples);
    }else if((*outputItemTopicVar) == 1){
        double temp;
        computeMeanVar(z_avg_mean, &temp, z_avg_outputVar, z_avg_sum, (*nItems), (*nTopics), *nSamples);
    }else error("outputItemTopicVar = %d should be only 0 or 1", *outputItemTopicVar);

    (*o_sum_adjvar) = 0;
    for(int k=0; k<*nObs; k++){
        user_i = user[k]; item_j = item[k];
        if(*debug > 0){ CHK_R_INDEX(user_i, *nUsers); CHK_R_INDEX(item_j, *nItems); }
        o_gamma_mean[k] = o_gamma_sum[k] / (*nSamples);
        (*o_sum_adjvar) += o_sos[k] / (*nSamples) - (o_gamma_mean[k] * o_gamma_mean[k]) / gamma2_mean[R_VEC(user_i)];
    }

    // add first part for eta_objval, lambda_objval
    if(*nTopics > 0){
        for(int k=0; k<*nEtaCandidates; k++)
            eta_objval[k] = (*nTopics)*((*nTerms) * lgammafn(eta_candidates[k]) - lgammafn((*nTerms) * eta_candidates[k])) + eta_objval[k] / (*nSamples);
        for(int k=0; k<*nLambdaCandidates; k++)
            lambda_objval[k] =  
                compute_LDA_lambda_objval_part1(lambda_candidates[k], nItems, nTopics, nRatingExponent, oiNum_item) + 
                lambda_objval[k] / (*nSamples);
    }
    if(*nTopics > 0 && *drawTopicSample != 0)
        finalize_phi(phi, nTopics, nTerms, nSamples, debug);
    
    // Free the allocated space
    //   The R Free() would only free a point if it is NOT NULL.
    if((*outputUserFactorVar) == 0){ Free(alpha_sos); Free(u_sos); }
    if((*outputItemFactorVar) == 0){ Free(beta_sos); Free(v_sos); }
    if((*outputUserTopicVar) == 0) Free(s_sos);
    if((*nVar_gamma) > 0){ Free(gamma_sum); Free(gamma_sos); }
    Free(alpha_sum);      Free(beta_sum);       Free(u_sum);
    Free(v_sum);          Free(s_sum);          Free(z_avg_sum);
    Free(o_sos);          Free(o_gamma_sum);    Free(rest);
    Free(obsIndex_user);  Free(oiStart_user);   Free(oiNum_user);
    Free(obsIndex_item);  Free(oiStart_item);   Free(oiNum_item);
    Free(cpsIndex);       Free(ciStart);        Free(ciNum);
    Free(cnt_item_topic); Free(cnt_topic_term); Free(cnt_topic);   Free(z_avg);
    
    if(*verbose > 0) Rprintf("END   MCEM_EStep.C\n");
}


//-----------------------------------------------------------------------------
//   Utility Functions
//-----------------------------------------------------------------------------
void fillInTopicCounts2(
    // OUTPUT
    double *cnt_item_topic, double *cnt_topic_term, double *cnt_topic, double *z_avg,
    // INPUT
    const double *corpus_topic, const int *corpus_item, const int *corpus_term, const double *corpus_weight,
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms, const int *nCorpusWeights,
    const int *debug
){
    for(int k=0; k<(*nItems)*(*nTopics); k++) cnt_item_topic[k] = 0;
    for(int k=0; k<(*nTopics)*(*nTerms); k++) cnt_topic_term[k] = 0;
    for(int k=0; k<(*nTopics)          ; k++) cnt_topic[k]      = 0;
    
    for(int w=0; w<*corpusSize; w++){
        int thisItem  = corpus_item[w];
        int thisTerm  = corpus_term[w];
        double thisWeight = 1;
        if(corpus_weight != NULL && *nCorpusWeights != 0) thisWeight = corpus_weight[w];
        if(*debug > 0){ 
            CHK_R_INDEX(thisItem,  *nItems); 
            CHK_R_INDEX(thisTerm,  *nTerms); 
        }
        for(int k=1; k<=(*nTopics); k++){
            double prob = corpus_topic[C_MAT(w,k-1,*corpusSize)];
            cnt_item_topic[R_MAT(thisItem,k,*nItems)]  += thisWeight * prob;
            cnt_topic_term[R_MAT(k,thisTerm,*nTopics)] += thisWeight * prob;
            cnt_topic[R_VEC(k)] += thisWeight * prob;
        }
    }
    // print_matrix("cnt_item_topic: ", cnt_item_topic, *nItems, *nTopics);
    compute_z_avg(z_avg, cnt_item_topic, nItems, nTopics);
}


//-----------------------------------------------------------------------------
//   VALIDATE DATA STRUCTURES
//-----------------------------------------------------------------------------

void validate_corpus_counts2(
    const double *corpus_topic, const double *cnt_item_topic, const double *cnt_topic_term, const double *cnt_topic,
    const int *corpus_item, const int *corpus_term, const double *corpus_weight, 
    const int *nItems, const int *corpusSize, const int *nTopics, const int *nTerms
){
    int buffer_size = (*nTopics) * MAX((*nItems), (*nTerms));
    double *buffer = (double*)Calloc(buffer_size, double);
    double *cnt_topic_test = (double*)Calloc(*nTopics, double);
    
    // Rprintf("validate_corpus_counts2: buffer[%d], cnt_topic_test[%d] (%d, %d), nTerms=%d\n", buffer_size, *nTopics, buffer, cnt_topic_test, *nTerms);
    
    for(int w=0; w<*nTopics; w++) cnt_topic_test[w] = 0;
    
    for(int w=0; w<(*nItems)*(*nTopics); w++) buffer[w] = 0;
    for(int w=0; w<*corpusSize; w++){
        int thisItem  = corpus_item[w];
        double thisWeight = 1;
        CHK_R_INDEX(thisItem,  *nItems);
        if(corpus_weight != NULL) thisWeight = corpus_weight[w];
        for(int k=1; k<=*nTopics; k++){
            double prob = corpus_topic[C_MAT(w,k-1,*corpusSize)];
            buffer[R_MAT(thisItem, k, *nItems)] += thisWeight * prob;
            cnt_topic_test[R_VEC(k)] += thisWeight * prob;
        }
    }
    // print_matrix("cnt_item_topic: ", cnt_item_topic, *nItems, *nTopics);
    // print_matrix("        buffer: ", buffer,         *nItems, *nTopics);

    for(int w=0; w<(*nItems)*(*nTopics); w++)
        CHK_SAME_NUMBER("buffer[w] != cnt_item_topic[w]", buffer[w], cnt_item_topic[w]);
        
    for(int w=0; w<(*nTopics); w++)
        CHK_SAME_NUMBER("cnt_topic_test[w] != cnt_topic[w]", cnt_topic_test[w], cnt_topic[w]);

    for(int w=0; w<(*nTopics)*(*nTerms); w++) buffer[w] = 0;
    for(int w=0; w<*corpusSize; w++){
        int thisTerm  = corpus_term[w];
        double thisWeight = 1;
        CHK_R_INDEX(thisTerm,  *nTerms);
        if(corpus_weight != NULL) thisWeight = corpus_weight[w];
        for(int k=1; k<=*nTopics; k++){
            double prob = corpus_topic[C_MAT(w,k-1,*corpusSize)];
            buffer[R_MAT(k, thisTerm, *nTopics)] += thisWeight * prob;
        }
    }
    for(int w=0; w<(*nTopics)*(*nTerms); w++)
        CHK_SAME_NUMBER("buffer[w] != cnt_topic_term[w]", buffer[w], cnt_topic_term[w]);
    
    // Rprintf("validate_corpus_counts2: start to free (%d, %d)\n", buffer, cnt_topic_test);
    Free(buffer);
    Free(cnt_topic_test);   
    // Rprintf("validate_corpus_counts2: end\n");
}

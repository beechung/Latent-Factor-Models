/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


/**
 * Draw a sample for a main effect (alpha) from the posterior distribution
 *
 *	  obs[k] ~ N(alpha[thisEffIndex[k]], obsVar)
 *  alpha[i] ~ N(priorMean[i], priorVar[i])
 *
 *	Compute E[alpha|obs], Var[alpha|obs] and/or draw a sample from the posterior
 *
 * NOTE:
 *  thisEffIndex:  nObs x 1   obs: nObs x 1
 *  priorMean: nThisEff x 1   multiplier: nObs x 1
 *
 */
void gaussianPosterior_mainEffect(
    // OUTPUT
    double* sample        /* nLevelsThisEff x 1 */,
    double* posteriorMean /* nLevelsThisEff x 1 */,
    double* posteriorVar  /* nLevelsThisEff x 1 */,
    //INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* thisEffIndex /* nObs x 1 */,
    const double* obs /* nObs x 1 */,
    const double* priorMean /* nLevelsThisEff x 1 */,
    const double* multiplier /*NULL or nObs x 1*/,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const int* nObs, const int* nLevelsThisEff,
    const int* nObsVar   /* 1 or nObs */,
    const int* nPriorVar /* 1 or nLevelsThisEff */,
    // OTHER
    const int* debug
);

/**
 * Draw a sample for an effect (u) in a two-way interaction term from the posterior distribution
 *
 *	  obs[k] ~ N(u[thisEffIndex[k]]' v[otherEffIndex[k]], obsVar)
 *      u[i] ~ N(priorMean[i], priorVar[i])
 *
 *		u[i]: nFactors x 1 ( this effect)
 *		v[j]: nFactors x 1 (other effect)
 *
 *	Compute E[u|obs], Var[u|obs] and/or draw a sample from the posterior
 *
 * NOTE:
 *  Observation index (consider, say, user i, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
 *  obsIndex[ oiStart[i]+0 ], ..., obsIndex[ oiStart[i]+oiNum[i]-1 ] are the indices of user i's observations
 *  in y, x, user, item
 */
void gaussianPosterior_2WayInteraction(
    // OUTPUT
    double* sample,        /* nLevelsThisEff x nFactors */
    double* posteriorMean, /* nLevelsThisEff x nFactors */
    double* posteriorVar,  /* nLevelsThisEff x nFactors x nFactors */
    // INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* thisEffIndex  /* nObs x 1 */,
    const int* otherEffIndex /* nObs x 1 */,
    const double* obs        /* nObs x 1 */,
    const double* priorMean   /* nLevelsThisEff  x nFactors */,
    const double* otherEffect /* nLevelsOtherEff x nFactors */,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const int* nObs, const int* nLevelsThisEff, const int* nLevelsOtherEff, const int* nFactors,
    const int* nObsVar   /* 1 or nObs */,
    const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    // OTHER
    const int* debug
);

/**
 * Draw a sample for an effect (v) in a self interaction term from the posterior distribution
 *
 *	  obs[k] ~ N(v[fromIndex[k]]' v[toIndex[k]], obsVar)
 *      v[i] ~ N(priorMean[i], priorVar[i])
 *
 *		v[j]: nFactors x 1
 *
 *	Compute E[v|obs], Var[v|obs] and/or draw a sample from the posterior
 *
 * NOTE:
 *  Observation index (consider, say, user i, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
 *   from_obsIndex[ from_oiStart[i]+0 ], ..., from_obsIndex[ from_oiStart[i]+from_oiNum[i]-1 ] are the indices of  voter i's observations
 *     to_obsIndex[   to_oiStart[i]+0 ], ...,   to_obsIndex[   to_oiStart[i]+  to_oiNum[i]-1 ] are the indices of author i's observations
 *   in y, x_dyad, from, to
 */
void gaussianPosterior_SelfInteraction(
	// IN/OUT
    double* sample /* the current sample: nLevelsThisEff x nFactors */,
    // OUTPUT
    double* posteriorMean, /* nLevelsThisEff x nFactors */
    double* posteriorVar,  /* nLevelsThisEff x nFactors x nFactors */
    // INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* fromIndex /* nObs x 1 */,
    const int*   toIndex /* nObs x 1 */,
    const double* obs    /* nObs x 1 */,
    const double* priorMean /* nLevelsThisEff x nFactors */,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const int* nObs, const int* nLevelsThisEff, const int* nFactors,
    const int* nObsVar   /* 1 or nObs */,
    const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
    const int*   to_obsIndex, const int*   to_oiStart, const int*   to_oiNum,
    const int* from_obsIndex, const int* from_oiStart, const int* from_oiNum,
    // OTHER
    const int* debug
);


/**
 * Draw a sample for an effect (u) in a three-way interaction term from the posterior distribution
 *
 *	  obs[k] ~ N(sum(u[thisEffIndex[k],] * v[otherEffIndex[k],] * w[thirdEffIndex[k],]), obsVar)
 *      u[i] ~ N(priorMean[i], priorVar[i])
 *
 *		u[i]: nFactors x 1 ( this effect)
 *		v[j]: nFactors x 1 (other effect)
 *		w[j]: nFactors x 1 (third effect)
 *
 *	Compute E[u|obs], Var[u|obs] and/or draw a sample from the posterior
 *
 *     priorVar can be: 1x1, nFactors x 1 (one per factor), or nLevelsThisEff x nFactors x nFactors
 *
 * NOTE:
 *  Observation index (consider, say, user i, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
 *  obsIndex[ oiStart[i]+0 ], ..., obsIndex[ oiStart[i]+oiNum[i]-1 ] are the indices of user i's observations
 *  in y, x, user, item
 */
void gaussianPosterior_3WayInteraction(
    // OUTPUT
    double* sample,        /* nLevelsThisEff x nFactors */
    double* posteriorMean, /* nLevelsThisEff x nFactors */
    double* posteriorVar,  /* nLevelsThisEff x nFactors x nFactors */
    // INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* thisEffIndex  /* nObs x 1 */,
    const int* otherEffIndex /* nObs x 1 */,
    const int* thirdEffIndex /* nObs x 1 */,
    const double* obs        /* nObs x 1 */,
    const double* priorMean   /* nLevelsThisEff  x nFactors */,
    const double* otherEffect /* nLevelsOtherEff x nFactors */,
    const double* thirdEffect /* nLevelsThirdEff x nFactors */,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const int* nObs, const int* nLevelsThisEff, const int* nLevelsOtherEff, const int* nLevelsThirdEff,
    const int* nFactors, const int* nObsVar   /* 1 or nObs */,
    const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    // OTHER
    const int* debug
);

/**
 * Draw a sample for an effect (v) in a self interaction term from the posterior distribution
 *
 *	  obs[k] ~ N(sum(v[fromIndex[k],] * v[toIndex[k],] * w[thirdEffIndex[k],]), obsVar)
 *      v[i] ~ N(priorMean[i], priorVar[i])
 *
 *		v[j]: nFactors x 1
 *		w[j]: nFactors x 1 (third effect)
 *
 *	Compute E[v|obs], Var[v|obs] and/or draw a sample from the posterior
 *
 *     priorVar can be: 1x1, nFactors x 1 (one per factor), or nLevelsThisEff x nFactors x nFactors
 *
 * NOTE:
 *  Observation index (consider, say, user i, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
 *   from_obsIndex[ from_oiStart[i]+0 ], ..., from_obsIndex[ from_oiStart[i]+from_oiNum[i]-1 ] are the indices of  voter i's observations
 *     to_obsIndex[   to_oiStart[i]+0 ], ...,   to_obsIndex[   to_oiStart[i]+  to_oiNum[i]-1 ] are the indices of author i's observations
 *   in y, x_dyad, from, to
 */
void gaussianPosterior_SelfPlusOneInteraction(
	// IN/OUT
    double* sample /* the current sample: nLevelsThisEff x nFactors */,
    // OUTPUT
    double* posteriorMean, /* nLevelsThisEff x nFactors */
    double* posteriorVar,  /* nLevelsThisEff x nFactors x nFactors */
    // INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* fromIndex /* nObs x 1 */,
    const int*   toIndex /* nObs x 1 */,
    const int* thirdEffIndex /* nObs x 1 */,
    const double* obs    /* nObs x 1 */,
    const double* priorMean /* nLevelsThisEff x nFactors */,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const double* thirdEffect /* nLevelsThirdEff x nFactors */,
    const int* nObs, const int* nLevelsThisEff, const int* nLevelsThirdEff, const int* nFactors,
    const int* nObsVar   /* 1 or nObs */,
    const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
    const int* from_obsIndex, const int* from_oiStart, const int* from_oiNum,
    const int*   to_obsIndex, const int*   to_oiStart, const int*   to_oiNum,
    // OTHER
    const int* debug
);

/**
 * Draw a sample for a main effect (alpha) from the posterior distribution
 *
 *	        obs[m] ~ N(alpha[i=thisEffIndex[m],k=context[m]], obsVar[m])
 *      alpha[i,k] ~ N(contextOffset[i,k] + q[k]*alpha_global[i], contextPriorVar[i,k])
 * alpha_global[i] ~ N(0, globalPriorVar[i])
 *
 *	Compute E[alpha|obs], Var[alpha|obs] and/or draw a sample from the posterior
 *
 * NOTE:
 *   * contextPosteriorMean, contextPosteriorVar, globalPosteriorMean, globalPosteriorVar
 *     can NOT be NULL!!
 *   * contextSample can be set the same as contextPosteriorMean
 *     globalSample  can be set the same as globalPosteriorMean
 *     If so, the output will be the sample only.
 */
void gaussianPosterior_mainEffect_2Levels(
    // OUTPUT
    double* contextSample /* nLevelsThisEff x nContexts */,
    double* globalSample  /* nLevelsThisEff x 1 */,
    double* contextPosteriorMean /* nLevelsThisEff x nContexts: This is the posterior given the globalSample */,
    double* contextPosteriorVar  /* nLevelsThisEff x nContexts: This is the posterior given the globalSample */,
    double* globalPosteriorMean /* nLevelsThisEff x 1 */,
    double* globalPosteriorVar  /* nLevelsThisEff x 1 */,
    //INPUT
    const int* thisEffIndex /* nObs x 1 */,
    const int* context /* nObs x 1 */,
    const double* obs /* nObs x 1 */,
    const double* q /* nContext x 1 */,
    const double* contextOffset, /* nLevelsThisEff x nContexts */
    const double* obsVar   /* nObsVar x 1 */,
	const double* contextPriorVar /* see nContextPriorVar */,
	const double* globalPriorVar /*  see nGlobalPriorVar */,
	const int* numObs, const int* numLevelsThisEff, const int* numContexts,
	const int* numObsVar   /* 1 or nObs */,
	const int* numContextPriorVar /* nContexts or nLevelsThisEff*nContexts */,
	const int* numGlobalPriorVar  /* 1 or nLevelsThisEff */,
	const int* numContextOffset /* 0 or nLevelsThisEff or nLevelsThisEff*nContexts */,
    // OTHER
    const int* debug, const int* verbose
);

/**
 * Draw a sample for an interaction factor vector (u) from the posterior distribution
 *
 *	        obs[m] ~ N(mean = sum(u[i=thisEffIndex[m],  , k=thisContext[m]] *
 *	                              v[j=otherEffIndex[m], , k=otherContext[m]]),
 *	                   var = obsVar[m])
 *         u[i,,k] ~ N(offset[i,,k] + Q[,,k] %*% u_global[i,],  localPriorVar[k])
 *    u_global[i,] ~ N(0, globalPriorVar)
 *
 *	Compute E[u_global | obs] and Var[u_global | obs], and draw a sample from the posterior
 *	Then,   E[u | obs, u_global] and Var[u | obs, u_global] and draw a sample
 *
 * NOTE:
 *   * localPostMean, localPostVar, globalPostMean, globalPostVar can be NULL.
 *   * Observation index (consider, say, user i, starting from 0); but obs indices are R indices (starting from 1, NOT 0)
 *     obsIndex[ oiStart[i]+0 ], ..., obsIndex[ oiStart[i]+oiNum[i]-1 ] are the indices of user i's observations
 *     in y, x, user, item
 *   * If nThisContexts > 1, then localPriorVar[k] can be 0.  In this case,
 *     u[,,k] = u_global without drawing any random number.
 *   * If nThisContexts == 1, then nGlobalFactors and Q_size must be 0 and
 *     Q, globalSample, globalPostMean, globalPostVar, globalPriorVar and thisContext will not be touched.
 *   * If nOtherContexts == 1, then otherContext will not be touched.
 */
void gaussianPosterior_2WayInteraction_2Levels(
    // OUTPUT
    double* localSample    /* nLevelsThisEff x nLocalFactors x nThisContexts (u) */,
    double* globalSample   /* nLevelsThisEff x nGlobalFactors (u_global) */,
    double* localPostMean  /* NULL or nLevelsThisEff x nLocalFactors x nThisContexts: This is the posterior given the globalSample */,
    double* localPostVar   /* NULL or nLevelsThisEff x nLocalFactors x nLocalFactors x nThisContexts: This is the posterior given the globalSample */,
    double* globalPostMean /* NULL or nLevelsThisEff x nGlobalFactors */,
    double* globalPostVar  /* NULL or nLevelsThisEff x nGlobalFactors x nGlobalFactors */,
    //INPUT
    const int* option /*1: output sample, 2: output mean & var, 3: output sample & Mean & Var*/,
    const int* thisEffIndex  /* nObs x 1 */,
    const int* otherEffIndex /* nObs x 1 */,
    const int* thisContext   /* nObs x 1 */,
    const int* otherContext  /* nObs x 1 */,
    const double* obs        /* nObs x 1 */,
    const double* Q /* nLocalFactors x nGlobalFactors x nThisContexts */,
    const double* offset,  /* NULL or nLevelsThisEff x nLocalFactors x nThisContexts */
    const double* obsVar   /* nObsVar x 1 */,
	const double* localPriorVar  /* nThisContexts x 1 */,
	const double* globalPriorVar /* 1x1 */,
	const double* otherEff /* nLevelsOtherEff x nLocalFactors x nOtherContexts */,
	const int* nObs_, const int* nObsVar_ /* 1 or nObs */,
	const int* nLevelsThisEff_,  const int* nThisContexts_,
	const int* nLevelsOtherEff_, const int* nOtherContexts_,
	const int* nLocalFactors_,   const int* nGlobalFactors_,
	const int* nOffset_ /* 0 or nLevelsThisEff*nLocalFactors*nThisContexts */,
	const int* Q_size_ /* 1 or nLocalFactors*nGlobalFactors*nThisContexts */,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    // OTHER
    const int* debug_, const int* verbose_
);


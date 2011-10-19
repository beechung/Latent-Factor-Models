/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/



/**
 * Observation:
 * 		(itemIndex[n], categoryIndex[n], obs[n], var_obs[n])
 *                         i                 k
 *      This tuple specifies that we have observation obs[n] for item i in category k
 *      with observation variance var_obs[n].
 *      An item can be in multiple categories.
 *
 *		IMPORTANT NOTE: Indices start from 1 (NOT 0)
 *
 * Model:
 * 		obs[n] ~ N(mean=b[i,k],    var=var_obs[n])
 *      b[i,k] ~ N(mean=q[k]*a[i], var=var_b[k])
 *      a[i]   ~ N(mean=0,         var=var_a)
 *
 * Output:
 *		a_mean[i] =   E[a[i] | obs]
 *		 a_var[i] = Var[a[i] | obs]
 *	  b_mean[i,k] =   E[b[i,k] | obs]
 *	   b_var[i,k] = Var[b[i,k] | obs]
 *     b_cov[i,k] = Cov[b[i,k], a[i] | obs]
 *
 * Nodes with no observations will not contribute to the posterior.
 */
void hierarchical_smoothing_1D_2Levels(
	// Output
	double *a_mean, double *a_var, // nItems x 1
	double *b_mean, double *b_var, double *b_cov, // nItems x nCategories
	// Input
	const int *itemIndex, const int *categoryIndex, const double *obs, const double *var_obs, // nObs x 1
	const double *q, const double *var_b, // nCategories x 1
	const double *var_a, // 1x1
	const int *numObs, const int *numItems, const int *numCategories,
	// Control
	const int *verbose, const int *debug
);

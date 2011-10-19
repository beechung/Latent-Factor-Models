/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include "util.h"
#include <stdio.h>

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
){
	int nObs = (*numObs), nItems = (*numItems), nCategories = (*numCategories);

	// Compute sufficient statistics
	// 		b_mean[i,k] = sum_j { o_{ijk}/var_obs_{ijk} } = C_ik
	//       b_var[i,k] = sum_j {       1/var_obs_{ijk} } = F_ik
	for(int n=0; n < nItems*nCategories; n++){ b_mean[n] = 0; b_var[n] = 0; }
	for(int j=0; j<nObs; j++){
		double o = obs[j];
		int i = itemIndex[j] - 1;
		int k = categoryIndex[j] - 1;
		CHK_C_INDEX(i,nItems); CHK_C_INDEX(k,nCategories);
		// double var_obs_plus_var_b = var_obs[j] + var_b[k];
		b_mean[C_MAT(i,k,nItems)] += o / var_obs[j];
		b_var[C_MAT(i,k,nItems)]  += 1 / var_obs[j];
	}

	// Compute E[a[i] | obs] and Var[a[i] | obs]
	for(int i=0; i<nItems; i++){
		double sum_for_mean = 0; // sum_k E[a_i | obs_ik]/Var[a_i | obs_ik]
		double sum_for_var  = 0; // sum_k (1/Var[a_i | obs_ik] - 1/var_a)
		for(int k=0; k<nCategories; k++){
			double F = b_var[C_MAT(i,k,nItems)];
			if(F == 0) continue;
			double C = b_mean[C_MAT(i,k,nItems)];
			double A = (1/var_b[k]) + F;
			// E[a_i | obs_ik]/Var[a_i | obs_ik] = (C * q[k]) / (A * var_b[k])
			// (1/Var[a_i | obs_ik] - 1/var_a)   = (q[k]^2 / var_b[k]) * (1 - 1/(A * var_b[k]))
			sum_for_mean += (C * q[k]) / (A * var_b[k]);
			sum_for_var  += (q[k]*q[k] / var_b[k]) * (1 - 1/(A * var_b[k]));
		}
		if((*verbose) >= 100) printf("   i=%d:  sum.for.a.mean=%f, sum.for.a.var=%f\n", i+1, sum_for_mean, sum_for_var);

		// Compute E[a[i] | obs] and Var[a[i] | obs]
		a_var[i]  = 1/( (1/var_a[0]) + sum_for_var );
		a_mean[i] = a_var[i] * sum_for_mean;
	}

	// Compute E[b[i,k] | obs], Var[b[i,k] | obs], Cov[b[i,k], a[i] | obs]
	for(int i=0; i<nItems; i++){ for(int k=0; k<nCategories; k++){
		int ik = C_MAT(i,k,nItems);
		double A = 1 / ( (1/var_b[k]) + b_var[ik] );
		double D = (A * q[k]) / var_b[k];
		b_mean[ik] = D * a_mean[i] + A * b_mean[ik];
		b_var[ik]  = A + D*D*a_var[i];
		b_cov[ik]  = D * a_var[i];
	}}
}

void hierarchical_smoothing_1D_2Levels_incorrect(
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
){
	int nObs = (*numObs), nItems = (*numItems), nCategories = (*numCategories);

	// Compute sufficient statistics
	// 		  a_mean[i] = sum{ (q_k * o_{ijk})/(var_obs_{ijk} + var_b_{k}) }
	//         a_var[i] = sum{ (q_k * q_k)/(var_obs_{ijk} + var_b_{k}) }
	// 		b_mean[i,k] = sum{ o_{ijk}/var_obs_{ijk} }
	//       b_var[i,k] = sum{       1/var_obs_{ijk} }
	for(int i=0; i<nItems; i++){ a_mean[i] = 0; a_var[i] = 0; }
	for(int n=0; n < nItems*nCategories; n++){ b_mean[n] = 0; b_var[n] = 0; }
	for(int j=0; j<nObs; j++){
		double o = obs[j];
		int i = itemIndex[j] - 1;
		int k = categoryIndex[j] - 1;
		CHK_C_INDEX(i,nItems); CHK_C_INDEX(k,nCategories);
		double var_obs_plus_var_b = var_obs[j] + var_b[k];
		a_mean[i] += (q[k]*o) / var_obs_plus_var_b;
		a_var[i]  += (q[k]*q[k]) / var_obs_plus_var_b;
		b_mean[C_MAT(i,k,nItems)] += o / var_obs[j];
		b_var[C_MAT(i,k,nItems)]  += 1 / var_obs[j];
	}

	// Compute E[a[i] | obs] and Var[a[i] | obs]
	for(int i=0; i<nItems; i++){
		if((*verbose) >= 100) printf("   i=%d:  sum.for.a.mean=%f, sum.for.a.var=%f\n", i+1, a_mean[i], a_var[i]);
		a_var[i]  = 1/( (1/var_a[0]) + a_var[i] );
		a_mean[i] = a_var[i] * a_mean[i];
	}

	// Compute E[b[i,k] | obs], Var[b[i,k] | obs], Cov[b[i,k], a[i] | obs]
	for(int i=0; i<nItems; i++){ for(int k=0; k<nCategories; k++){
		int ik = C_MAT(i,k,nItems);
		double A = 1 / ( (1/var_b[k]) + b_var[ik] );
		double D = (A * q[k]) / var_b[k];
		b_mean[ik] = D * a_mean[i] + A * b_mean[ik];
		b_var[ik]  = A + D*D*a_var[i];
		b_cov[ik]  = D * a_var[i];
	}}
}

/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#ifndef DEFINITION_HPP_
#define DEFINITION_HPP_

#include <assert.h>
#include <stdio.h>

class ObsTable {
private:                  // Feature table:     Response table:
	const int *user_;     //   feature$user       response$user
	const int *app_;      //   feature$app        response$app
	const int *index_;    //   feature$index      response$item
	const double *value_; //   feature$x          response$y
	const double *w_;     //   feature$w          response$w
	const int nrow_;
	const int has_w_; // 0 or 1

	ReverseIndex userRevIndex;

public:

	/**
	 * Input from a R table: input = data.frame(user, app, index, x, w)
	 * where input$user, input$app, input$index are indices starting from 1 (not 0)
	 */
	ObsTable(
		const int *user, const int *app,  const int *index,
		const double *x, const double *w,
		const int nrow,  const int has_w,
		const int nUsers
	) : user_(user), app_(app), index_(index), value_(x), w_(w), nrow_(nrow),
	    has_w_(has_w), userRevIndex(user, nrow, nUsers)
	{
		assert(nrow_ >= 0);
		assert(has_w_ == 0 || has_w_ == 1);
	}

	// table$user[k]: user of the kth row; user starts from 0 and k starts from 0 (not 1)
	inline int user(int k) const { assert(k >= 0 && k<nrow_); return user_[k]-1; }

	// table$app[k]: app of the kth row; app starts from 0 and k starts from 0 (not 1)
	inline int app(int k) const { assert(k >= 0 && k<nrow_); return app_[k]-1; }

	// table$index[k]: index of the kth row; index starts from 0 and k starts from 0 (not 1)
	inline int index(int k) const { assert(k >= 0 && k<nrow_); return index_[k]-1; }

	inline double value(int k) const { assert(k >= 0 && k<nrow_); return value_[k];}

	inline double w(int k) const {
		if(!has_w_) return 1;
		assert(k >= 0 && k<nrow_);
		return w_[k];
	}

	inline bool has_w(void) const { return (has_w_ == 1); }

	inline int nrow() const { return nrow_; }

	// number of rows that user i has; user starts from 0 (not 1)
	inline int nRowsOfUser(int i) const {
		assert(i >= 0 && i < userRevIndex.size);
		return userRevIndex.num[i];
	}

	// get the row ID of the nth row of user i; row ID, user and n start from 0 (not 1)
	inline int getRowIdOfUser(int i, int n) const {
		assert(i >= 0 && i < userRevIndex.size);
		int id = userRevIndex.revIndex[i][n];
		assert(user_[id]-1 == i);
		return id;
	}

	/**
	 * Print the table for debugging
	 */
	void print(char* prefix, FILE* fp=stdout) const {
		for(int m=0; m<nrow(); m++){
			fprintf(fp,"%s %10d %10d %10d %11.7f %10.6f\n", prefix, user(m), app(m), index(m), value(m), w(m));
		}
	}
};

/**
 * Run one E-step and one M-step
 * INPUT:
 *    - Observations: feature and response
 *    - Parameters: A, B, b, alpha, beta, var_x, var_y, var_z, var_u
 *    - Options
 *    	  option & 0x01 == 0: In E-step, ignore z_{ik} if user i has no observation in app k
 *                         1:            include all z_{ik}
 *        option & 0x02 == 0: In M-step, ignore z_{ik} if user i has no observation in app k
 *                         1:            include all z_{ik}
 * OUTPUT:
 *    - Updated parameters: The contents of A, B, b, ... are replaced
 *      with the new values.
 *    - The posterior mean of the factors: They are stored in u and z.
 *    - xMean_var[ikm] = Var[B_{k,m} z_{ik}]
 *    - yMean_var[ijk] = Var[beta_{jk} z_{ik}]
 */
void run_EM_one_iteration(
	// Parameters: INPUT & OUTPUT
	ListOfMatrices& A,     ListOfMatrices& B,    ListOfMatrices& b,
	ListOfMatrices& alpha, ListOfMatrices& beta,
	double var_x[],  double var_y[],  double var_z[],  // vectors of length: nApps
	double& var_u, // scalar
	// Posterior mean of factors: OUTPUT
	Matrix_ColumnMajor& u,  ListOfMatrices& z,
	// Posterior variance for logistic regression
	double *xMean_var, const int xMean_var_length, // 0 or feature.nrow()
	double *yMean_var, const int yMean_var_length, // 0 or response.nrow()
	// Observation tables: INPUT
	ObsTable& feature,  ObsTable& response,
	// Size information: INPUT
	const int nApps, const int nUsers, const int nGlobalFactors,
	const int nFeatures[], const int nItems[], const int nLocalFactors[], // vectors of length: nApps
	// Others
	const double lambda_A, const double lambda_B, const double lambda_beta,
	const int option,  const int verbose,  const int debug
);


#endif /* DEFINITION_HPP_ */

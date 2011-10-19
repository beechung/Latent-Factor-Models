/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <assert.h>
#include "../../C/utils.hpp"
#include "definition.hpp"

/**
 * This is just a wrapper for the R-to-C interface.
 * The actual function is run_EM_one_iteration() (see definition.hpp)
 */
extern "C" void EM_one_iteration(
	// Parameters: INPUT and OUTPUT
	double *A_data, const int *A_dim,  double *B_data, const int *B_dim,  double *b_data, const int *b_dim,
	double *alpha_data, const int *alpha_dim,  double *beta_data, const int *beta_dim,
	double *var_x,  double *var_y,  double *var_z,  // vectors of length: nApps
	double *var_u, // scalar
	// Posterior mean of factors: OUTPUT
	double *u, // matrix of size nUsers x nGlobalFactors
	double *z_data, const int *z_dim,
	// Posterior variance for logistic regression: OUTPUT
	double *xMean_var, const int *xMean_var_length, // 0 or feature_nrow
	double *yMean_var, const int *yMean_var_length, // 0 or response_nrow
	// Feature table: INPUT
	const int *feature_user, const int *feature_app,  const int *feature_index,
	const double *feature_x, const double *feature_w,
	const int *feature_nrow, const int *feature_has_w,
	// Response table: INPUT
	const int *response_user, const int *response_app,  const int *response_item,
	const double *response_y, const double *response_w,
	const int *response_nrow, const int *response_has_w,
	// Ridge regression parameters: INPUT
	const double *ridge_lambda, const int *ridge_lambda_length, // ridge_lambda for A, B, beta
	// Size information: INPUT
	const int *nApps, const int *nUsers, const int *nGlobalFactors,    // scalars
	const int *nFeatures, const int *nItems, const int *nLocalFactors, // vectors of length: nApps
	// Others
	const int *option,  const int *verbose,  const int *debug // scalars
){
	if((*xMean_var_length) != 0 && (*xMean_var_length) != (*feature_nrow))
		STOP2("length(xMean_var) = %d, nrow(feature) = %d: They should be the same", *xMean_var_length, *feature_nrow);
	if((*yMean_var_length) != 0 && (*yMean_var_length) != (*response_nrow))
		STOP2("length(yMean_var) = %d, nrow(response) = %d: They should be the same", *yMean_var_length, *response_nrow);

	ListOfMatrices A(A_data, A_dim), B(B_data, B_dim), b(b_data, b_dim),
			       alpha(alpha_data, alpha_dim), beta(beta_data, beta_dim),
			       z(z_data, z_dim);
	Matrix_ColumnMajor u_matrix;
	u_matrix.wrap(u, *nUsers, *nGlobalFactors);
	ObsTable feature(feature_user,  feature_app,  feature_index, feature_x,  feature_w,  *feature_nrow,  *feature_has_w,  *nUsers),
			response(response_user, response_app, response_item, response_y, response_w, *response_nrow, *response_has_w, *nUsers);

	if(*verbose >= 100){
		// Print the input data for debugging
		printf("---------------------------------------------------------------\n"
			   "INPUT DATA to EM_one_iteration (C function)\n"
			   "(option:%d  verbose:%d  debug:%d\n"
			   "---------------------------------------------------------------\n",
			   *option, *verbose, *debug);
		printf("nApps=%d, nUsers=%d, nGlobalFeatures=%d\n", *nApps, *nUsers, *nGlobalFactors);
		printf("nFeatures:\n");       print_intVector("   ", nFeatures,     *nApps);
		printf("nItems:\n");          print_intVector("   ", nItems,        *nApps);
		printf("nLocalFactors:\n");   print_intVector("   ", nLocalFactors, *nApps);
		printf("Feature table:\n");   feature.print("   ");
		printf("Response table:\n");  response.print("   ");
		printf("param$A:\n");         A.print();
		printf("param$B:\n");         B.print();
		printf("param$b:\n");         b.print();
		printf("param$alpha:\n");     alpha.print();
		printf("param$beta:\n");      beta.print();
		printf("param$var_x:\n");	  print_vector("   ", var_x, *nApps);
		printf("param$var_y:\n");	  print_vector("   ", var_y, *nApps);
		printf("param$var_z:\n");	  print_vector("   ", var_z, *nApps);
		printf("param$var_u:\n");	  print_vector("   ", var_u, 1);
		printf("factor$u:\n");	      u_matrix.print("   ");
		printf("factor$z:\n");	      z.print();
		printf("ridge.lambda:\n");	  print_vector("   ", ridge_lambda, *ridge_lambda_length);
	}

	if(*ridge_lambda_length != 3) STOP1("length(redige.lambda) = %d (should be 3)", *ridge_lambda_length);

	run_EM_one_iteration(
		A, B, b, alpha, beta, var_x, var_y, var_z, *var_u,
		u_matrix, z,
		xMean_var, *xMean_var_length,
		yMean_var, *yMean_var_length,
		feature, response,
		*nApps, *nUsers, *nGlobalFactors, nFeatures, nItems, nLocalFactors,
		ridge_lambda[0], ridge_lambda[1], ridge_lambda[2],
		*option, *verbose, *debug
	);
}


//---------------------------------------------------------------------------
//  FUNCTIONS FOR DEBUGGING or AS EXAMPLES
//---------------------------------------------------------------------------
/**
 * Example of how to use ListOfMatrices and Matrix_ColumnMajor
 * Input:  list(X_1, X_2, X_3, X_4)
 * Output: X_1 %*% X_2 %*% X_3 %*% X_4
 */
extern "C" void product_listOfMatrices(
	// Output
	double* out,
	// Input
	const int* nrow_out, const int* ncol_out, // nrow(out), ncol(out)
	double* data, const int* dim, // packed list of matrices using src/R/utils.R/pack.list.of.matrices()
	const int* verbose
){
	ListOfMatrices list(data, dim);
	if(*verbose > 0){
		printf("Input List:\n");
		list.print();
	}
	Matrix_ColumnMajor A, temp;
	A = list.get(0); // copy the first matrix to A
	if(*verbose > 0){
		printf("Initial matrix:\n");
		A.print("  ");
	}
	for(int k=1; k<list.length(); k++){
		temp.product(A, list.get(k));
		A = temp; // copy the current result to A
		if(*verbose > 0){
			printf("After %d matrices:\n", k+1);
			A.print("  ");
		}
	}
	temp.wrap(out, *nrow_out, *ncol_out); // user temp to wrap the output matrix
	temp = A; // copy A to temp
}

/**
 * Example of how to use ListOfMatrices and Matrix_ColumnMajor
 * Input:  x = list(X_1, X_2, X_3)
 *         y = list(Y_1, Y_2, Y_3)
 * Output: list(...  (t(X_k)%*%X_k + diag(lambda_k))^-1 %*% t(X_k)%*%Y_k  ...);
 * NOTE: x, y and output are packed list of matrices using src/R/utils.R/pack.list.of.matrices()
 */
extern "C" void ridgeRegression_listOfMatrices(
	// Output
	double* output_data, const int* output_dim,
	// Input
	double* x_data,      const int* x_dim,
	double* y_data,      const int* y_dim,
	double* lambda_data, const int* lambda_dim,
	const int* verbose
){
	ListOfMatrices out(output_data, output_dim), x(x_data, x_dim), y(y_data, y_dim),
			       lambda(lambda_data, lambda_dim);
	Matrix_ColumnMajor A, B, temp;

	assert(out.length() == x.length() && out.length() == y.length() && out.length() == lambda.length());

	for(int k=0; k<out.length(); k++){
		B.transpose(x.get(k)); // B = t(X_k)
		temp.product(B, x.get(k));
		A = temp; // A = t(X_k)%*%X_k
		A.addDiagonal(lambda.get(k).getData(), lambda.get(k).length());
			// A = t(X_k)%*%X_k + diag(lambda)
		A.sym_invert(1); // A = (t(X_k)%*%X_k + diag(lambda))^-1
		temp.product(B, y.get(k));
		B = temp; // B = t(X_k)%*%Y_k
		temp.product(A, B);
		out.get(k) = temp;
	}
}

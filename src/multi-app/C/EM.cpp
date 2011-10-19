/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include "../../C/utils.hpp"
#include "definition.hpp"

inline bool isIdentity(const Matrix_ColumnMajor& A){
	if(A.length() == 1 && A(0,0) == 1) return true;
	return false;
}

// Get the data of user i from the observation table: obsTable,
//   which may be the response or feature table
//   OUTPUT:
//    * (value[k](m,0), index[k](m,0), w[k](m,0)) is the m-th observation
//         of user i in app k
//    * d[k], D[k] and C[k] are the parts of d_{ik}, D_{ik} and C_{ik}
//         for this obsTable.
void getUserDataForEStep_byApp(
	Matrix_ColumnMajor d[],     Matrix_ColumnMajor D[], Matrix_ColumnMajor C[],
	Matrix_ColumnMajor value[], Matrix_ColumnMajor w[], Matrix_ColumnMajor index[],
	const ListOfMatrices& inParam, const ListOfMatrices& inBias, const double inVar[],
	const ObsTable& obsTable, const int i,
	const int nApps, const int nLocalFactors[], const int nItems[]
){
	int *num = new int[nApps]; for(int k=0; k<nApps; k++) num[k] = 0;
	int n = obsTable.nRowsOfUser(i);
	for(int m=0; m<n; m++){
		int row_id = obsTable.getRowIdOfUser(i, m);
		int k = obsTable.app(row_id);
		assert(k >= 0 && k < nApps);
		num[k]++;
	}
	assert(inParam.length() == nApps);
	assert(inBias.length()  == nApps);
	for(int k=0; k<nApps; k++){
		d[k].resize(num[k], 1);
		D[k].resize(num[k], 1);
		value[k].resize(num[k], 1);
		w[k].resize(num[k], 1);
		index[k].resize(num[k], 1);

		if(isIdentity(inParam.get(k))){
			C[k].resize(num[k], nLocalFactors[k]);
		}else{
			assert(inParam.get(k).ncol() == nLocalFactors[k]);
			C[k].resize(num[k], nLocalFactors[k]);
		}
	}

	for(int k=0; k<nApps; k++) num[k] = 0;
	for(int m=0; m<n; m++){
		int row_id = obsTable.getRowIdOfUser(i, m);
		int k = obsTable.app(row_id);
		int j = obsTable.index(row_id);
		double y_ijk = obsTable.value(row_id);
		double w_ijk = obsTable.w(row_id);
		assert(j >= 0 && j < nItems[k]);

		value[k](num[k],0) = y_ijk;
		index[k](num[k],0) = j;
		w[k](num[k],0) = w_ijk;
		d[k](num[k],0) = y_ijk - inBias.get(k)(j,0);
		D[k](num[k],0) = w_ijk * inVar[k];
		if(isIdentity(inParam.get(k))){
			for(int a=0; a<nLocalFactors[k]; a++) C[k](num[k],a) = (a == j ? 1 : 0);
		}else{
			for(int a=0; a<nLocalFactors[k]; a++) C[k](num[k],a) = inParam.get(k)(j,a);
		}

		num[k]++;
	}

	delete[] num;
}

// Compute xMean_var[ikm] = Var[B_{k,m} z_{ik}]
//      or yMean_var[ijk] = Var[beta_{jk} z_{ik}]
// Variable names are based on yMean_var[ijk]
//
void compute_varOfScore(
	// output
	double *yMean_var, const int yMean_var_length, // 0 or response.nrow()
	const ListOfMatrices& beta, const Matrix_ColumnMajor *Var_z_i, // length: nApps
	const ObsTable& response, const int i,
	const int nApps, const int nLocalFactors[], const int nItems[]
){
	if(yMean_var_length == 0) return;
	if(yMean_var_length != response.nrow()) STOP("yMean_var_length != response.nrow()");
	int n = response.nRowsOfUser(i);

	for(int m=0; m<n; m++){
		int row_id = response.getRowIdOfUser(i, m);
		int k = response.app(row_id);
		int j = response.index(row_id);
		assert(j >= 0 && j < nItems[k]);
		double var = 0;
		if(isIdentity(beta.get(k))){
			var = Var_z_i[k].length() == 1 ? Var_z_i[k](0,0) : Var_z_i[k](j,j);
		}else{
			for(int a=0; a<nLocalFactors[k]; a++) for(int b=0; b<nLocalFactors[k]; b++){
				var += beta.get(k)(j,a) * Var_z_i[k](a,b) * beta.get(k)(j,b);
			}
		}
		yMean_var[row_id] = var;
	}
}

/**
 * Update the sufficient stats for either features or response
 * called in run_E_step()
 * If we work with features,
 * 		XX =  XX_x[k] defined in the comment before run_E_step()
 * 		XY =  XY_z[k] defined in the comment before run_E_step()
 * 		Y2 =  Y2_z[k] defined in the comment before run_E_step()
 * 	   num = num_z[k] defined in the comment before run_E_step()
 *    (x_ik[a], w_ik[a], index[a]) = the value, weight and feature index of the a-th obs
 */
void updateSufficientStats(
	Matrix_ColumnMajor XX[], Matrix_ColumnMajor XY[], double& Y2, int& num,
	const Matrix_ColumnMajor& E_z_ik, const Matrix_ColumnMajor& Var_z_ik,
	const Matrix_ColumnMajor& x_ik, const Matrix_ColumnMajor& w_ik, const Matrix_ColumnMajor& index,
	const Matrix_ColumnMajor& b_k, const int verbose
){
	int N = x_ik.nrow();
	assert((XX != NULL && XY != NULL) || (XX == NULL && XY == NULL));
	assert(x_ik.ncol() == 1 && w_ik.ncol() == 1 && index.ncol() == 1);
	assert(x_ik.nrow() == N && w_ik.nrow() == N && index.nrow() == N);

	if(N == 0) return;

	// Let z1_{ik} = c(1, z_{ik})
	//       E_z1 = E[z1_{ik}]
	//     Var_z1 = Var[z1_{ik}]
	Matrix_ColumnMajor E_z1, Var_z1, temp, temp2;
	E_z1.resize(E_z_ik.nrow()+1,1);
	E_z1(0,0) = 1;
	for(int m=0; m<E_z_ik.nrow(); m++) E_z1(m+1,0) = E_z_ik(m,0);
	Var_z1.resize(E_z_ik.nrow()+1, E_z_ik.nrow()+1);
	Var_z1.setToZero();
	if(Var_z_ik.length() == 1){
		Var_z1.addDiagonal(Var_z_ik);
		Var_z1(0,0) = 0;
	}else{
		for(int i=1; i<Var_z1.nrow(); i++){
			for(int j=1; j<Var_z1.ncol(); j++) Var_z1(i,j) = Var_z_ik(i-1,j-1);
		}
	}

	for(int a=0; a<N; a++){
		int m = (int)index(a,0);
		double x_ikm = x_ik(a,0);
		double w_ikm = w_ik(a,0);
		if(XX != NULL && XY != NULL){
			// XX_x[k][m] = sum_{i in I_{km}} { (Var[z1_{ik}] + E[z1_{ik}]E[z1_{ik}]') / w_{x,ikm} }
			temp.transpose(E_z1);
			temp2.product(E_z1, temp);
			temp2.add(Var_z1);
			temp2.scale(1/w_ikm);
			XX[m].add(temp2);
			// XY_x[k][m] = sum_{i in I_{km}} { E[z1_{ik}] x_{ik,m} / w_{x,ikm} }
			temp = E_z1;
			temp.scale(x_ikm/w_ikm);
			XY[m].add(temp);
			// Y2_x[k] = sum_{m} sum_{i in I_{km}} { x_{ik,m}^2 / w_{x,ikm} }
			Y2 += (x_ikm * x_ikm / w_ikm);
		}else{
			// Y2_x[k] = sum_{m} sum_{i in I_{km}} { ((x_{ik,m} - E[z_{ik,m}] - b_{k,m})^2 + Var[z_{ik,m}]) / w_{x,ikm} }
			double err = x_ikm - E_z_ik(m,0) - b_k(m,0);
			Y2 += ((err*err + Var_z_ik(m,m)) / w_ikm);
		}
		// num_x[k] = sum_{m} sum_{i in I_{km}} 1
		num++;
	}
}

/**
 * To save space, the EM algorithm is based on the following sufficient statistics.
 * 	   The E-step computes the sufficient stats
 *     The M-step uses the sufficient stats
 * Sufficient stats for var_u
 * 	   sos_u = sum_{i in I} sum_{m} { E[u_{i,m}]^2 + Var[u_{i,m}] }
 * 	   num_u = sum_{i in I} sum_{m} 1
 * Sufficient stats for A, var_z
 *   If A_k != identity:
 *     XX_z[k]     = sum_{i in I_k} { (Var[u_i] + E[u_i]E[u_i]') }
 *     XY_z[k](,m) = sum_{i in I_k} { E[u_i]E[z_{ik,m}] + Cov[z_{ik,m},u_i] }
 *     Y2_z[k]     = sum_{i in I_k} sum_{m} { E[z_{ik,m}]^2 + Var[z_{ik,m}] }
 *    num_z[k]     = sum_{i in I_k} sum_{m} 1
 *     index: m-th local factor in app k
 *   If A_k == identity:
 *     XX_z[k] = XY_z[k] = zero-length matrix
 *     Y2_z[k] = sum_{i in I_k} sum_{m} { (E[z_{ik,m}] - E[u_{i,m}])^2 + Var[z_{ik,m}] + Var[u_{i,m}] - 2*Cov[z_{ik,m},u_{i,m}] }
 * Sufficient stats for B, b, var_x
 *   If B_k != identity:
 *     Let z1_{ik} = c(1, z_{ik})
 *     XX_x[k][m]  = sum_{i in I_{km}} { (Var[z1_{ik}] + E[z1_{ik}]E[z1_{ik}]') / w_{x,ikm} }
 *     XY_x[k][m]  = sum_{i in I_{km}} { E[z1_{ik}] x_{ik,m} / w_{x,ikm} }
 *     Y2_x[k]     = sum_{m} sum_{i in I_{km}} { x_{ik,m}^2 / w_{x,ikm} }
 *    num_x[k]     = sum_{m} sum_{i in I_{km}} 1
 *     index: m-th feature in app k
 *   If B_k == identity:
 *     XX_x[k] = XY_x[k] = NULL
 *     Y2_x[k] = sum_{m} sum_{i in I_{km}} { ((x_{ik,m} - E[z_{ik,m}] - b_{k,m})^2 + Var[z_{ik,m}]) / w_{x,ikm} }
 * Sufficient stats for beta, alpha, var_y
 *     Let z1_{ik} = c(1, z_{ik})
 *     XX_y[k][j]  = sum_{i in I_{jk}} { (Var[z1_{ik}] + E[z1_{ik}]E[z1_{ik}]') / w_{y,ijk} }
 *     XY_y[k][j]  = sum_{i in I_{jk}} { E[z1_{ik}] y_{ijk} / w_{y,ijk} }
 *     Y2_y[k]     = sum_{j} sum_{i in I_{jk}} { y_{ijk}^2 / w_{y,ijk} }
 *    num_y[k]     = sum_{j} sum_{i in I_{jk}} 1
 *     index: item j in app k
 */
void run_E_step(
	// Posterior mean of factors: OUTPUT
	Matrix_ColumnMajor& u,  ListOfMatrices& z,
	// Sufficient stats: OUTPUT
	double& sos_u, double& num_u,
	Matrix_ColumnMajor  *XX_z, Matrix_ColumnMajor  *XY_z, double Y2_z[], int num_z[],
	Matrix_ColumnMajor **XX_x, Matrix_ColumnMajor **XY_x, double Y2_x[], int num_x[],
	Matrix_ColumnMajor **XX_y, Matrix_ColumnMajor **XY_y, double Y2_y[], int num_y[],
	// Posterior variance for logistic regression
	double *xMean_var, const int xMean_var_length, // 0 or feature.nrow()
	double *yMean_var, const int yMean_var_length, // 0 or response.nrow()
	// Parameters: INPUT
	const ListOfMatrices& A,     const ListOfMatrices& B,   const ListOfMatrices& b,
	const ListOfMatrices& alpha, const ListOfMatrices& beta,
	const double var_x[], const double var_y[], const double var_z[],  // vectors of length: nApps
	const double& var_u, // scalar
	// Observation tables: INPUT
	const ObsTable& feature, const ObsTable& response,
	// Size information: INPUT
	const int nApps, const int nUsers, const int nGlobalFactors,
	const int nFeatures[], const int nItems[], const int nLocalFactors[], // vectors of length: nApps
	// Others
	const int option,  const int verbose,  const int debug
){
	// --------------------------------------------------------------
	// Initialization
	// --------------------------------------------------------------
	if(verbose >= 5) printf("   Initialize variables\n");
	sos_u = 0;  num_u = 0;
	for(int k=0; k<nApps; k++){
		XX_z[k].setToZero();  XY_z[k].setToZero();  Y2_z[k] = 0;  num_z[k] = 0;
		for(int m=0; m<nFeatures[k]; m++){
			if(XX_x[k] != NULL) XX_x[k][m].setToZero();
			if(XX_x[k] != NULL) XY_x[k][m].setToZero();
		}
		Y2_x[k] = 0; num_x[k] = 0;
		for(int m=0; m<nItems[k]; m++){ XX_y[k][m].setToZero(); XY_y[k][m].setToZero(); }
		Y2_y[k] = 0; num_y[k] = 0;
	}
	Matrix_ColumnMajor *d = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *C = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *D = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *H = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *R = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *d_y = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *C_y = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *D_y = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *d_x = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *C_x = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *D_x = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *x_i = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *y_i = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *w_x = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *w_y = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *index_x = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *index_y = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *Sigma     = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *Sigma_inv = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *z_ik_ik   = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *Gamma_ik_ik = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *Gamma_i_ik  = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *Gamma_i_ik_inv = new Matrix_ColumnMajor[nApps];

	Matrix_ColumnMajor F, D_inv, CD, CDC, CDd, SigmaCDC, SigmaCDCF, GammaHtGamma, GammaHGamma, temp, temp2,
	                   E_u_i, Var_u_i, E_z_ik, Var_z_ik, Cov_z_u;

	Matrix_ColumnMajor *Var_z_i = NULL;
	if(xMean_var_length > 0 || yMean_var_length > 0) Var_z_i = new Matrix_ColumnMajor[nApps];

	// Prepare Sigma, H, R
	if(verbose >= 5) printf("   Prepare Sigma, H, R\n");
	for(int k=0; k<nApps; k++){
		if(isIdentity(A.get(k))){
			Sigma[k].resize(1,1);
			Sigma[k](0,0) = var_u + var_z[k];
			Sigma_inv[k].resize(1,1);
			Sigma_inv[k](0,0) = 1.0/(var_u + var_z[k]);
			H[k].resize(1,1);
			H[k](0,0) = var_u * Sigma_inv[k](0,0);
			R[k].resize(1,1);
			R[k](0,0) = var_u - var_u * var_u * Sigma_inv[k](0,0);
		}else{
			temp.transpose(A.get(k));       // temp = t(A_k)
			Sigma[k].product(A.get(k),temp); // Sigma[k] = A_k %*% t(A_k)
			Sigma[k].scale(var_u);           // Sigma[k] = var_u * A_k %*% t(A_k)
			Sigma[k].addDiagonal(&(var_z[k]), 1);
			// Now, Sigma[k] = var_u * A_k %*% t(A_k) + var_z[k] * I
			Sigma_inv[k] = Sigma[k];
			Sigma_inv[k].sym_invert(debug);
			H[k].product(temp, Sigma_inv[k]);
			H[k].scale(var_u);
			// Now, H[k] = var_u * t(A_k) %*% Sigma[k]^-1
			R[k].product(H[k], A.get(k));
			R[k].scale(var_u);
			R[k].negate();
			R[k].addDiagonal(&var_u, 1);
			// Now, R[k] = var_u * I - var_u^2 * t(A_k) %*% Sigma[k]^-1 %*% A_k
		}
	}

	// --------------------------------------------------------------
	// Loop through each user i
	// --------------------------------------------------------------
	if(verbose >= 5) printf("   Loop through all users\n");
	for(int i=0; i<nUsers; i++){
		if(verbose >= 10) printf("   USER %d\n", i);
		// Prepare d, C, D
		if(verbose >= 10) printf("      Prepare d, C, D\n");
		getUserDataForEStep_byApp(d_x, D_x, C_x, x_i, w_x, index_x, B,    b,     var_x, feature,  i, nApps, nLocalFactors, nFeatures);
		getUserDataForEStep_byApp(d_y, D_y, C_y, y_i, w_y, index_y, beta, alpha, var_y, response, i, nApps, nLocalFactors, nItems);
		for(int k=0; k<nApps; k++){
			d[k].rbind(d_x[k], d_y[k]);
			D[k].rbind(D_x[k], D_y[k]);
			C[k].rbind(C_x[k], C_y[k]);
		}
		// Go up-tree
		if(verbose >= 10) printf("      Go up-tree\n");
		bool user_i_has_no_data = true;
		for(int k=0; k<nApps; k++){
			if(d[k].length() == 0){
				z_ik_ik[k].resize(nLocalFactors[k],1);
				z_ik_ik[k].setToZero();
				Gamma_ik_ik[k] = Sigma[k];
			}else{
				user_i_has_no_data = false;
				D_inv = D[k];
				D_inv.elementwiseInverse();  // D_inv = (D_{ik})^-1
				temp.transpose(C[k]);        // temp = t(C_{ik})
				CD.product_2ndDiagonal(temp, D_inv);
				CDC.product(CD, C[k]); // CDC = t(C_{ik}) %*% (D_{ik})^-1 %*% C_{ik}
				CDd.product(CD, d[k]); // CDd = t(C_{ik}) %*% (D_{ik})^-1 %*% d_{ik}
				F = CDC;
				F.add(Sigma_inv[k], true);
				F.sym_invert(debug); // F = (F_{ik})^-1
				SigmaCDC.product(Sigma[k], CDC, true);
				SigmaCDCF.product(SigmaCDC, F);
				// compute z_ik_ik[k] = z_{ik|ik}
				z_ik_ik[k].product(Sigma[k], CDd, true);
				temp.product(SigmaCDCF, CDd);
				z_ik_ik[k].subtract(temp);
				// compute Gamma_ik_ik[k] = Gamma_{ik|ik}
				Gamma_ik_ik[k].product(SigmaCDC, Sigma[k], true);
				temp.product(SigmaCDCF,CDC);
				temp2.product(temp, Sigma[k], true);
				Gamma_ik_ik[k].subtract(temp2);
				Gamma_ik_ik[k].negate();
				Gamma_ik_ik[k].add(Sigma[k], true);
			}
			// compute Gamma_i_ik_inv[k] = Gamma_{i|ik}^-1
			temp.product(H[k], Gamma_ik_ik[k], true);
			temp2.transpose(H[k]);
			Gamma_i_ik[k].product(temp, temp2, true);
			Gamma_i_ik[k].add(R[k], true);
			Gamma_i_ik_inv[k] = Gamma_i_ik[k];
			Gamma_i_ik_inv[k].sym_invert(debug);  //TODO: is this symmetric??
		}
		// compute E[u_i] and Var[u_i]
		if(verbose >= 10) printf("      Compute E[u_i] and Var[u_i]\n");
		double one_over_var_u = 1.0/var_u;
		Var_u_i.resize(nGlobalFactors, nGlobalFactors);
		Var_u_i.setToZero();
		Var_u_i.addDiagonal(&one_over_var_u,1);
		E_u_i.resize(nGlobalFactors,1);
		E_u_i.setToZero();
		for(int k=0; k<nApps; k++){
			if((option & 0x01) == 0 && d[k].length() == 0) continue;
			temp = Gamma_i_ik_inv[k];
			double neg_one_over_var_u = -1.0/var_u;
			temp.addDiagonal(&neg_one_over_var_u,1);
			Var_u_i.add(temp);
			temp.product(Gamma_i_ik_inv[k], H[k], true);
			temp2.product(temp, z_ik_ik[k]);
			E_u_i.add(temp2);
		}
		Var_u_i.sym_invert(debug); // TODO: is this symmetric??
		temp = E_u_i;
		E_u_i.product(Var_u_i, temp);
		// Now, Var_u_i = Var[u_i]
		//        E_u_i =   E[u_i]
		for(int m=0; m<nGlobalFactors; m++) u(i,m) = E_u_i(m,0);

		// update sos_u and num_u
		if(verbose >= 10) printf("      Update sos_u and num_u\n");
		if(!user_i_has_no_data){
			for(int m=0; m<nGlobalFactors; m++){
				 // sos_u = sum_{i in I} sum_{m} { E[u_{i,m}]^2 + Var[u_{i,m}] }
				 // num_u = sum_{i in I} sum_{m} 1
				sos_u += E_u_i(m,0) * E_u_i(m,0) + Var_u_i(m,m);
				num_u ++;
			}
		}

		// Go down-tree
		if(verbose >= 10) printf("      Go down-tree\n");
		for(int k=0; k<nApps; k++){
			if(verbose >= 11) printf("         APPLICATION k = %d\n", k);
			temp.transpose(H[k]); // temp = t(H_{ik})
			temp2.product(Gamma_ik_ik[k], temp, true);
			GammaHtGamma.product(temp2, Gamma_i_ik_inv[k], true);
			GammaHGamma.transpose(GammaHtGamma);
			// compute E_z_ik = E[z_{ik}]
			if(verbose >= 11) printf("         compute E[z_{ik}]\n");
			temp.product(H[k], z_ik_ik[k], true);
			temp.negate();
			temp.add(E_u_i);
			E_z_ik.product(GammaHtGamma, temp, true);
			E_z_ik.add(z_ik_ik[k]);
			// compute Var_z_ik = Var[z_{ik}]
			if(verbose >= 11) printf("         compute Var[z_{ik}]\n");
			temp = Var_u_i;
			temp.subtract(Gamma_i_ik[k], true);
			temp2.product(GammaHtGamma,temp, true);
			Var_z_ik.product(temp2, GammaHGamma, true);
			Var_z_ik.add(Gamma_ik_ik[k], true);
			// compute Cov_z_u = Cov[z_{ik}, u_i]
			if(verbose >= 11) printf("         compute Cov[z_{ik}, u_i]\n");
			Cov_z_u.product(GammaHtGamma, Var_u_i, true);

			// update z
			for(int m=0; m<nLocalFactors[k]; m++) z.get(k)(i,m) = E_z_ik(m,0);

			if(Var_z_i != NULL) Var_z_i[k] = Var_z_ik;

			// update XX_z, XY_z, Y2_z, num_z
			if(option & 0x02 == 1 || d[k].length() > 0){
				if(XX_z[k].length() == 0 && XY_z[k].length() == 0){
					//  Y2_z[k] = sum_{i in I_k} sum_{m} { (E[z_{ik,m}] - E[u_{i,m}])^2 + Var[z_{ik,m}] + Var[u_{i,m}] - 2*Cov[z_{ik,m},u_{i,m}] }
					// num_z[k] = sum_{i in I_k} sum_{m} 1
					if(verbose >= 11) printf("         update Y2_z[%d] and num_z[%d]\n",k,k);
					for(int m=0; m<nLocalFactors[k]; m++){
						double E_z_ikm = E_z_ik(m,0), E_u_im = E_u_i(m,0);
						double V_z_ikm = Var_z_ik.length() == 1 ? Var_z_ik(0,0) : Var_z_ik(m,m);
						double V_u_im  = Var_u_i(m,m);
						double Cov     = Cov_z_u(m,m);
						Y2_z[k] += (E_z_ikm-E_u_im)*(E_z_ikm-E_u_im) + V_z_ikm + V_u_im - 2 * Cov;
						num_z[k]++;
					}
				}else{
					// XX_z[k] = sum_{i in I_k} { (Var[u_i] + E[u_i]E[u_i]') }
					if(verbose >= 11) printf("         update XX_z[%d]\n",k);
					temp.transpose(E_u_i);
					temp2.product(E_u_i, temp);
					temp2.add(Var_u_i, true);
					XX_z[k].add(temp2);
					// XY_z[k](,m) = sum_{i in I_k} { E[u_i]E[z_{ik,m}] + Cov[z_{ik,m},u_i] }
					// XY_z[k]     = sum_{i in I_k} { E[u_i]%*%t(E[z_{ik}]) + t(Cov[z_{ik},u_i]) }
					if(verbose >= 11) printf("         update XY_z[%d]\n",k);
					temp.transpose(E_z_ik);
					temp2.product(E_u_i, temp);
					temp.transpose(Cov_z_u);
					temp2.add(temp);
					XY_z[k].add(temp2);
					//  Y2_z[k] = sum_{i in I_k} sum_{m} { E[z_{ik,m}]^2 + Var[z_{ik,m}] }
					// num_z[k] = sum_{i in I_k} sum_{m} 1
					if(verbose >= 11) printf("         update Y2_z[%d] and num_z[%d]\n",k,k);
					for(int m=0; m<nLocalFactors[k]; m++){
						double E_z_ikm = E_z_ik(m,0);
						double V_z_ikm = Var_z_ik.length() == 1 ? Var_z_ik(0,0) : Var_z_ik(m,m);
						Y2_z[k] += E_z_ikm*E_z_ikm + V_z_ikm;
						num_z[k]++;
					}
				}
			}
			// update XX_x, XY_x, Y2_x, num_x
			if(verbose >= 11) printf("         update XX_x[%d], XY_x[%d], Y2_x[%d], num_x[%d]\n",k,k,k,k);
			updateSufficientStats(
				XX_x[k], XY_x[k], Y2_x[k], num_x[k],
				E_z_ik, Var_z_ik, x_i[k], w_x[k], index_x[k], b.get(k), verbose
			);
			// update XX_y, XY_y, Y2_y, num_y
			if(verbose >= 11) printf("         update XX_y[%d], XY_y[%d], Y2_y[%d], num_y[%d]\n",k,k,k,k);
			updateSufficientStats(
				XX_y[k], XY_y[k], Y2_y[k], num_y[k],
				E_z_ik, Var_z_ik, y_i[k], w_y[k], index_y[k], alpha.get(k), verbose
			);
		}
		compute_varOfScore(
			xMean_var, xMean_var_length, B, Var_z_i, feature, i, nApps, nLocalFactors, nFeatures
		);
		compute_varOfScore(
			yMean_var, yMean_var_length, beta, Var_z_i, response, i, nApps, nLocalFactors, nItems
		);
	}
	if(verbose >= 5) printf("   Finish processing all users\n");

	delete[] d;  delete[] C;  delete[] D;  delete[] H;  delete[] R;
	delete[] d_y;  delete[] C_y;  delete[] D_y;  delete[] index_y;  delete[] y_i;  delete[] w_x;
	delete[] d_x;  delete[] C_x;  delete[] D_x;  delete[] index_x;  delete[] x_i;  delete[] w_y;
	delete[] Sigma;    delete[] Sigma_inv;    delete[] Gamma_i_ik;
	delete[] z_ik_ik;  delete[] Gamma_ik_ik;  delete[] Gamma_i_ik_inv;
	if(Var_z_i != NULL) delete[] Var_z_i;

	if(verbose >= 5) printf("   Finish deallocating variables\n");
}

void estimateParamVar_oneApp(
	// OUTPUT
	Matrix_ColumnMajor& B, Matrix_ColumnMajor& b, double& var_x,
	// INPUT
	Matrix_ColumnMajor *XX_x, Matrix_ColumnMajor *XY_x, const double Y2_x, const int num_x,
	const double lambda, const int nFeatures, const int nLocalFactors,
	const int verbose, const int debug
){
	if(num_x == 0) return;
	Matrix_ColumnMajor XX_inv, temp, temp2, eta;
	if(isIdentity(B)){
		assert(XX_x == NULL && XY_x == NULL);
		var_x = Y2_x / num_x;
	}else{
		double sos = 0;
		for(int m=0; m<nFeatures; m++){
			// estimate B_{k,m}, b_{k,m}
			XX_inv = XX_x[m];
			XX_inv.addDiagonal(&lambda,1);
			XX_inv.sym_invert(debug);
			eta.product(XX_inv, XY_x[m]);
			b(m,0) = eta(0,0);
			for(int a=0; a<nLocalFactors; a++) B(m,a) = eta(a+1,0);
			// compute sos
			temp.product(XX_x[m], eta);
			temp2.transpose(eta); eta = temp2;
			temp2.product(eta, temp);
			temp.product(eta, XY_x[m]);
			temp.scale(2);
			temp2.subtract(temp);
			assert(temp2.length() == 1);
			sos += temp2(0,0);
		}
		sos += Y2_x;
		var_x = sos / num_x;
	}
}

void run_M_step(
	// Parameters: OUTPUT
	ListOfMatrices& A,      ListOfMatrices& B,    ListOfMatrices& b,
	ListOfMatrices& alpha,  ListOfMatrices& beta,
	double var_x[],  double var_y[],  double var_z[],  // vectors of length: nApps
	double& var_u, // scalar
	// Posterior mean of factors: INPUT
	const Matrix_ColumnMajor& u,  const ListOfMatrices& z,
	// Sufficient stats: INPUT
	const double sos_u, const double num_u,
	const Matrix_ColumnMajor  *XX_z, const Matrix_ColumnMajor  *XY_z, const double Y2_z[], const int num_z[],
	Matrix_ColumnMajor **XX_x, Matrix_ColumnMajor **XY_x, const double Y2_x[], const int num_x[],
	Matrix_ColumnMajor **XX_y, Matrix_ColumnMajor **XY_y, const double Y2_y[], const int num_y[],
	// Observation tables: INPUT
	const ObsTable& feature, const ObsTable& response,
	// Size information: INPUT
	const int nApps, const int nUsers, const int nGlobalFactors,
	const int nFeatures[], const int nItems[], const int nLocalFactors[], // vectors of length: nApps
	// Others
	const double lambda_A, const double lambda_B, const double lambda_beta,
	const int option,  const int verbose,  const int debug
){
	Matrix_ColumnMajor XX_inv, temp, temp2, eta, XY_m;

	// Estimate var_u
	var_u = sos_u / num_u;

	for(int k=0; k<nApps; k++){
		// Estimate A_k and var_{z,k}
		if(isIdentity(A.get(k))){
			assert(XX_z[k].length() == 0 && XY_z[k].length() == 0);
			var_z[k] = Y2_z[k] / num_z[k];
			if(verbose >= 5) printf("   A_%d = 1 and var_z[%d] = %f\n", k, k, var_z[k]);
		}else{
			// estimate A_k
			if(verbose >= 5) printf("   Estimate A_%d\n", k);
			XX_inv = XX_z[k];
			XX_inv.addDiagonal(&lambda_A, 1);
			XX_inv.sym_invert(debug);
			temp.product(XX_inv, XY_z[k]);
			A.get(k).transpose(temp);
			// estimate var_{z,k}
			if(verbose >= 5) printf("   Estimate var_{z,%d}", k);
			double sos = 0;
			XY_m.resize(nGlobalFactors,1);
			for(int m=0; m<nLocalFactors[k]; m++){
				eta.resize(nGlobalFactors, 1);
				for(int a=0; a<nGlobalFactors; a++)  eta(a,0) = A.get(k)(m,a);
				for(int a=0; a<nGlobalFactors; a++) XY_m(a,0) = XY_z[k](a,m);
				temp.product(XX_z[k], eta);
				temp2.transpose(eta); eta = temp2;
				temp2.product(eta, temp);
				temp.product(eta, XY_m);
				temp.scale(2);
				temp2.subtract(temp);
				assert(temp2.length() == 1);
				sos += temp2(0,0);
			}
			sos += Y2_z[k];
			var_z[k] = sos / num_z[k];
			if(verbose >= 5) printf(" = %f\n", var_z[k]);
		}
		// Estimate B_k, b_k, var_{x,k}
		estimateParamVar_oneApp(
			B.get(k), b.get(k), var_x[k],
			XX_x[k], XY_x[k], Y2_x[k], num_x[k],
			lambda_B, nFeatures[k], nLocalFactors[k], verbose, debug
		);
		// Estimate beta_k, alpha_k, var_{y,k}
		estimateParamVar_oneApp(
			beta.get(k), alpha.get(k), var_y[k],
			XX_y[k], XY_y[k], Y2_y[k], num_y[k],
			lambda_beta, nItems[k], nLocalFactors[k], verbose, debug
		);
	}
}


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
){
	// --------------------------------------------------------------
	// Prepare the space for the sufficient stats
	//    Computed in the E-step and used in the M-step
	// --------------------------------------------------------------
	if(verbose >= 2) printf("Prepare the space for the sufficient statistics\n");
	Matrix_ColumnMajor *XX_z = new Matrix_ColumnMajor[nApps];
	Matrix_ColumnMajor *XY_z = new Matrix_ColumnMajor[nApps];
	double *Y2_z = new double[nApps];
	int   *num_z = new int[nApps];
	for(int k=0; k<nApps; k++){
		if(isIdentity(A.get(k))){
			XX_z[k].resize(0,0);
			XY_z[k].resize(0,0);
		}else{
			XX_z[k].resize(nGlobalFactors, nGlobalFactors);
			XY_z[k].resize(nGlobalFactors, nLocalFactors[k]);
		}
	}

	Matrix_ColumnMajor **XX_x = new Matrix_ColumnMajor*[nApps];
	Matrix_ColumnMajor **XY_x = new Matrix_ColumnMajor*[nApps];
	double *Y2_x = new double[nApps];
	int   *num_x = new int[nApps];
	for(int k=0; k<nApps; k++){
		if(isIdentity(B.get(k))){
			XX_x[k] = NULL;
			XY_x[k] = NULL;
		}else{
			XX_x[k] = new Matrix_ColumnMajor[nFeatures[k]];
			XY_x[k] = new Matrix_ColumnMajor[nFeatures[k]];
			for(int m=0; m<nFeatures[k]; m++) XX_x[k][m].resize(nLocalFactors[k]+1, nLocalFactors[k]+1);
			for(int m=0; m<nFeatures[k]; m++) XY_x[k][m].resize(nLocalFactors[k]+1, 1);
		}
	}

	Matrix_ColumnMajor **XX_y = new Matrix_ColumnMajor*[nApps];
	Matrix_ColumnMajor **XY_y = new Matrix_ColumnMajor*[nApps];
	double *Y2_y = new double[nApps];
	int   *num_y = new int[nApps];
	for(int k=0; k<nApps; k++){
		XX_y[k] = new Matrix_ColumnMajor[nItems[k]];
		XY_y[k] = new Matrix_ColumnMajor[nItems[k]];
		for(int m=0; m<nItems[k]; m++) XX_y[k][m].resize(nLocalFactors[k]+1, nLocalFactors[k]+1);
		for(int m=0; m<nItems[k]; m++) XY_y[k][m].resize(nLocalFactors[k]+1, 1);
	}

	double sos_u, num_u;

	if(verbose >= 2) printf("E-step starts\n");
	run_E_step(
		// Posterior mean of factors: OUTPUT
		u, z,
		// Sufficient stats: OUTPUT
		sos_u, num_u,
		XX_z, XY_z, Y2_z, num_z,
		XX_x, XY_x, Y2_x, num_x,
		XX_y, XY_y, Y2_y, num_y,
		// Posterior variance for logistic regression
		xMean_var, xMean_var_length,
		yMean_var, yMean_var_length,
		// Parameters: INPUT
		A, B, b, alpha, beta, var_x, var_y, var_z, var_u,
		// Observation tables: INPUT
		feature, response,
		// Size information: INPUT
		nApps, nUsers, nGlobalFactors, nFeatures, nItems, nLocalFactors,
		// Others
		option, verbose, debug
	);

	if(verbose >= 2) printf("M-step starts\n");
	run_M_step(
		// Parameters: OUTPUT
		A, B, b, alpha, beta, var_x, var_y, var_z, var_u,
		// Posterior mean of factors: INPUT
		u, z,
		// Sufficient stats: INPUT
		sos_u, num_u,
		XX_z, XY_z, Y2_z, num_z,
		XX_x, XY_x, Y2_x, num_x,
		XX_y, XY_y, Y2_y, num_y,
		// Observation tables: INPUT
		feature, response,
		// Size information: INPUT
		nApps, nUsers, nGlobalFactors, nFeatures, nItems, nLocalFactors,
		// Others
		lambda_A, lambda_B, lambda_beta,
		option, verbose, debug
	);
	// printf("A = \n"); A.print();

	delete[] XX_z;
	delete[] XY_z;
	for(int k=0; k<nApps; k++){
		if(XX_x[k] != NULL) delete[] XX_x[k];  if(XY_x[k] != NULL) delete[] XY_x[k];
		if(XX_y[k] != NULL) delete[] XX_y[k];  if(XY_y[k] != NULL) delete[] XY_y[k];
	}
	delete[] XX_x;   delete[] XY_x;
	delete[] XX_y;   delete[] XY_y;
	delete[] Y2_z;   delete[] num_z;
	delete[] Y2_x;   delete[] num_x;
	delete[] Y2_y;   delete[] num_y;
}


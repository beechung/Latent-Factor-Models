/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <time.h>
#include "utils.hpp"

// u is used to denote this  effect
// v                   other effect
// If nThisContexts == 1, then nGlobalFactors and Q_size must be 0 and
//     Q, globalSample, globalPostMean, globalPostVar, globalPriorVar and thisContext will not be touched.
// If nOtherContexts == 1, then otherContext will not be touched.
extern "C" void gaussianPosterior_2WayInteraction_2Levels(
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
){
	int nObs=(*nObs_), nObsVar=(*nObsVar_), nLevelsThisEff=(*nLevelsThisEff_),
		nThisContexts=(*nThisContexts_), nLevelsOtherEff=(*nLevelsOtherEff_),
		nOtherContexts=(*nOtherContexts_), nLocalFactors=(*nLocalFactors_),
		nGlobalFactors=(*nGlobalFactors_), nOffset=(*nOffset_),
		Q_size=(*Q_size_), debug=(*debug_), verbose=(*verbose_);

	if(nObsVar != 1 && nObsVar != nObs) STOP1("nObsVar = %d", nObsVar);
	if(nThisContexts == 1 && nGlobalFactors != 0) STOP1("nThisContexts = 1, but nGlobalFactors = %d", nGlobalFactors);
	if(nThisContexts == 1 && Q_size != 0) STOP1("nThisContexts = 1, but Q_size = %d", Q_size);
	if(Q_size != 1 && Q_size != nLocalFactors*nGlobalFactors*nThisContexts) STOP1("Q_size = %d", Q_size);
	if(Q_size == 1 && Q[0] != 1 && Q[0] != 0) STOP1("Q = ", Q[0]);

	int outputMeanVar, drawSample;

	if(*option == 1){
        outputMeanVar = 0;  drawSample = 1;
    }else if(*option == 2){
        outputMeanVar = 1;  drawSample = 0;
    }else if(*option == 3){
        outputMeanVar = 1;  drawSample = 1;
    }else STOP1("Unknown option: %d", *option);

	// Initialize objects
	int max_nFactors = MAX(nLocalFactors, nGlobalFactors);
	int *num = (int*)Calloc(nThisContexts, int); // num[k]: #obs of user i in context k
	int N = (nThisContexts == 1 ? nLocalFactors : nGlobalFactors);
	Matrix_ColumnMajor temp1(max_nFactors, max_nFactors),
					   temp2(max_nFactors, max_nFactors),
					   temp3(max_nFactors, max_nFactors),
			           globalMean(N,1), globalVar(N,N);
	Matrix_ColumnMajor *E_i = new Matrix_ColumnMajor[nThisContexts];
	Matrix_ColumnMajor *S_i = new Matrix_ColumnMajor[nThisContexts];
	Matrix_ColumnMajor *z_i = new Matrix_ColumnMajor[nThisContexts];
	Matrix_ColumnMajor *Q_  = new Matrix_ColumnMajor[nThisContexts];
	Matrix_ColumnMajor *Qt_ = new Matrix_ColumnMajor[nThisContexts];
	for(int k=0; k<nThisContexts; k++){
		E_i[k].resize(nLocalFactors, nLocalFactors);
		S_i[k].resize(nLocalFactors, nLocalFactors);
		z_i[k].resize(nLocalFactors, 1);
		if(Q_size > 1){
			Q_[k].resize(nLocalFactors, nGlobalFactors);
			for(int a=0; a<nLocalFactors; a++) for(int b=0; b<nGlobalFactors; b++)
				Q_[k](a,b) = Q[C_3DA(a,b,k,nLocalFactors,nGlobalFactors)];
			Qt_[k].transpose(Q_[k]);
		}
	}
	double *v_jk = (double*)Calloc(nLocalFactors, double);
	// o[ijk] = obs[ijk] - v[j,,k]' offset[i,,k]
	double *o;
	if(nOffset == 0){
		o = (double*)obs;
	}else if(nOffset == nLevelsThisEff*nLocalFactors*nThisContexts){
		o = (double*)Calloc(nObs, double);
		computeMultiResponseUV(
			o, offset, otherEff, thisEffIndex, otherEffIndex, thisContext, otherContext,
			nObs_, nLocalFactors_, nLevelsThisEff_, nThisContexts_,
			nLevelsOtherEff_, nOtherContexts_, debug_
		);
		for(int m=0; m<nObs; m++) o[m] = obs[m] - o[m];
	}else STOP1("nOffset = %d", nOffset);
	const double one=1;  const int zero=0;
	int workspace_size = workspace_size_for_draw_multivar_gaussian(max_nFactors);
	double *workspace  = (double*)Calloc(workspace_size, double);
	double *sample     = (double*)Calloc(max_nFactors, double);
	double *sample2    = (double*)Calloc(max_nFactors, double);
	double *eigen_val  = (double*)Calloc(max_nFactors, double);
	double *eigen_vec  = (double*)Calloc(max_nFactors*max_nFactors, double);
	double *temp_a     = (double*)Calloc(max_nFactors, double);
	double *temp_b     = (double*)Calloc(max_nFactors*max_nFactors, double);
	int *is_Q_zero     = (int*)Calloc(nThisContexts, int);
	int *is_Q_identity = (int*)Calloc(nThisContexts, int);
	if(nThisContexts == 1){
		is_Q_zero[0] = 0;
		is_Q_identity[0] = 1;
	}else{
		for(int k=0; k<nThisContexts; k++){
			if(Q_size == 1){
				if(Q[0] == 1){
					is_Q_zero[k] = 0;
					is_Q_identity[k] = 1;
				}else if(Q[0] == 0){
					is_Q_zero[k] = 1;
					is_Q_identity[k] = 0;
				}else STOP1("Q[0] = ", Q[0]);
				continue;
			}
			is_Q_zero[k] = 1;
			is_Q_identity[k] = (nLocalFactors == nGlobalFactors ? 1 : 0);
			for(int a=0; a<nLocalFactors; a++){ for(int b=0; b<nGlobalFactors; b++){
					double Q_ab = Q[C_3DA(a,b,k,nLocalFactors,nGlobalFactors)];
					if(Q_ab != 0) is_Q_zero[k] = 0;
					if(a == b){ if(Q_ab != 1) is_Q_identity[k] = 0; }
					else      { if(Q_ab != 0) is_Q_identity[k] = 0; }
			} }
		}
	}
	for(int k=0; k<nThisContexts && nThisContexts > 1; k++){
		if(is_Q_identity[k] && nGlobalFactors != nLocalFactors) STOP("Q is identity, but nGlobalFactors != nLocalFactors");
	}

	double var_global = (nThisContexts == 1 ? localPriorVar[0] : globalPriorVar[0]);

	GetRNGstate();

	for(int i=0; i<nLevelsThisEff; i++){
		if(verbose >= 10) printf("process user %d:\n", i);
		// For each user i,
		//   S_i[k] = sum_j  (v[j,,k] v[j,,k]') / obsVar[ijk]
		//   z_i[k] = sum_j  v[j,,k] o[ijk] / obsVar[ijk]
		//   E_i[k] = (I/localPriorVar[k] + S_i[k])^-1
		for(int k=0; k<nThisContexts; k++){
			S_i[k].setToZero(); z_i[k].setToZero(); num[k] = 0;
		}
		int start_i = oiStart[i]-1;
		for(int q=0; q<oiNum[i]; q++){
			// For each observation of user i
			int oIndex = obsIndex[start_i+q] - 1; // convert to C index
			if(oIndex < 0 || oIndex > nObs) STOP4("oIndex = %d, i=%d, start_i=%d, num_i=%d", oIndex, i, start_i, oiNum[i]);
			double o_ijk = o[oIndex];
			int j  = otherEffIndex[oIndex]-1;
			int k1 = (nThisContexts ==1 ? 0 : thisContext[oIndex]-1);
			int k2 = (nOtherContexts==1 ? 0 : otherContext[oIndex]-1);
			double var_ijk = (nObsVar == 1 ? obsVar[0] : obsVar[oIndex]);
			if(debug > 0){
				if(thisEffIndex[oIndex]-1 != i) STOP3("thisEffIndex[%d]-1=%d, i=%d", oIndex, thisEffIndex[oIndex]-1, i);
				CHK_C_INDEX(j,nLevelsOtherEff); CHK_C_INDEX(k1,nThisContexts); CHK_C_INDEX(k2,nOtherContexts);
			}
			for(int a=0; a<nLocalFactors; a++) v_jk[a] = otherEff[C_3DA(j,a,k2,nLevelsOtherEff,nLocalFactors)];
			for(int a=0; a<nLocalFactors; a++){
				z_i[k1](a,0) += v_jk[a] * o_ijk / var_ijk;
				for(int b=0; b<nLocalFactors; b++)
					S_i[k1](a,b) += v_jk[a] * v_jk[b] / var_ijk;
			}
			num[k1]++;
		}
		// Now, S_i[k] = sum_j  (v[j,,k] v[j,,k]') / obsVar[ijk]
		//      z_i[k] = sum_j  v[j,,k] o[ijk] / obsVar[ijk]
		//      num[k] = number of obs of user i in context k
		globalMean.setToZero(); globalVar.setToZero();
		for(int k=0; k<nThisContexts; k++){
			// At the end,
			//    temp1 = Q_k' (I - S_i[k] E_i[k]) S_i[k] Q_k
			//    temp2 = Q_k' (I - S_i[k] E_i[k]) z_i[k]
			//   E_i[k] = (I/localPriorVar[k] + S_i[k])^-1
			if(localPriorVar[k] == 0 || nThisContexts == 1){
				E_i[k].setToZero();
				if(is_Q_zero[k]) continue;
				temp1 = S_i[k];
				temp2 = z_i[k];
			}else if(localPriorVar[k] > 0){
				if(num[k] == 0){
					E_i[k].setToZero();
					E_i[k].addDiagonal(localPriorVar[k]);
				}else{
					E_i[k] = S_i[k];
					E_i[k].addDiagonal(1.0/localPriorVar[k]);
					E_i[k].sym_invert(debug);
					// E_i[k] is exactly symmetric at this point
				}
				if(is_Q_zero[k]) continue;
				temp3.product(S_i[k], E_i[k]);
				temp3.negate();
				temp3.addDiagonal(1.0);       // temp3 = (I - S_i[k] E_i[k])
				temp1.product(temp3, S_i[k]); // temp1 = (I - S_i[k] E_i[k]) S_i[k]
				temp1.symmetrize(debug, 1e-10, 1e-6);
				temp2.product(temp3, z_i[k]); // temp2 = (I - S_i[k] E_i[k]) z_i[k]
			}else STOP("variance < 0");
			if(!is_Q_identity[k]){
				temp3.product(Qt_[k], temp1); // temp3 = Q_k' (I - S_i[k] E_i[k]) S_i[k]
				temp1.product(temp3, Q_[k]);  // temp1 = Q_k' (I - S_i[k] E_i[k]) S_i[k] Q_k
				temp1.symmetrize(debug, 1e-10, 1e-6);
				temp3 = temp2;                // temp3 = (I - S_i[k] E_i[k]) z_i[k]
				temp2.product(Qt_[k], temp3); // temp2 = Q_k' (I - S_i[k] E_i[k]) z_i[k]
			}
			globalMean.add(temp2);
			globalVar.add(temp1);
		}

		// Now, E_i[k] = (I/localPriorVar[k] + S_i[k])^-1
		//      globalMean = sum_k Q_k' (I - S_i[k] E_i[k]) z_i[k]
		//      globalVar  = sum_k Q_k' (I - S_i[k] E_i[k]) S_i[k] Q_k
		globalVar.addDiagonal(1.0/var_global);
		globalVar.sym_invert(debug);
		temp1 = globalMean;
		globalMean.product(globalVar, temp1);
		// Now, globalMean =   E[u_{i0} | o_i]
		//      globalVar  = Var[u_{i0} | o_i]

		if(nThisContexts == 1 && nOffset != 0){
			// special case: The mean and variance to be computed are actually for the local factor vector
			for(int a=0; a<nLocalFactors; a++)
				globalMean(a,0) += offset[C_3DA(i,a,0,nLevelsThisEff,nLocalFactors)];
		}

		// Draw a sample for u_{i0}, the global factor vector
		draw_multivar_gaussian(
			// OUTPUT
			sample,
			// INPUT/OUTPUT
			eigen_val, eigen_vec,
			// INPUT
			globalMean.getData(), globalVar.getData(),
			&N, &zero, &debug,
			// WORK/TEMP SPACE
			workspace, &workspace_size, temp_a, temp_b
		);
		if(verbose >= 15){
			print_vector("   eigen value = ", eigen_val, N);
			printf("   eigen vectors =\n");
			print_matrix("      ", eigen_vec, N, N);
		}
		if(drawSample){
			if(nThisContexts == 1){
				// special case: The sample drawn is actually the local factor vector
				for(int a=0; a<nLocalFactors; a++) localSample[C_3DA(i,a,0,nLevelsThisEff,nLocalFactors)] = sample[a];
			}else{
				for(int a=0; a<nGlobalFactors; a++) globalSample[C_MAT(i,a,nLevelsThisEff)] = sample[a];
			}
		}
		if(outputMeanVar){
			if(nThisContexts == 1){
				// special case: The mean and variance computed are actually for the local factor vector
				for(int a=0; a<nLocalFactors; a++){
					localPostMean[C_3DA(i,a,0,nLevelsThisEff,nLocalFactors)] = globalMean(a,0);
					for(int b=0; b<nLocalFactors; b++)
						localPostVar[C_4DA(i,a,b,0,nLevelsThisEff,nLocalFactors,nLocalFactors)] = globalVar(a,b);
				}
			}else{
				for(int a=0; a<nGlobalFactors; a++){
					globalPostMean[C_MAT(i,a,nLevelsThisEff)] = globalMean(a,0);
					for(int b=0; b<nGlobalFactors; b++)
						globalPostVar[C_3DA(i,a,b,nLevelsThisEff,nGlobalFactors)] = globalVar(a,b);
				}
			}
		}
		if(nThisContexts == 1) continue;

		for(int k=0; k<nThisContexts; k++){
			// Draw a sample for u_{ik}, the local factor vector in context k
			if(localPriorVar[k] > 0){
				if(is_Q_zero[k]){
					temp1.resize(nLocalFactors,1);
					temp1.setToZero();
				}else if(!is_Q_identity[k]){
					temp2.resize(nGlobalFactors,1);
					for(int a=0; a<nGlobalFactors; a++) temp2(a,0) = sample[a] / localPriorVar[k];
					temp1.product(Q_[k], temp2);
				}else{
					temp1.resize(nLocalFactors,1);
					for(int a=0; a<nLocalFactors; a++) temp1(a,0) = sample[a] / localPriorVar[k];
				}
				// Now, temp1 = (Q_k u_{i0}) / localPriorVar[k]
				if(num[k] > 0) temp1.add(z_i[k]);
				// Now, temp1 = (Q_k u_{i0}) / localPriorVar[k] + z_{ik}
				temp2.product(E_i[k], temp1);
				if(nOffset != 0){
					for(int a=0; a<nLocalFactors; a++)
						temp2(a,0) += offset[C_3DA(i,a,k,nLevelsThisEff,nLocalFactors)];
				}
				// Now, temp2 =   E[u_{ik} | u_{i0}, o_i]
				//     E_i[k] = Var[u_{ik} | u_{i0}, o_i]
				draw_multivar_gaussian(
					// OUTPUT
					sample2,
					// INPUT/OUTPUT
					eigen_val, eigen_vec,
					// INPUT
					temp2.getData(), E_i[k].getData(),
					&nLocalFactors, &zero, &debug,
					// WORK/TEMP SPACE
					workspace, &workspace_size, temp_a, temp_b
				);
			}else{
				// localPriorVar[k] == 0
				if(is_Q_zero[k]){
					temp2.resize(nLocalFactors,1);
					temp2.setToZero();
				}else if(!is_Q_identity[k]){
					temp1.resize(nLocalFactors,1);
					for(int a=0; a<nLocalFactors; a++) temp1(a,0) = sample[a];
					temp2.product(Q_[k], temp1);
				}else{
					temp2.resize(nLocalFactors,1);
					for(int a=0; a<nLocalFactors; a++) temp2(a,0) = sample[a];
				}
				if(nOffset != 0){
					for(int a=0; a<nLocalFactors; a++)
						temp2(a,0) += offset[C_3DA(i,a,k,nLevelsThisEff,nLocalFactors)];
				}
				// Now, temp2 =   E[u_{ik} | u_{i0}, o_i]
				//     E_i[k] = Var[u_{ik} | u_{i0}, o_i]
				for(int a=0; a<nLocalFactors; a++) sample2[a] = temp2(a,0);
				if(debug){
					// make sure the random numbers drawn are the same as the R testing code
					for(int a=0; a<nLocalFactors; a++) norm_rand();
				}
			}
			if(drawSample){
				for(int a=0; a<nLocalFactors; a++) localSample[C_3DA(i,a,k,nLevelsThisEff,nLocalFactors)] = sample2[a];
			}
			if(outputMeanVar){
				for(int a=0; a<nLocalFactors; a++){
					localPostMean[C_3DA(i,a,k,nLevelsThisEff,nLocalFactors)] = temp2(a,0);
					for(int b=0; b<nLocalFactors; b++)
						localPostVar[C_4DA(i,a,b,k,nLevelsThisEff,nLocalFactors,nLocalFactors)] = E_i[k](a,b);
				}
			}
		}
	}

	PutRNGstate();

	Free(num);  Free(v_jk);  Free(workspace);  Free(sample);  Free(sample2);
	Free(eigen_val);  Free(eigen_vec);  Free(temp_a);  Free(temp_b);
	Free(is_Q_zero);  Free(is_Q_identity);
	delete[] E_i;  delete[] z_i;  delete[] S_i;
	delete[] Q_;   delete[] Qt_;
	if(nOffset != 0) Free(o);
}


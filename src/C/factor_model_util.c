/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#include <time.h>
#include "util.h"


void gaussianPosterior_mainEffect(
    // OUTPUT
    double* outSample, double* outMean, double* outVar,
    //INPUT
    const int* option /*1:Sample, 2:Mean&Var, 3:Sample&Mean&Var*/,
    const int* thisEffIndex /*user or item*/,
    const double* rest /*o in the paper*/,
    const double* fittedEff /*g0*x_user or d0*x_item*/,
    const double* multiplier /*NULL*/,
    const double* var_y, const double* var_eff /*var_alpha or var_beta*/,
    const int* nObs, const int* nThisEff, const int* nVar_y, const int* nVar_eff,
    // OTHER
    const int* debug
){
    int outputSample=0, outputMeanVar=0;

    if(*option == 1){
        outputSample = 1; outputMeanVar = 0;
    }else if(*option == 2){
        outputSample = 0; outputMeanVar = 1;
    }else if(*option == 3){
        outputSample = 1; outputMeanVar = 1;
    }else error("Unknown option: %d", *option);

    double *sum_ivar = (double*)Calloc(*nThisEff, double); // sum_ivar[effect_i] = sum_{k having effect_i} 1/var_y[k]
    double *sum_o    = (double*)Calloc(*nThisEff, double); // sum_o[effect_i] = sum_{k having effect_i} o[k]/var_y[k]

    for(int k=0; k<*nObs; k++){
        const int thisIndex  = thisEffIndex[k]-1;
        if(*debug > 0) CHK_C_INDEX(thisIndex,  *nThisEff);

        double o = rest[k]; // for alpha and beta: o.  for gamma: o * x_dyad * b
        double square = 1;  // for alpha and beta: 1.  for gamma: (x_dyad * b)^2
        if(multiplier != NULL){
            o *= multiplier[k];
            square = multiplier[k] * multiplier[k];
        }

        if((*nVar_y) == 1){
            sum_ivar[thisIndex] += square / var_y[0];
            sum_o[   thisIndex] +=      o / var_y[0];
        }else if((*nVar_y) == (*nObs)){
            sum_ivar[thisIndex] += square / var_y[k];
            sum_o[   thisIndex] +=      o / var_y[k];
        }else error("nVar_y = %d, nObs = %d", *nVar_y, *nObs);
    }

    if(outputSample) GetRNGstate();
    for(int i=0; i<*nThisEff; i++){
    	double mean = 0, var=0;
        if((*nVar_eff) == 1){
            var  = 1.0 / (sum_ivar[i] + (1.0 / var_eff[0]));
            mean = var * (sum_o[i] + (fittedEff[i] / var_eff[0]));
        }else if((*nVar_eff) == (*nThisEff)){
            var  = 1.0 / (sum_ivar[i] + (1.0 / var_eff[i]));
            mean = var * (sum_o[i] + (fittedEff[i] / var_eff[i]));
        }else error("nVar_eff = %d, nThisEff = %d", *nVar_eff, *nThisEff);

        if(outputMeanVar){
            outMean[i] = mean;
            outVar[i]  = var;
        }
        if(outputSample){
            outSample[i] = rnorm(mean, sqrt(var));
        }
    }
    if(outputSample) PutRNGstate();

    Free(sum_ivar);
    Free(sum_o);
}


void gaussianPosterior_2WayInteraction(
    // OUTPUT
    double* sample, double* posteriorMean, double* posteriorVar,
    // INPUT
    const int* option /*1:Sample, 2:Mean&Var, 3:Sample&Mean&Var*/,
    const int* thisEffIndex /*user or item*/, const int* otherEffIndex /*item or user*/, const double* obs /*o in the paper*/,
    const double* priorMean /*Gw or Dz*/, const double* otherEff /*v or u*/,
    const double* obsVar, const double* priorVar /*var_u or var_v*/,
    const int* nObs, const int* nLevelsThisEff, const int* nLevelsOtherEff,
    const int* nFactors, const int* nObsVar, const int *nPriorVar,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    // OTHER
    const int* debug
){
    double *sum_vv /*nFactors x nFactors*/, *sum_ov /*nFactors x 1*/, o, *vj /*nFactors x 1*/, *Gwi /*nFactors x 1*/,
           *work, size, *mean, *var, *temp, *rnd;
    int i,j,k,m, thisIndex, otherIndex, outputSample, outputMeanVar, oIndex, lwork=-1, info, symCheck;
    char jobz = 'V', uplo = 'L';

    symCheck = (*debug) - 2;

    if(*option == 1){
        outputSample = 1; outputMeanVar = 0;
    }else if(*option == 2){
        outputSample = 0; outputMeanVar = 1;
    }else if(*option == 3){
        outputSample = 1; outputMeanVar = 1;
    }else error("Unknown option: %d", *option);

    vj     = (double*)Calloc(*nFactors, double);
    Gwi    = (double*)Calloc(*nFactors, double);
    sum_ov = (double*)Calloc(*nFactors, double);
    sum_vv = (double*)Calloc((*nFactors)*(*nFactors), double);
    mean   = (double*)Calloc(*nFactors, double);
    var    = (double*)Calloc((*nFactors)*(*nFactors), double);
    temp   = (double*)Calloc((*nFactors)*(*nFactors), double);

    if(outputSample) GetRNGstate();

    for(i=0; i<*nLevelsThisEff; i++){

        thisIndex = i+1;
        for(j=0; j<*nFactors; j++) sum_ov[j] = 0;
        for(j=0; j<(*nFactors)*(*nFactors); j++) sum_vv[j] = 0;

        for(j=0; j<oiNum[i]; j++){

            oIndex = obsIndex[R_VEC(oiStart[i]+j)];
            otherIndex = otherEffIndex[R_VEC(oIndex)];

            if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
            if(*debug > 0) CHK_R_INDEX(otherIndex, *nLevelsOtherEff);
            if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");

            o = obs[R_VEC(oIndex)];
            for(k=1; k<=*nFactors; k++) vj[R_VEC(k)] = otherEff[R_MAT(otherIndex,k,*nLevelsOtherEff)];

            double var_y_thisObs = 0;
            if((*nObsVar) == 1)            var_y_thisObs = obsVar[0];
            else if((*nObsVar) == (*nObs)) var_y_thisObs = obsVar[R_VEC(oIndex)];
            else error("nVar_y = %d, nObs = %d", *nObsVar, *nObs);

            for(k=0; k<*nFactors; k++) sum_ov[k] += (o * vj[k]) / var_y_thisObs;
            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += (vj[k] * vj[m]) / var_y_thisObs;
        }

        for(k=1; k<=*nFactors; k++) Gwi[R_VEC(k)] = priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)];

        if((*nPriorVar) == 1){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[0]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[0];
        }else if((*nPriorVar) == (*nLevelsThisEff)*(*nFactors)*(*nFactors)){

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) temp[C_MAT(k,m,*nFactors)] = priorVar[C_3DA(i,k,m,*nLevelsThisEff,*nFactors)];
            sym_inv_byCholesky(temp, nFactors, &symCheck);
            // Now, temp is var_u[i]^-1

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += temp[C_MAT(k,m,*nFactors)];

            // reuse the space of vj
            for(k=0; k<*nFactors; k++) vj[k] = Gwi[k];
            for(k=0; k<*nFactors; k++){
                Gwi[k] = 0;
                for(m=0; m<*nFactors; m++) Gwi[k] += temp[C_MAT(k,m,*nFactors)] * vj[m];
            }

        }else error("nVar_eff = %d, nThisEff = %d", *nPriorVar, *nLevelsThisEff);

        // Now, sum_vv = var^-1
        //      Gwi = var_u[i]^-1 G wi

        if((*nFactors) == 1){
            var[0]  = 1/sum_vv[0];
            mean[0] = var[0] * (sum_ov[0] + Gwi[0]);
            if(outputSample){
                sample[i] = rnorm(mean[0], sqrt(var[0]));
            }
            if(outputMeanVar){
                posteriorMean[i] = mean[0];
                posteriorVar[i]  = var[0];
            }
        }else{
            double *eigen_val, *eigen_vec;
            eigen_val = vj;
            eigen_vec = sum_vv;

            if(*debug > 2) CHK_SYMMETRIC(eigen_vec, *nFactors);

            //
            // Compute the variance-covariance matrix
            //
            // Allocate workspace for eigen value decomposition (only allocate once)
            if(lwork == -1){
                F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, &size, &lwork, &info);
                if(info != 0) error("error in dsyev(...)");
                lwork = (int)size;
                work  = (double*)Calloc(lwork, double);
            }

            F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, work, &lwork, &info);
            if(info != 0) error("error in dsyev(...)");
            for(j=0; j<*nFactors; j++) eigen_val[j] = 1/eigen_val[j];
            // Now, eigen_val, eigen_vec are the eigen values and vectors of var

            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];
            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++){
                    var[C_MAT(j,k,*nFactors)] = 0;
                    for(m=0; m<*nFactors; m++) var[C_MAT(j,k,*nFactors)] += temp[C_MAT(j,m,*nFactors)] * eigen_vec[C_MAT(k,m,*nFactors)];
                }
            // Now, var is the variance-covariance matrix

            //
            // Compute the mean vector
            //

            // print_vector("sum_ov: ", sum_ov, *nFactors);
            // print_vector("Gwi: ", Gwi, *nFactors);

            for(j=0; j<*nFactors; j++) sum_ov[j] += Gwi[j];
            for(j=0; j<*nFactors; j++){
                mean[j] = 0;
                for(k=0; k<*nFactors; k++) mean[j] += var[C_MAT(j,k,*nFactors)] * sum_ov[k];
            }

            if(outputMeanVar){
                for(j=0; j<*nFactors; j++) posteriorMean[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) posteriorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)] = var[C_MAT(j,k,*nFactors)];
            }

            // DEBUG CODE
            // if(i==38){
        	//     printf("i = %d\n",i);
            //     print_vector("mean: ", mean, *nFactors);
            //     print_matrix(" var: ", var, *nFactors, *nFactors);
            // }

            if(*debug >= 2 && oiNum[i] == 0){
                for(k=1; k<=*nFactors; k++){
                    CHK_SAME_NUMBER("mean[k] != fittedEff[k]", mean[R_VEC(k)], priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)]);
                }
                if((*nPriorVar) == 1){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[0]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else{
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)]);
                        }
                }
            }
            //
            //  Generate the random vector
            //
            if(outputSample){
                // reuse the space allocated for Gwi
                rnd = Gwi;

                // DEBUG CODE
                // if(i==38){
                //     print_vector("eigenvalue: ", eigen_val, *nFactors);
                //     Rprintf("eigenvector:\n");  print_matrix("    ", eigen_vec, *nFactors, *nFactors);
                // }

                for(j=0; j<*nFactors; j++) eigen_val[j] = sqrt(eigen_val[j]);
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];

                for(j=0; j<*nFactors; j++) rnd[j] = norm_rand();

                // DEBUG CODE
                // if(i==38){
                //     Rprintf("temp:\n");
                //     print_matrix("    ", temp, *nFactors, *nFactors);
                //     print_vector("rnd: ", rnd, *nFactors);
                // }

                for(j=0; j<*nFactors; j++){
                    sample[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                    for(k=0; k<*nFactors; k++) sample[C_MAT(i,j,*nLevelsThisEff)] += temp[C_MAT(j,k,*nFactors)] * rnd[k];
                }
            }
        }
    }
    if(outputSample) PutRNGstate();

    Free(sum_ov);
    Free(sum_vv);
    Free(vj);
    Free(Gwi);
    Free(mean);
    Free(var);
    Free(temp);
    if(lwork > 0) Free(work);
}


void gaussianPosterior_SelfInteraction(
	// IN/OUT
    double* sample /* v */,
    // OUTPUT
    double* posteriorMean, double* posteriorVar,
    // INPUT
    const int* option /*1:Sample only, 2:Sample&Mean&Var, 3:Sample&Mean&Var*/,
    const int* fromIndex, const int* toIndex, const double* obs /*o in the paper*/,
    const double* priorMean /*Gx_user*/,
    const double* obsVar, const double* priorVar /* var_v */,
    const int* nObs, const int* nLevelsThisEff, const int* nFactors,
    const int* nObsVar, const int *nPriorVar,
    const int* author_obsIndex, const int* author_oiStart, const int* author_oiNum,
    const int*  voter_obsIndex, const int*  voter_oiStart, const int*  voter_oiNum,
    // OTHER
    const int* debug
){
    double *sum_vv /*nFactors x nFactors*/, *sum_ov /*nFactors x 1*/, o, *vj /*nFactors x 1*/, *Gwi /*nFactors x 1*/,
           *work, size, *mean, *var, *temp, *rnd;
    int i,j,k,m, thisIndex, otherIndex, outputMeanVar, oIndex, lwork=-1, info, symCheck, role;
    char jobz = 'V', uplo = 'L';

    symCheck = (*debug) - 2;

    if(*option == 1){
        outputMeanVar = 0;
    }else if(*option == 2){
        outputMeanVar = 1;
    }else if(*option == 3){
        outputMeanVar = 1;
    }else error("Unknown option: %d", *option);

    vj     = (double*)Calloc(*nFactors, double);
    Gwi    = (double*)Calloc(*nFactors, double);
    sum_ov = (double*)Calloc(*nFactors, double);
    sum_vv = (double*)Calloc((*nFactors)*(*nFactors), double);
    mean   = (double*)Calloc(*nFactors, double);
    var    = (double*)Calloc((*nFactors)*(*nFactors), double);
    temp   = (double*)Calloc((*nFactors)*(*nFactors), double);

    GetRNGstate();

    for(i=0; i<*nLevelsThisEff; i++){

        thisIndex = i+1;
        for(j=0; j<*nFactors; j++) sum_ov[j] = 0;
        for(j=0; j<(*nFactors)*(*nFactors); j++) sum_vv[j] = 0;
        int nObs_for_i = 0;

        // begin: compute sum_ov and sum_vv
        for(role=0; role<2; role++){
        	int *obsIndex = NULL, *oiNum = NULL, *oiStart=NULL, *thisEffIndex=NULL, *otherEffIndex=NULL;
        	if(role == 0){
        		// user i is an author
        		obsIndex = (int*)author_obsIndex;  oiNum = (int*)author_oiNum;  oiStart = (int*)author_oiStart;
        		thisEffIndex = (int*)toIndex;  otherEffIndex = (int*)fromIndex;
        	}else if(role == 1){
        		// user i is an voter
        		obsIndex = (int*)voter_obsIndex;   oiNum = (int*)voter_oiNum;   oiStart = (int*)voter_oiStart;
        		thisEffIndex = (int*)fromIndex; otherEffIndex = (int*)toIndex;
        	}else DIE_HERE;

			for(j=0; j<oiNum[i]; j++){

				nObs_for_i++;

				oIndex = obsIndex[R_VEC(oiStart[i]+j)];
				otherIndex = otherEffIndex[R_VEC(oIndex)];
				if(otherIndex == i+1) error("self votes are not allowed");

				if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
				if(*debug > 0) CHK_R_INDEX(otherIndex, *nLevelsThisEff);
				if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");

				o = obs[R_VEC(oIndex)];
				for(k=1; k<=*nFactors; k++) vj[R_VEC(k)] = sample[R_MAT(otherIndex,k,*nLevelsThisEff)];

				double var_y_thisObs = 0;
				if((*nObsVar) == 1)            var_y_thisObs = obsVar[0];
				else if((*nObsVar) == (*nObs)) var_y_thisObs = obsVar[R_VEC(oIndex)];
				else error("nVar_y = %d, nObs = %d", *nObsVar, *nObs);

				for(k=0; k<*nFactors; k++) sum_ov[k] += (o * vj[k]) / var_y_thisObs;
				for(k=0; k<*nFactors; k++)
					for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += (vj[k] * vj[m]) / var_y_thisObs;
			}
        }
        // end: compute sum_ov and sum_vv

        for(k=1; k<=*nFactors; k++) Gwi[R_VEC(k)] = priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)];

        if((*nPriorVar) == 1){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[0]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[0];
        }else if((*nPriorVar) == (*nLevelsThisEff)*(*nFactors)*(*nFactors)){

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) temp[C_MAT(k,m,*nFactors)] = priorVar[C_3DA(i,k,m,*nLevelsThisEff,*nFactors)];
            sym_inv_byCholesky(temp, nFactors, &symCheck);
            // Now, temp is var_u[i]^-1

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += temp[C_MAT(k,m,*nFactors)];

            // reuse the space of vj
            for(k=0; k<*nFactors; k++) vj[k] = Gwi[k];
            for(k=0; k<*nFactors; k++){
                Gwi[k] = 0;
                for(m=0; m<*nFactors; m++) Gwi[k] += temp[C_MAT(k,m,*nFactors)] * vj[m];
            }

        }else error("nVar_eff = %d, nThisEff = %d", *nPriorVar, *nLevelsThisEff);

        // Now, sum_vv = var^-1
        //      Gwi = var_u[i]^-1 G wi

        if((*nFactors) == 1){
            var[0]  = 1/sum_vv[0];
            mean[0] = var[0] * (sum_ov[0] + Gwi[0]);
            sample[i] = rnorm(mean[0], sqrt(var[0]));
            if(outputMeanVar){
                posteriorMean[i] = mean[0];
                posteriorVar[i]  = var[0];
            }
        }else{
            double *eigen_val, *eigen_vec;
            eigen_val = vj;
            eigen_vec = sum_vv;

            if(*debug > 2) CHK_SYMMETRIC(eigen_vec, *nFactors);

            //
            // Compute the variance-covariance matrix
            //
            // Allocate workspace for eigen value decomposition (only allocate once)
            if(lwork == -1){
                F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, &size, &lwork, &info);
                if(info != 0) error("error in dsyev(...)");
                lwork = (int)size;
                work  = (double*)Calloc(lwork, double);
            }

            F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, work, &lwork, &info);
            if(info != 0) error("error in dsyev(...)");
            for(j=0; j<*nFactors; j++) eigen_val[j] = 1/eigen_val[j];
            // Now, eigen_val, eigen_vec are the eigen values and vectors of var

            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];
            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++){
                    var[C_MAT(j,k,*nFactors)] = 0;
                    for(m=0; m<*nFactors; m++) var[C_MAT(j,k,*nFactors)] += temp[C_MAT(j,m,*nFactors)] * eigen_vec[C_MAT(k,m,*nFactors)];
                }
            // Now, var is the variance-covariance matrix

            //
            // Compute the mean vector
            //

            // print_vector("sum_ov: ", sum_ov, *nFactors);
            // print_vector("Gwi: ", Gwi, *nFactors);

            for(j=0; j<*nFactors; j++) sum_ov[j] += Gwi[j];
            for(j=0; j<*nFactors; j++){
                mean[j] = 0;
                for(k=0; k<*nFactors; k++) mean[j] += var[C_MAT(j,k,*nFactors)] * sum_ov[k];
            }

            if(outputMeanVar){
                for(j=0; j<*nFactors; j++) posteriorMean[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) posteriorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)] = var[C_MAT(j,k,*nFactors)];
            }

            if(*debug >= 2 && nObs_for_i == 0){
                for(k=1; k<=*nFactors; k++){
                    CHK_SAME_NUMBER("mean[k] != fittedEff[k]", mean[R_VEC(k)], priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)]);
                }
                if((*nPriorVar) == 1){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[0]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else{
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)]);
                        }
                }
            }
            //
            //  Generate the random vector
            //
			// reuse the space allocated for Gwi
			rnd = Gwi;

			for(j=0; j<*nFactors; j++) eigen_val[j] = sqrt(eigen_val[j]);
			for(j=0; j<*nFactors; j++)
				for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];

			for(j=0; j<*nFactors; j++) rnd[j] = norm_rand();

			for(j=0; j<*nFactors; j++){
				sample[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
				for(k=0; k<*nFactors; k++) sample[C_MAT(i,j,*nLevelsThisEff)] += temp[C_MAT(j,k,*nFactors)] * rnd[k];
			}
        }
    }
    PutRNGstate();

    Free(sum_ov);  Free(sum_vv);  Free(vj);    Free(Gwi);
    Free(mean);    Free(var);     Free(temp);
    if(lwork > 0) Free(work);
}

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
    const double* otherEff /* nLevelsOtherEff x nFactors */,
    const double* thirdEff /* nLevelsThirdEff x nFactors */,
    const double* obsVar   /* nObsVar x 1 */,
    const double* priorVar /* nPriorVar x 1 */,
    const int* nObs, const int* nLevelsThisEff, const int* nLevelsOtherEff, const int* nLevelsThirdEff,
    const int* nFactors, const int* nObsVar   /* 1 or nObs */,
    const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
    const int* obsIndex, const int* oiStart, const int* oiNum,
    // OTHER
    const int* debug
){
    double *sum_vv /*nFactors x nFactors*/, *sum_ov /*nFactors x 1*/, o, *vj /*nFactors x 1*/, *Gwi /*nFactors x 1*/,
           *work, size, *mean, *var, *temp, *rnd;
    int i,j,k,m, thisIndex, otherIndex, outputSample, outputMeanVar, oIndex, lwork=-1, info, symCheck;
    char jobz = 'V', uplo = 'L';

    symCheck = (*debug) - 2;

    if(*option == 1){
        outputSample = 1; outputMeanVar = 0;
    }else if(*option == 2){
        outputSample = 0; outputMeanVar = 1;
    }else if(*option == 3){
        outputSample = 1; outputMeanVar = 1;
    }else error("Unknown option: %d", *option);

    vj     = (double*)Calloc(*nFactors, double);
    Gwi    = (double*)Calloc(*nFactors, double);
    sum_ov = (double*)Calloc(*nFactors, double);
    sum_vv = (double*)Calloc((*nFactors)*(*nFactors), double);
    mean   = (double*)Calloc(*nFactors, double);
    var    = (double*)Calloc((*nFactors)*(*nFactors), double);
    temp   = (double*)Calloc((*nFactors)*(*nFactors), double);

    if(outputSample) GetRNGstate();

    for(i=0; i<*nLevelsThisEff; i++){

        thisIndex = i+1;
        for(j=0; j<*nFactors; j++) sum_ov[j] = 0;
        for(j=0; j<(*nFactors)*(*nFactors); j++) sum_vv[j] = 0;

        for(j=0; j<oiNum[i]; j++){

            oIndex = obsIndex[R_VEC(oiStart[i]+j)];
            otherIndex = otherEffIndex[R_VEC(oIndex)];

            if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
            if(*debug > 0) CHK_R_INDEX(otherIndex, *nLevelsOtherEff);
            if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");

            o = obs[R_VEC(oIndex)];

            if((*nLevelsThirdEff) > 0){
                int thirdIndex = thirdEffIndex[R_VEC(oIndex)];
                if(*debug > 0) CHK_R_INDEX(thirdIndex, *nLevelsThirdEff);
            	for(k=1; k<=*nFactors; k++)
            		vj[R_VEC(k)] = otherEff[R_MAT(otherIndex,k,*nLevelsOtherEff)] *
            		               thirdEff[R_MAT(thirdIndex,k,*nLevelsThirdEff)];
            }else{
            	for(k=1; k<=*nFactors; k++) vj[R_VEC(k)] = otherEff[R_MAT(otherIndex,k,*nLevelsOtherEff)];
            }

            double var_y_thisObs = 0;
            if((*nObsVar) == 1)            var_y_thisObs = obsVar[0];
            else if((*nObsVar) == (*nObs)) var_y_thisObs = obsVar[R_VEC(oIndex)];
            else error("nVar_y = %d, nObs = %d", *nObsVar, *nObs);

            for(k=0; k<*nFactors; k++) sum_ov[k] += (o * vj[k]) / var_y_thisObs;
            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += (vj[k] * vj[m]) / var_y_thisObs;
        }

        for(k=1; k<=*nFactors; k++) Gwi[R_VEC(k)] = priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)];

        if((*nPriorVar) == 1){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[0]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[0];
        }else if((*nPriorVar) == (*nFactors)){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[k]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[k];
        }else if((*nPriorVar) == (*nLevelsThisEff)*(*nFactors)*(*nFactors)){

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) temp[C_MAT(k,m,*nFactors)] = priorVar[C_3DA(i,k,m,*nLevelsThisEff,*nFactors)];
            sym_inv_byCholesky(temp, nFactors, &symCheck);
            // Now, temp is var_u[i]^-1

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += temp[C_MAT(k,m,*nFactors)];

            // reuse the space of vj
            for(k=0; k<*nFactors; k++) vj[k] = Gwi[k];
            for(k=0; k<*nFactors; k++){
                Gwi[k] = 0;
                for(m=0; m<*nFactors; m++) Gwi[k] += temp[C_MAT(k,m,*nFactors)] * vj[m];
            }

        }else error("nVar_eff = %d, nThisEff = %d", *nPriorVar, *nLevelsThisEff);

        // Now, sum_vv = var^-1
        //      Gwi = var_u[i]^-1 G wi

        if((*nFactors) == 1){
            var[0]  = 1/sum_vv[0];
            mean[0] = var[0] * (sum_ov[0] + Gwi[0]);
            if(outputSample){
                sample[i] = rnorm(mean[0], sqrt(var[0]));
            }
            if(outputMeanVar){
                posteriorMean[i] = mean[0];
                posteriorVar[i]  = var[0];
            }
        }else{
            double *eigen_val, *eigen_vec;
            eigen_val = vj;
            eigen_vec = sum_vv;

            if(*debug > 2) CHK_SYMMETRIC(eigen_vec, *nFactors);

            //
            // Compute the variance-covariance matrix
            //
            // Allocate workspace for eigen value decomposition (only allocate once)
            if(lwork == -1){
                F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, &size, &lwork, &info);
                if(info != 0) error("error in dsyev(...)");
                lwork = (int)size;
                work  = (double*)Calloc(lwork, double);
            }

            F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, work, &lwork, &info);
            if(info != 0) error("error in dsyev(...)");
            for(j=0; j<*nFactors; j++) eigen_val[j] = 1/eigen_val[j];
            // Now, eigen_val, eigen_vec are the eigen values and vectors of var

            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];
            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++){
                    var[C_MAT(j,k,*nFactors)] = 0;
                    for(m=0; m<*nFactors; m++) var[C_MAT(j,k,*nFactors)] += temp[C_MAT(j,m,*nFactors)] * eigen_vec[C_MAT(k,m,*nFactors)];
                }
            // Now, var is the variance-covariance matrix

            //
            // Compute the mean vector
            //
            for(j=0; j<*nFactors; j++) sum_ov[j] += Gwi[j];
            for(j=0; j<*nFactors; j++){
                mean[j] = 0;
                for(k=0; k<*nFactors; k++) mean[j] += var[C_MAT(j,k,*nFactors)] * sum_ov[k];
            }

            if(outputMeanVar){
                for(j=0; j<*nFactors; j++) posteriorMean[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) posteriorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)] = var[C_MAT(j,k,*nFactors)];
            }

            if(*debug >= 2 && oiNum[i] == 0){
				// printf("i = %d\n",i);
				// print_vector("     mean: ", mean, *nFactors);
				// print_matrix("      var: ", var, *nFactors, *nFactors);
				// print_vector("prior_var: ", priorVar, *nPriorVar);
                for(k=1; k<=*nFactors; k++){
                    CHK_SAME_NUMBER("mean[k] != fittedEff[k]", mean[R_VEC(k)], priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)]);
                }
                if((*nPriorVar) == 1){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[0]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else if((*nPriorVar) == (*nFactors)){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[k]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else{
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)]);
                        }
                }
            }
            //
            //  Generate the random vector
            //
            if(outputSample){
                // reuse the space allocated for Gwi
                rnd = Gwi;

                for(j=0; j<*nFactors; j++) eigen_val[j] = sqrt(eigen_val[j]);
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];

                for(j=0; j<*nFactors; j++) rnd[j] = norm_rand();

                for(j=0; j<*nFactors; j++){
                    sample[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                    for(k=0; k<*nFactors; k++) sample[C_MAT(i,j,*nLevelsThisEff)] += temp[C_MAT(j,k,*nFactors)] * rnd[k];
                }
            }
        }
    }
    if(outputSample) PutRNGstate();

    Free(sum_ov);
    Free(sum_vv);
    Free(vj);
    Free(Gwi);
    Free(mean);
    Free(var);
    Free(temp);
    if(lwork > 0) Free(work);
}

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
	const double* thirdEff /* nLevelsThirdEff x nFactors */,
	const int* nObs, const int* nLevelsThisEff, const int* nLevelsThirdEff, const int* nFactors,
	const int* nObsVar   /* 1 or nObs */,
	const int *nPriorVar /* 1 or nLevelsThisEff*nFactors*nFactors */,
	const int* from_obsIndex, const int* from_oiStart, const int* from_oiNum,
	const int*   to_obsIndex, const int*   to_oiStart, const int*   to_oiNum,
	// OTHER
	const int* debug
){
    double *sum_vv /*nFactors x nFactors*/, *sum_ov /*nFactors x 1*/, o, *vj /*nFactors x 1*/, *Gwi /*nFactors x 1*/,
           *work, size, *mean, *var, *temp, *rnd;
    int i,j,k,m, thisIndex, otherIndex, outputMeanVar, oIndex, lwork=-1, info, symCheck, role;
    char jobz = 'V', uplo = 'L';

    symCheck = (*debug) - 2;

    if(*option == 1){
        outputMeanVar = 0;
    }else if(*option == 2){
        outputMeanVar = 1;
    }else if(*option == 3){
        outputMeanVar = 1;
    }else error("Unknown option: %d", *option);

    vj     = (double*)Calloc(*nFactors, double);
    Gwi    = (double*)Calloc(*nFactors, double);
    sum_ov = (double*)Calloc(*nFactors, double);
    sum_vv = (double*)Calloc((*nFactors)*(*nFactors), double);
    mean   = (double*)Calloc(*nFactors, double);
    var    = (double*)Calloc((*nFactors)*(*nFactors), double);
    temp   = (double*)Calloc((*nFactors)*(*nFactors), double);

    GetRNGstate();

    for(i=0; i<*nLevelsThisEff; i++){

        thisIndex = i+1;
        for(j=0; j<*nFactors; j++) sum_ov[j] = 0;
        for(j=0; j<(*nFactors)*(*nFactors); j++) sum_vv[j] = 0;
        int nObs_for_i = 0;

        // begin: compute sum_ov and sum_vv
        for(role=0; role<2; role++){
        	int *obsIndex = NULL, *oiNum = NULL, *oiStart=NULL, *thisEffIndex=NULL, *otherEffIndex=NULL;
        	if(role == 0){
        		// user i is an author
        		obsIndex = (int*)to_obsIndex;  oiNum = (int*)to_oiNum;  oiStart = (int*)to_oiStart;
        		thisEffIndex = (int*)toIndex;  otherEffIndex = (int*)fromIndex;
        	}else if(role == 1){
        		// user i is an voter
        		obsIndex = (int*)from_obsIndex; oiNum = (int*)from_oiNum;   oiStart = (int*)from_oiStart;
        		thisEffIndex = (int*)fromIndex; otherEffIndex = (int*)toIndex;
        	}else DIE_HERE;

			for(j=0; j<oiNum[i]; j++){

				nObs_for_i++;

				oIndex = obsIndex[R_VEC(oiStart[i]+j)];
				otherIndex = otherEffIndex[R_VEC(oIndex)];
				if(otherIndex == i+1) error("self votes are not allowed");

				if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
				if(*debug > 0) CHK_R_INDEX(otherIndex, *nLevelsThisEff);
				if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");

				o = obs[R_VEC(oIndex)];

	            if((*nLevelsThirdEff) > 0){
	                int thirdIndex = thirdEffIndex[R_VEC(oIndex)];
	                if(*debug > 0) CHK_R_INDEX(thirdIndex, *nLevelsThirdEff);
	            	for(k=1; k<=*nFactors; k++)
	            		vj[R_VEC(k)] = sample[R_MAT(otherIndex,k,*nLevelsThisEff)] *
	            		               thirdEff[R_MAT(thirdIndex,k,*nLevelsThirdEff)];
	            }else{
	            	for(k=1; k<=*nFactors; k++) vj[R_VEC(k)] = sample[R_MAT(otherIndex,k,*nLevelsThisEff)];
	            }

				double var_y_thisObs = 0;
				if((*nObsVar) == 1)            var_y_thisObs = obsVar[0];
				else if((*nObsVar) == (*nObs)) var_y_thisObs = obsVar[R_VEC(oIndex)];
				else error("nVar_y = %d, nObs = %d", *nObsVar, *nObs);

				for(k=0; k<*nFactors; k++) sum_ov[k] += (o * vj[k]) / var_y_thisObs;
				for(k=0; k<*nFactors; k++)
					for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += (vj[k] * vj[m]) / var_y_thisObs;
			}
        }
        // end: compute sum_ov and sum_vv

        for(k=1; k<=*nFactors; k++) Gwi[R_VEC(k)] = priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)];

        if((*nPriorVar) == 1){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[0]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[0];
        }else if((*nPriorVar) == (*nFactors)){
            for(k=0; k<*nFactors; k++) sum_vv[C_MAT(k,k,*nFactors)] += (1/priorVar[k]);
            for(k=0; k<*nFactors; k++) Gwi[k] /= priorVar[k];
        }else if((*nPriorVar) == (*nLevelsThisEff)*(*nFactors)*(*nFactors)){

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) temp[C_MAT(k,m,*nFactors)] = priorVar[C_3DA(i,k,m,*nLevelsThisEff,*nFactors)];
            sym_inv_byCholesky(temp, nFactors, &symCheck);
            // Now, temp is var_u[i]^-1

            for(k=0; k<*nFactors; k++)
                for(m=0; m<*nFactors; m++) sum_vv[C_MAT(k,m,*nFactors)] += temp[C_MAT(k,m,*nFactors)];

            // reuse the space of vj
            for(k=0; k<*nFactors; k++) vj[k] = Gwi[k];
            for(k=0; k<*nFactors; k++){
                Gwi[k] = 0;
                for(m=0; m<*nFactors; m++) Gwi[k] += temp[C_MAT(k,m,*nFactors)] * vj[m];
            }

        }else error("nVar_eff = %d, nThisEff = %d", *nPriorVar, *nLevelsThisEff);

        // Now, sum_vv = var^-1
        //      Gwi = var_u[i]^-1 G wi

        if((*nFactors) == 1){
            var[0]  = 1/sum_vv[0];
            mean[0] = var[0] * (sum_ov[0] + Gwi[0]);
            sample[i] = rnorm(mean[0], sqrt(var[0]));
            if(outputMeanVar){
                posteriorMean[i] = mean[0];
                posteriorVar[i]  = var[0];
            }
        }else{
            double *eigen_val, *eigen_vec;
            eigen_val = vj;
            eigen_vec = sum_vv;

            if(*debug > 2) CHK_SYMMETRIC(eigen_vec, *nFactors);

            //
            // Compute the variance-covariance matrix
            //
            // Allocate workspace for eigen value decomposition (only allocate once)
            if(lwork == -1){
                F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, &size, &lwork, &info);
                if(info != 0) error("error in dsyev(...)");
                lwork = (int)size;
                work  = (double*)Calloc(lwork, double);
            }

            F77_NAME(dsyev)(&jobz, &uplo, nFactors, eigen_vec, nFactors, eigen_val, work, &lwork, &info);
            if(info != 0) error("error in dsyev(...)");
            for(j=0; j<*nFactors; j++) eigen_val[j] = 1/eigen_val[j];
            // Now, eigen_val, eigen_vec are the eigen values and vectors of var

            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];
            for(j=0; j<*nFactors; j++)
                for(k=0; k<*nFactors; k++){
                    var[C_MAT(j,k,*nFactors)] = 0;
                    for(m=0; m<*nFactors; m++) var[C_MAT(j,k,*nFactors)] += temp[C_MAT(j,m,*nFactors)] * eigen_vec[C_MAT(k,m,*nFactors)];
                }
            // Now, var is the variance-covariance matrix

            //
            // Compute the mean vector
            //
            for(j=0; j<*nFactors; j++) sum_ov[j] += Gwi[j];
            for(j=0; j<*nFactors; j++){
                mean[j] = 0;
                for(k=0; k<*nFactors; k++) mean[j] += var[C_MAT(j,k,*nFactors)] * sum_ov[k];
            }

            if(outputMeanVar){
                for(j=0; j<*nFactors; j++) posteriorMean[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
                for(j=0; j<*nFactors; j++)
                    for(k=0; k<*nFactors; k++) posteriorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)] = var[C_MAT(j,k,*nFactors)];
            }

            if(*debug >= 2 && nObs_for_i == 0){
                for(k=1; k<=*nFactors; k++){
                    CHK_SAME_NUMBER("mean[k] != fittedEff[k]", mean[R_VEC(k)], priorMean[R_MAT(thisIndex,k,*nLevelsThisEff)]);
                }
                if((*nPriorVar) == 1){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[0]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else if((*nPriorVar) == (*nFactors)){
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            if(j==k){ CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[k]);}
                            else{     CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], 0);}
                        }
                }else{
                    for(j=0; j<*nFactors; j++)
                        for(k=0; k<*nFactors; k++){
                            CHK_SAME_NUMBER("var != var_eff", var[C_MAT(j,k,*nFactors)], priorVar[C_3DA(i,j,k,*nLevelsThisEff,*nFactors)]);
                        }
                }
            }
            //
            //  Generate the random vector
            //
			// reuse the space allocated for Gwi
			rnd = Gwi;

			for(j=0; j<*nFactors; j++) eigen_val[j] = sqrt(eigen_val[j]);
			for(j=0; j<*nFactors; j++)
				for(k=0; k<*nFactors; k++) temp[C_MAT(j,k,*nFactors)] = eigen_vec[C_MAT(j,k,*nFactors)] * eigen_val[k];

			for(j=0; j<*nFactors; j++) rnd[j] = norm_rand();

			for(j=0; j<*nFactors; j++){
				sample[C_MAT(i,j,*nLevelsThisEff)] = mean[j];
				for(k=0; k<*nFactors; k++) sample[C_MAT(i,j,*nLevelsThisEff)] += temp[C_MAT(j,k,*nFactors)] * rnd[k];
			}
        }
    }
    PutRNGstate();

    Free(sum_ov);  Free(sum_vv);  Free(vj);    Free(Gwi);
    Free(mean);    Free(var);     Free(temp);
    if(lwork > 0) Free(work);
}

/**
 * contextSample can be set the same as contextPosteriorMean
 * globalSample  can be set the same as globalPosteriorMean
 * If so, the output will be the sample only.
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
){

	const int nObs = (*numObs), nItems = (*numLevelsThisEff), nCategories = (*numContexts),
		      nObsVar = (*numObsVar), nContextPriorVar = (*numContextPriorVar), nGlobalPriorVar = (*numGlobalPriorVar),
		      nCatOffset = (*numContextOffset);
	const int *itemIndex = thisEffIndex, *categoryIndex=context;
	double *b_mean=contextPosteriorMean, *b_var=contextPosteriorVar,
		   *a_mean=globalPosteriorMean,  *a_var=globalPosteriorVar;
	const double *b_offset=contextOffset, *var_obs=obsVar, *var_b=contextPriorVar, *var_a=globalPriorVar;
	// Model:
	// 		obs[n] ~ N(mean=b[i,k],                    var=var_obs[n])
	//      b[i,k] ~ N(mean=b_offset[i,k] + q[k]*a[i], var=var_b[i,k])
	//      a[i]   ~ N(mean=0,                         var=var_a[i])

	if(nObsVar != 1 && nObsVar != nObs) error("nObsVar != 1 && nObsVar != nObs");
	if(nContextPriorVar != nCategories && nContextPriorVar != nItems*nCategories) error("nContextPriorVar != nCategories && nContextPriorVar != nItem*nCategories");
	if(nGlobalPriorVar != 1 && nGlobalPriorVar != nItems) error("nContextPriorVar != nCategories && nContextPriorVar != nItem*nCategories");

	// Compute sufficient statistics
	// 		b_mean[i,k] = sum_j { o_{ijk}/var_obs_{ijk} } = C_ik
	//       b_var[i,k] = sum_j {       1/var_obs_{ijk} } = F_ik
	for(int n=0; n < nItems*nCategories; n++){ b_mean[n] = 0; b_var[n] = 0; }
	for(int j=0; j<nObs; j++){
		int i = itemIndex[j] - 1;     CHK_C_INDEX(i,nItems);
		int k = categoryIndex[j] - 1; CHK_C_INDEX(k,nCategories);
		int ik = C_MAT(i,k,nItems);
		double var = (nObsVar==1 ? var_obs[0] : var_obs[j]);

		double o = obs[j];
		if(nCatOffset == nItems)                  o -= b_offset[i];
		else if(nCatOffset == nItems*nCategories) o -= b_offset[ik];
		else if(nCatOffset != 0) error("nCatOffset = %d", nCatOffset);

		b_mean[ik] += o / var;
		b_var[ik]  += 1 / var;
	}

	GetRNGstate();

	// Compute E[a[i] | obs] and Var[a[i] | obs]
	for(int i=0; i<nItems; i++){
		double sum_for_mean = 0; // sum_k E[a_i | obs_ik]/Var[a_i | obs_ik]
		double sum_for_var  = 0; // sum_k (1/Var[a_i | obs_ik] - 1/var_a)
		for(int k=0; k<nCategories; k++){
			double var_b_this = (nContextPriorVar==nCategories ? var_b[k] : var_b[C_MAT(i,k,nItems)]);
			double F = b_var[C_MAT(i,k,nItems)];
			if(F == 0) continue;
			double C = b_mean[C_MAT(i,k,nItems)];
			if(var_b_this > 0){
				double A = (1/var_b_this) + F;
				// E[a_i | obs_ik]/Var[a_i | obs_ik] = (C * q[k]) / (A * var_b[k])
				// (1/Var[a_i | obs_ik] - 1/var_a)   = (q[k]^2 / var_b[k]) * (1 - 1/(A * var_b[k]))
				sum_for_mean += (C * q[k]) / (A * var_b_this);
				sum_for_var  += (q[k]*q[k] / var_b_this) * (1 - 1/(A * var_b_this));
			}else if(var_b_this == 0){
				STOP("not yet support prior var = 0");
			}else STOP("var_b_this < 0");
		}
		if((*verbose) >= 100) printf("   i=%d:  sum.for.a.mean=%f, sum.for.a.var=%f\n", i+1, sum_for_mean, sum_for_var);

		// Compute E[a[i] | obs] and Var[a[i] | obs]
		double var_a_this = (nGlobalPriorVar==1 ? var_a[0] : var_a[i]);
		a_var[i]  = 1/( (1/var_a_this) + sum_for_var );
		a_mean[i] = a_var[i] * sum_for_mean;

		// draw for the globalSample
		globalSample[i] = rnorm(a_mean[i], sqrt(a_var[i]));
	}

	// Compute E[b[i,k] | a, obs] and Var[b[i,k] | a, obs]
	for(int i=0; i<nItems; i++){ for(int k=0; k<nCategories; k++){
		int ik = C_MAT(i,k,nItems);
		double var_b_this = (nContextPriorVar==nCategories ? var_b[k] : var_b[C_MAT(i,k,nItems)]);
		double A = 1 / ( (1/var_b_this) + b_var[ik] );
		double D = (A * q[k]) / var_b_this;

		b_mean[ik] = D * globalSample[i] + A * b_mean[ik];
		b_var[ik]  = A;

		if(nCatOffset == nItems)                  b_mean[ik] += b_offset[i];
		else if(nCatOffset == nItems*nCategories) b_mean[ik] += b_offset[ik];
		else if(nCatOffset != 0) error("nCatOffset = %d", nCatOffset);

		contextSample[ik] = rnorm(b_mean[ik], sqrt(b_var[ik]));
	}}

    PutRNGstate();
}


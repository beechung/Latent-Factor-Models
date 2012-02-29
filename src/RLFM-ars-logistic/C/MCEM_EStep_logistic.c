/*
        Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
        Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Liang Zhang
*/

/*
To Compile:
R CMD SHLIB util.c MCEM_EStep.c -o MCEM_EStep.so

*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "util.h"
//#include "logistic.h"
#include "arsspline.h"

///// centered logistic ars code starts from here
// FUNCTION: mainEffect_condMeanVar
//  thisEffIndex:  nObs x 1   otherEffIndex: nObs x 1   y_minus_xb_minus_uv: nObs x 1
//  fittedEff: nThisEff x 1   otherEff: nOtherEff x 1
// For alpha, call
//  mainEffect_condMeanVar(1, user, item, y-xb-uv, g0w, beta, var_y, var_alpha, nObs, nUsers, nItems, sample, NULL, NULL, ...);
// For beta, call
//  mainEffect_condMeanVar(1, item, user, y-xb-uv, d0z, alpha, var_y, var_beta, nObs, nItems, nUsers, sample, NULL, NULL, ...);
// Observation index (consider, say, user i) - R indices (starting from 1, NOT 0)
//  obsIndex[ oiStart[i]+0 ], ..., obsIndex[ oiStart[i]+oiNum[i]-1 ] are the indices of user i's observations
//  in y, x, user, item
void mainEffect_condMeanVarSample_arscid(
				       // OUTPUT
				       double* outSample, //int* naccepts,
				       double* ars_XI, /* ars_ninit X ncov X nThisEff*/
				       //INPUT
				       const int* thisEffIndex /*user or item*/, const int* otherEffIndex /*item or user*/,
				       const double* y, const double* offset, /*xb+uv*/
				       const double* fittedEff /*g0w or d0z*/,
				       const double* otherEff /*beta or alpha*/,
				       const double* var_eff /*var_alpha or var_beta*/,
				       const int* nObs, const int* nThisEff, const int* nOtherEff,
				       const int* obsIndex, const int* oiStart, const int* oiNum,
				       const int* ars_ninit, const double* ars_qcent, const double* ars_xl, const double* ars_xu, const double* ars_alpha,
				       // OTHER
    const int* debug

				       ){
  int i,j,k,m, thisIndex, otherIndex, oIndex;
  //*naccepts = 0;
  GetRNGstate();
  int neval = 0;
  for(i=0; i<*nThisEff; i++){
    thisIndex = i+1;
    if(oiNum[i]>0){
      double* y_thisEff = (double*)calloc(oiNum[i], sizeof(double)); // the observations for this user
      double* offset_thisEff = (double*)calloc(oiNum[i], sizeof(double)); // the offsets for this user
      double* X = (double*)calloc(oiNum[i], sizeof(double)); // the design matrix nobs*ncov, ncov=1 for main effects
      for (j=0;j<oiNum[i];j++)
	X[j] = 1;
      for(j=0; j<oiNum[i]; j++){
	oIndex = obsIndex[R_VEC(oiStart[i]+j)];
	otherIndex = otherEffIndex[R_VEC(oIndex)];
	if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
	if(*debug > 0) CHK_R_INDEX(otherIndex, *nOtherEff);
	if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");
	y_thisEff[j] = y[R_VEC(oIndex)];
	offset_thisEff[j] = offset[R_VEC(oIndex)] + otherEff[R_VEC(otherIndex)];
      }
      //double new_thisEff = 0;
      //int accept = 0;
      int ncov = 1;
      //MHlogistic(&outSample[i], X, y_thisEff, offset_thisEff, &fittedEff[i], var_eff, &oiNum[i], &ncov, &new_thisEff, &accept);
      ARSLOGISTICSPLINE(offset_thisEff, &fittedEff[i], var_eff, X, y_thisEff, &ncov, &oiNum[i], &outSample[i], ars_qcent, ars_ninit, ars_ninit, ars_xl, ars_xu, &ars_XI[i*(*ars_ninit)], ars_alpha, &neval);
      //printf("%f\n",outSample[i]);
      //outSample[i] = new_thisEff;
      //if (accept==1) *naccepts = *naccepts + 1;
      free(y_thisEff);
      free(offset_thisEff);
      free(X);
    }
  }
  PutRNGstate();
}
void factor_condMeanVarSample_arscid(
				   // OUTPUT
				   double* outSample, //int* naccepts,
				   double* ars_XI, /* ars_ninit X ncov X nThisEff*/
				   // INPUT
				   const int* thisEffIndex /*user or item*/, const int* otherEffIndex /*item or user*/,
				   const double* y, const double* offset, /*xb+alpha+beta*/
				   const double* fittedEff /*Gw or Dz*/,
				   const double* otherEff /*v or u*/,
				   const double* var_eff /*var_u or var_v*/,
				   const int* nObs, const int* nThisEff, const int* nOtherEff, const int* nFactors,
				   const int* obsIndex, const int* oiStart, const int* oiNum,
				   const int* ars_ninit, const double* ars_qcent, const double* ars_xl, const double* ars_xu, const double* ars_alpha,
				   // OTHER
    const int* debug
				   ){
  int i,j,k,m, thisIndex, otherIndex, oIndex;
  int* neval = (int*)calloc(*nFactors,sizeof(double));
  //*naccepts = 0;
  GetRNGstate();
  for(i=0; i<*nThisEff; i++){
    thisIndex = i+1;
    if(oiNum[i]>0)
      {
	double* y_thisEff = (double*)calloc(oiNum[i], sizeof(double)); // the observations for this user
	double* offset_thisEff = (double*)calloc(oiNum[i], sizeof(double)); // the offsets for this user
	double* vj = (double*)calloc(oiNum[i]*(*nFactors), sizeof(double)); // the design matrix nobs*ncov, ncov=nfactors
	double* thisEff_i = (double*)calloc((*nFactors),sizeof(double));
	//double* new_thisEff_i = (double*)calloc((*nFactors),sizeof(double));
	double* fittedEff_i = (double*)calloc((*nFactors),sizeof(double));
	double* var_this_eff = (double*)calloc((*nFactors),sizeof(double));

	for (j=0;j<*nFactors;j++)
	  {
	    thisEff_i[j] = outSample[C_MAT(i,j,*nThisEff)];
	    fittedEff_i[j] = fittedEff[C_MAT(i,j,*nThisEff)];
	    var_this_eff[j] = var_eff[j];
	  }
	for(j=0; j<oiNum[i]; j++){
	  oIndex = obsIndex[R_VEC(oiStart[i]+j)];
	  otherIndex = otherEffIndex[R_VEC(oIndex)];
	  if(*debug > 0) CHK_R_INDEX(oIndex, *nObs);
	  if(*debug > 0) CHK_R_INDEX(otherIndex, *nOtherEff);
	  if(*debug > 1) if(thisEffIndex[R_VEC(oIndex)] != i+1) error("error in obsIndex, oiStart, oiNum\n");
	  y_thisEff[j] = y[R_VEC(oIndex)];
	  offset_thisEff[j] = offset[R_VEC(oIndex)];
	  for(k=1; k<=*nFactors; k++) {
	    vj[R_MAT((j+1),k,oiNum[i])] = otherEff[R_MAT(otherIndex,k,*nOtherEff)];
	  }
	}
	//double new_thisEff = 0;
	//  int accept = 0;
	//  MHlogistic(thisEff_i, vj, y_thisEff, offset_thisEff, fittedEff_i, var_eff, &oiNum[i], nFactors, new_thisEff_i, &accept);
	//for (j=0;j<*nFactors;j++)
	//  printf("%f ",thisEff_i[j]);
	//printf("\n");
	ARSLOGISTICSPLINE(offset_thisEff, fittedEff_i, var_this_eff, vj, y_thisEff, nFactors, &oiNum[i], thisEff_i, ars_qcent, ars_ninit, ars_ninit, ars_xl, ars_xu, &ars_XI[i*(*ars_ninit)*(*nFactors)], ars_alpha, neval);

	for (j=0;j<*nFactors;j++)
	  {
	    //printf("%f ",thisEff_i[j]);
	    outSample[C_MAT(i,j,*nThisEff)] = thisEff_i[j];
	    //if (accept==1) *naccepts = *naccepts + 1;
	  }
	//printf("\n");
	free(y_thisEff);
	free(offset_thisEff);
	free(vj);
	free(thisEff_i);
	//free(new_thisEff_i);
	free(fittedEff_i);
	free(var_this_eff);
      }
  }
  PutRNGstate();
  free(neval);
}


void computeMeanSumvar_arsc(
			    // OUTPUT
			    double *mean, double *sumvar,
			    // INPUT
			    const double *sum, const double *sos /* sum of squares */, const int length, const int nSamples
			    ){
  int k;
  double ratio = ((double)nSamples) / (nSamples - 1.0);
  sumvar[0] = 0;
  for(k=0; k<length; k++){
    mean[k] = sum[k] / nSamples;
    sumvar[0] += (sos[k] / (nSamples - 1.0)) - (ratio * mean[k] * mean[k]);
  }
}

void computeMeanSumvarFactor_arsc(
			     // OUTPUT
			     double *mean, double *sumvar /* nFactors x 1*/,
			     // INPUT
			     const double *sum, const double *sos /* sum of squares */, const int length, const int nFactors, const int nSamples
			     ){
  int k,l;
  double ratio = ((double)nSamples) / (nSamples - 1.0);
  for (k=0;k<nFactors;k++)
    sumvar[k] = 0;
  for (l=0; l<nFactors;l++)
    for(k=0; k<length; k++){
      mean[C_MAT(k,l,length)] = sum[C_MAT(k,l,length)] / nSamples;
      sumvar[l] += (sos[C_MAT(k,l,length)] / (nSamples - 1.0)) - (ratio * mean[C_MAT(k,l,length)] * mean[C_MAT(k,l,length)]);
    }
}

void InitializeXI_arsc(double* xi, const int ninit, const int ncov, const int nthisEff, const double xu, const double xl)
{
  int i,j,k;
  for (i=0;i<nthisEff;i++)
    for (j=0;j<ncov;j++)
      for (k=0;k<ninit;k++)
        {
          double tmp = xl + (k + 2.0)*(xu - xl)/(ninit + 1.0);
          if (tmp>=xu) tmp = tmp - .1;
          if (tmp<=xl) tmp = tmp + .1;
          xi[i*(ncov*ninit)+j*ninit+k] = tmp;
        }
}
void MCEM_EStep_arscid(
		     // OUTPUT
		     double* o_mean /*nObs x 1*/,
		     double* alpha_mean/*nUsers x 1*/,    double* alpha_sumvar/*1x1*/,
		     double* beta_mean/*nItems x 1*/,     double* beta_sumvar/*1x1*/,
		     double* u_mean/*nUsers x nFactors*/, double* u_sumvar/*nFactors x 1*/,
		     double* v_mean/*nItems x nFactors*/, double* v_sumvar/*nFactors x 1*/,
		     double* ars_XI_alpha, double* ars_XI_beta,
		     double* ars_XI_u, double* ars_XI_v,
		     //double* acceptrate_maineff, // the overall acceptance rate for main effects
		     //double* acceptrate_fact, // the overall acceptance rate for factors
		     // INPUT
		     const int* nSamples, const int* nBurnin,
		     const int* user/*nObs x 1*/,           const int* item/*nObs x 1*/,
		     const double* y/*nObs x 1*/,           const double* xb/*nObs x 1*/,
		     const double* g0w/*nUsers x 1*/,       const double* d0z/*nItems x 1*/,
		     const double* Gw/*nUsers x nFactors*/, const double* Dz/*nItems x nFactors*/,
		     const double* alpha_in/*nUsers x 1*/,    const double* beta_in/*nItems x 1*/,
		     const double* u_in/*nUsers x nFactors*/, const double* v_in/*nItems x nFactors*/,
		     const double* var_alpha, const double* var_beta, 
		     const double* var_u, const double* var_v, /*nFactors x 1*/
		     const int* nObs, const int* nUsers, const int* nItems, const int* nFactors,
		     const int* ars_ninit, const double* ars_qcent, const double* ars_xl, const double* ars_xu, const double* ars_alpha,
		     // OTHER
		     const int* debug, const int* main_effects, const int* beta_int, const int* center

		     ){
  int *obsIndex_user, *oiStart_user, *oiNum_user,
    *obsIndex_item, *oiStart_item, *oiNum_item,
    s, k, f, user_i, item_j, naccepts;
  double *alpha, *beta, *u, *v,
    *alpha_sum, *alpha_sos, *beta_sum, *beta_sos, *u_sum, *u_sos, *v_sum, *v_sos,
    *temp, *xb_plus_uv, *xb_plus_alpha_beta, uv;
  double *ars_xl_v = (double*)calloc(1, sizeof(double));
  ars_xl_v[0] = 0;

  //    double *ars_XI_alpha, *ars_XI_beta, *ars_XI_u, *ars_XI_v;
  //double accept_maineff_denom = 0, accept_maineff_numeri = 0,
  //             accept_fact_denom = 0, accept_fact_numeri = 0;


  // Allocate space for sum and sum-of-squares
  alpha_sum = (double*)calloc(*nUsers, sizeof(double));
  beta_sum  = (double*)calloc(*nItems, sizeof(double));
  u_sum     = (double*)calloc((*nUsers)*(*nFactors), sizeof(double));
  v_sum     = (double*)calloc((*nItems)*(*nFactors), sizeof(double));
  alpha_sos = (double*)calloc(*nUsers, sizeof(double));
  beta_sos  = (double*)calloc(*nItems, sizeof(double));
  u_sos     = (double*)calloc((*nUsers)*(*nFactors), sizeof(double));
  v_sos     = (double*)calloc((*nItems)*(*nFactors), sizeof(double));
  // Allocate space for the observation indices
  obsIndex_user = (int*)calloc(*nObs, sizeof(int));
  oiStart_user  = (int*)calloc(*nUsers, sizeof(int));
  oiNum_user    = (int*)calloc(*nUsers, sizeof(int));
  obsIndex_item = (int*)calloc(*nObs,   sizeof(int));
  oiStart_item  = (int*)calloc(*nItems, sizeof(int));
  oiNum_item    = (int*)calloc(*nItems, sizeof(int));
  // Allocate space for XI ars
  //ars_XI_alpha = (double*)calloc((*ars_ninit)*(*nUsers), sizeof(double));
  //ars_XI_beta = (double*)calloc((*ars_ninit)*(*nItems), sizeof(double));
  //ars_XI_u = (double*)calloc((*ars_ninit)*(*nFactors)*(*nUsers), sizeof(double));
  //ars_XI_v = (double*)calloc((*ars_ninit)*(*nFactors)*(*nItems), sizeof(double));

  // Allocate temp space
  temp = (double*)calloc(*nObs, sizeof(double));

  memset(o_mean,0,(*nObs)*sizeof(double));
  // Use the memory space of the output to store the current alpha, beta, u and v
  alpha = alpha_mean; beta = beta_mean; u = u_mean; v = v_mean;
  // Use the temp space for both xb_plus_uv and xb_plus_alpha_beta
  xb_plus_uv = temp;
  xb_plus_alpha_beta = temp;


  // Create Observation indices for users and items
  generateObsIndex(obsIndex_user, oiStart_user, oiNum_user, user, nObs, nUsers, debug);
  generateObsIndex(obsIndex_item, oiStart_item, oiNum_item, item, nObs, nItems, debug);

  // Initialize alpha, beta, u, v
  for(k=0; k<*nUsers; k++) alpha[k] = alpha_in[k];
  for(k=0; k<*nItems; k++) beta[k]  = beta_in[k];
  for(k=0; k<(*nUsers)*(*nFactors); k++) u[k] = u_in[k];
  for(k=0; k<(*nItems)*(*nFactors); k++) v[k] = v_in[k];

  // Initialize ars_XI_alpha, etc.
  InitializeXI_arsc(ars_XI_alpha, *ars_ninit, 1, *nUsers, *ars_xu, *ars_xl);
  InitializeXI_arsc(ars_XI_beta, *ars_ninit, 1, *nItems, *ars_xu, *ars_xl);
  InitializeXI_arsc(ars_XI_u, *ars_ninit, *nFactors, *nUsers, *ars_xu, *ars_xl);
  InitializeXI_arsc(ars_XI_v, *ars_ninit, *nFactors, *nItems, *ars_xu, *ars_xl_v);

  for(s=0; s<(*nSamples+*nBurnin); s++){

    // Compute xb+uv
    for(k=0; k<*nObs; k++){
      user_i = user[k]; item_j = item[k];
      if(*debug > 0) CHK_R_INDEX(user_i, *nUsers);
      if(*debug > 0) CHK_R_INDEX(item_j, *nItems);
      uv = 0;
      for(f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
      xb_plus_uv[k] = xb[k] + uv;
    }
    // Sample alpha
    //mainEffect_condMeanVarSample_arsc(alpha, &naccepts, user, item, y, xb_plus_uv, g0w, beta, var_alpha, nObs, nUsers, nItems, obsIndex_user, oiStart_user, oiNum_user, debug);
    //accept_maineff_denom += *nUsers; accept_maineff_numeri += naccepts;
    mainEffect_condMeanVarSample_arscid(alpha, ars_XI_alpha, user, item, y, xb_plus_uv, g0w, beta, var_alpha, nObs, nUsers, nItems, obsIndex_user, oiStart_user, oiNum_user, ars_ninit, ars_qcent, ars_xl, ars_xu,ars_alpha,debug);
    //center
    if(*center==1){ center_array(alpha, *nUsers); }
    // Sample beta
    //mainEffect_condMeanVarSample_arsc(beta, &naccepts, item, user, y, xb_plus_uv, d0z, alpha, var_beta, nObs, nItems, nUsers, obsIndex_item, oiStart_item, oiNum_item, debug);
    //accept_maineff_denom += *nItems; accept_maineff_numeri += naccepts;
    mainEffect_condMeanVarSample_arscid(beta, ars_XI_beta, item, user, y, xb_plus_uv, d0z, alpha, var_beta, nObs, nItems, nUsers, obsIndex_item, oiStart_item, oiNum_item, ars_ninit, ars_qcent, ars_xl, ars_xu,ars_alpha,debug);
    //subtract mean from the betas ...
    if(*center==1 && *beta_int==0){  center_array(beta, *nItems); }

    // Compute y - (xb + alpha + beta)
    for(k=0; k<*nObs; k++){
      user_i = user[k]; item_j = item[k];
      if(*debug > 0) CHK_R_INDEX(user_i, *nUsers);
      if(*debug > 0) CHK_R_INDEX(item_j, *nItems);
      xb_plus_alpha_beta[k] = xb[k] + alpha[R_VEC(user_i)] + beta[R_VEC(item_j)];
    }

    if(*main_effects==0){
      // Sample u
      //factor_condMeanVarSample_arsc(u, &naccepts, user, item, y, xb_plus_alpha_beta, Gw, v, var_u, nObs, nUsers, nItems, nFactors, obsIndex_user, oiStart_user, oiNum_user, debug);
      //accept_fact_denom += *nUsers; accept_fact_numeri += naccepts;
      factor_condMeanVarSample_arscid(u, ars_XI_u, user, item, y, xb_plus_alpha_beta, Gw, v, var_u, nObs, nUsers, nItems, nFactors, obsIndex_user, oiStart_user, oiNum_user, ars_ninit, ars_qcent, ars_xl, ars_xu,ars_alpha,debug);
      if(*center==1){ center_array_2d(u, *nUsers, *nFactors, 2); }
      // Sample v
      //factor_condMeanVarSample_arsc(v, &naccepts, item, user, y, xb_plus_alpha_beta, Dz, u, var_v, nObs, nItems, nUsers, nFactors, obsIndex_item, oiStart_item, oiNum_item, debug);
      //accept_fact_denom += *nItems; accept_fact_numeri += naccepts;
      factor_condMeanVarSample_arscid(v, ars_XI_v, item, user, y, xb_plus_alpha_beta, Dz, u, var_v, nObs, nItems, nUsers, nFactors, obsIndex_item, oiStart_item, oiNum_item, ars_ninit, ars_qcent, ars_xl_v, ars_xu,ars_alpha,debug);
      //if(*center==1){ center_array_2d(v, *nItems, *nFactors, 2); }
    }
    // Ignore the first several samples
    if(s >= *nBurnin){
      // update o
      for(k=0; k<*nObs; k++){
	user_i = user[k]; item_j = item[k];
	if(*debug > 0) CHK_R_INDEX(user_i, *nUsers);
	if(*debug > 0) CHK_R_INDEX(item_j, *nItems);
	uv = 0;
	for(f=1; f<=*nFactors; f++) uv += u[R_MAT(user_i,f,*nUsers)] * v[R_MAT(item_j,f,*nItems)];
	double o = alpha[R_VEC(user_i)] + beta[R_VEC(item_j)] + uv;
	o_mean[k] += o/(*nSamples);
      }
      // update alpha
      for(k=0; k<*nUsers; k++){
	alpha_sum[k] += alpha[k];
	alpha_sos[k] += alpha[k]*alpha[k];
      }
      // update beta
      for(k=0; k<*nItems; k++){
	beta_sum[k] += beta[k];
	beta_sos[k] += beta[k]*beta[k];
      }
      // update u
      for(k=0; k<(*nUsers)*(*nFactors); k++){
	u_sum[k] += u[k];
	u_sos[k] += u[k]*u[k];
      }
      // update v
      for(k=0; k<(*nItems)*(*nFactors); k++){
	v_sum[k] += v[k];
	v_sos[k] += v[k]*v[k];
      }
    }
  }

  computeMeanSumvar_arsc(alpha_mean, alpha_sumvar, alpha_sum, alpha_sos, *nUsers, *nSamples);
  computeMeanSumvar_arsc(beta_mean,  beta_sumvar,  beta_sum,  beta_sos,  *nItems, *nSamples);
  computeMeanSumvarFactor_arsc(u_mean, u_sumvar, u_sum, u_sos, *nUsers, *nFactors, *nSamples);
  computeMeanSumvarFactor_arsc(v_mean, v_sumvar, v_sum, v_sos, *nItems, *nFactors, *nSamples);
  //*acceptrate_maineff = accept_maineff_numeri/accept_maineff_denom;
  //*acceptrate_fact = accept_fact_numeri/accept_fact_denom;

  Free(alpha_sos);
  Free(beta_sos);
  Free(u_sos);
  Free(v_sos);
  Free(alpha_sum);
  Free(beta_sum);
  Free(u_sum);
  Free(v_sum);

  Free(obsIndex_user);
  Free(oiStart_user);
  Free(oiNum_user);
  Free(obsIndex_item);
  Free(oiStart_item);
  Free(oiNum_item);
  //Free(ars_XI_alpha);
  //Free(ars_XI_beta);
  //Free(ars_XI_u);
  //Free(ars_XI_v);
  Free(temp);
}


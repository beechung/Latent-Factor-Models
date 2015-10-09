/*
 * Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
 * Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
 * Author: Deepak Agarwal
 */

#include "utilR.h"
#include "arms.h"

// main function
double logexp(double x){
  double ans;
  if(x > 0)ans = x + log1p(exp(-x)); else ans = log1p(exp(x));
  return ans;
}

double lgtsplinefn(double b, void *regressionProblem){
  REGS *problem; // WW specifies the regression problem
  int i;
  double sum, etanew, a, v, earg, p, alpha;
  problem = (REGS*)regressionProblem;
  alpha = problem->alpha;
  sum = 0.0;

  for(i=1;i <= problem->nobs;++i){
    etanew = problem->eta[R_VEC(i)] +
             problem->X[R_MAT(i,problem->id,problem->nobs)] * (b - problem->betacurr[R_VEC(problem->id)]);
    if(etanew > 0.0){
      earg = exp(-2.0 * alpha * etanew);
      p = 2.0*alpha - 1.0 + 2.0*(1.0-alpha)/(1.0+earg);
    } else {
      earg = exp(2.0*(1-alpha)*etanew);
      p = 2.0*alpha*earg/(1.0 + earg);
    }
    if(problem->Y[R_VEC(i)] > 0) sum += log(p); else sum += log(1-p);
  }
  a =(b - problem->beta0[R_VEC(problem->id)]);
  v = problem->varbeta[R_VEC(problem->id)];
  sum -= .5*a*a/v;
  return sum;
}

void print_REGS(REGS* problem){
  Rprintf("  Regression Problem (REGS object):\n");
  Rprintf("    betacurr[%d]=%f, beta0[%d]=%f, varbeta[%d]=%f\n", problem->id, problem->betacurr[R_VEC(problem->id)], problem->id, problem->beta0[R_VEC(problem->id)], problem->id, problem->varbeta[R_VEC(problem->id)]);
  for(int i=1;i <= problem->nobs; ++i){
    Rprintf("    y[%d]=%f, eta[%d]=%f, X[%d,%d]=%f\n", i, problem->Y[R_VEC(i)], i, problem->eta[R_VEC(i)], i, problem->id, problem->X[R_MAT(i,problem->id,problem->nobs)]);
  }
}

// beta_sample is a nSamples X nFactors result matrix
void fitLogistic(const double *offset, double *beta_mean, double *beta_var,
  double *X, double *y, const int *nFactors, const int *nObs, double* beta_sample,
  const double *qcent, const int *ncent, const int *ninit, const double *x_lower,
  const double *x_upper, const int* nSamples) {
	const double alpha = 0.5;
	int nEvalBuckets = (*nSamples) / 100;
	int* neval = (int*)calloc(nEvalBuckets*(*nFactors),sizeof(int));
	int* nevalCurrent = (int*)calloc(*nFactors,sizeof(int));
	double* one_sample = (double*)calloc(*nFactors, sizeof(double));
    double* x_init = (double*)calloc((*nFactors)*(*ninit), sizeof(double));
	for (int i=1;i<=*nFactors;i++) {
		one_sample[R_VEC(i)] = 0;
	}
	for (int i=1;i<=*ninit;i++) {
		for (int j=1;j<=*nFactors;j++) {
			x_init[R_MAT(i,j,*ninit)] = *x_lower + i*((*x_upper)-(*x_lower))/((*ninit)+1.0);
			if (x_init[R_MAT(i,j,*ninit)] >= (*x_upper)) x_init[R_MAT(i,j,*ninit)] -= 0.1;
			if (x_init[R_MAT(i,j,*ninit)] <= (*x_upper)) x_init[R_MAT(i,j,*ninit)] += 0.1;
		}
	}
	for (int i=1; i<=*nSamples; i++) {
		ARSLOGISTICSPLINE(offset, beta_mean, beta_var, X, y, nFactors, nObs, one_sample,
				          qcent, ncent, ninit, x_lower, x_upper, x_init, &alpha, nevalCurrent);
		for (int j=1;j<=*nFactors;j++) {
			beta_sample[R_MAT(i,j,*nSamples)] = one_sample[R_VEC(j)];
			int ind = (i-1)/100;
			neval[R_MAT(ind+1,j,nEvalBuckets)] += nevalCurrent[R_VEC(j)];
		}
	}
	for (int i=1;i<=nEvalBuckets;i++)
		Rprintf("%d ", neval[R_MAT(i,1,nEvalBuckets)]);
	Rprintf("\n");
	Rprintf("Last xinit values:\n");
	for (int i=1;i<=*ninit;i++) {
		for (int j=1;j<=*nFactors;j++) {
			Rprintf("%f ", x_init[R_MAT(i,j,*ninit)]);
		}
		Rprintf("\n");
	}

	free(neval); free(one_sample);free(x_init);free(nevalCurrent);
}
// See the arsspline.h for documentation of this function
void ARSLOGISTICSPLINE(
  const double *offset, double *beta_mean, double *beta_var,
  double *X, double *y, const int *nFactors, const int *nObs,
  double *beta_sample, const double *qcent, const int *ncent,
  const int *ninit, const double *x_lower, const double *x_upper,
  double *x_init, const double *alpha, int *neval
){
  double xsamp[1], convex, lin, XL, XU, *XI, *eta, *xcent;
  int npoint, nsamp, i, j, flag, count, id;
  REGS R;

  if((*ncent) != (*ninit)){ error("ERROR in ARS: ncent (%d) != ninit (%d)",*ncent,*ninit); }

  convex=1.0; npoint=500; nsamp=1; flag=1; count=0;
  XL = *x_lower; XU= *x_upper;
  eta = (double*)Calloc(*nObs,double);
  XI = (double*)Calloc(*ninit,double);
  xcent = (double*)Calloc(*ncent,double);

  // initializing eta
  for(i = 1; i <= *nObs; ++i){
    eta[R_VEC(i)] = offset[R_VEC(i)];
    lin = 0.0;
    for(j = 1; j <= *nFactors; ++j){
      lin += X[R_MAT(i,j,*nObs)] * beta_sample[R_VEC(j)];
    }
    eta[R_VEC(i)] += lin;
  }

  R.eta = eta; R.Y=y; R.X=X; R.betacurr=beta_sample;
  R.beta0 = beta_mean; R.varbeta=beta_var; R.nobs=*nObs; R.alpha=*alpha;

  for(j = 1; j <= *nFactors; ++j){
    R.id = j;
    for(i = 1; i <= *ninit; ++i) XI[R_VEC(i)] = x_init[R_VEC((*ninit)*(j-1)+i)];

    flag = arms(XI, *ninit, &XL, &XU, lgtsplinefn, &R, &convex, npoint, 0,
    		    &beta_sample[R_VEC(j)], xsamp, nsamp, qcent, xcent, *ncent, &neval[R_VEC(j)]);

    if(flag > 0){ error("Error in ARS: error code = %d (returned from arms)",flag); }
    if(ISNAN(xsamp[0])){
      Rprintf("ARS returns NaN (flag=%d):  lower=%f, upper=%f, convex=%f, npoint=%d\n", flag, XL, XU, convex, npoint);
      Rprintf("  beta_sample=%f, neval=%d\n", beta_sample[R_VEC(j)], neval[R_VEC(j)]);
      Rprintf("  XI ="); for(i = 0; i<*ninit; i++) Rprintf(" %f", XI[i]); Rprintf("\n");
      Rprintf("  qcent ="); for(i = 0; i<*ncent; i++) Rprintf(" %f", qcent[i]); Rprintf("\n");
      Rprintf("  xcent ="); for(i = 0; i<*ncent; i++) Rprintf(" %f", xcent[i]); Rprintf("\n");
      print_REGS(&R);
      error("In %s at line %d\n", __FILE__, __LINE__);
    }

    for(i=1; i <= *ncent; ++i) x_init[R_VEC((*ncent)*(j-1)+i)] = xcent[R_VEC(i)];

    // update etas, beta_sample and R.betacurr
    for(i = 1; i <= *nObs; ++i) eta[R_VEC(i)] += X[R_MAT(i,j,*nObs)] * (xsamp[0] - beta_sample[R_VEC(j)]);
    beta_sample[R_VEC(j)] = xsamp[0];
    R.betacurr[R_VEC(j)]  = xsamp[0];
  }
  Free(XI); Free(eta); Free(xcent);
}


double lgtfn(double b, void *W){
  REG *WW;
  int i,npos,obsid;
  double sum,etanew,a,v;
  WW = W;
  //alpha = WW->alpha;
  sum = 0.0;
  npos = WW->X[R_VEC(WW->id)].nidx;
  for(i=1;i <= npos;++i){
    obsid = (WW->X[R_VEC(WW->id)]).obsidx[R_VEC(i)];
    etanew = WW->eta[R_VEC(obsid)] + ((WW->X[R_VEC(WW->id)]).fval[R_VEC(i)])*(b - WW->betacurr[R_VEC(WW->id)]);
    //printf("id=%d\teta=%f\teta=%f\tfval=%f\n",obsid,etanew,WW->eta[R_VEC(obsid)],WW->X[R_VEC(WW->id)].fval[R_VEC(i)]);
    if(WW->Y[R_VEC(obsid)] > 0) sum -= logexp(-etanew); else sum -= logexp(etanew);
  }

  a=(b - WW->beta0[R_VEC(WW->id)]);
  v= WW->varbeta[R_VEC(WW->id)];
  sum -= .5*a*a/v;
  return sum;
}

// function to initialize array of FIDX.

FIDX *SPSTRUCT(int *obsid,double *featval,int *nfobs,int nf){
  // fid:obsid:featval (featureid, obsid, featureval). Specifying design matrix in sparse format. obsid is sorted by feature index from 1:M,
  //nfobs gives the number of entries per feature.
  // nobs (featureid, nobs)
  // nf: number of features
  int cursor,i,k;
  FIDX *S;
  S=(FIDX *)malloc(nf*sizeof(FIDX));
  cursor=0;
  for(i=1;i<= nf;++i){
    S[R_VEC(i)].nidx = nfobs[R_VEC(i)];
    S[R_VEC(i)].obsidx = (int *)Calloc(nfobs[R_VEC(i)],int);
    S[R_VEC(i)].fval = (double *)Calloc(nfobs[R_VEC(i)],double);

    for(k=0;k < nfobs[R_VEC(i)];++k){
      (S[R_VEC(i)].obsidx)[k] = obsid[cursor + k]; (S[R_VEC(i)].fval)[k] = featval[cursor + k];
    }
    cursor += nfobs[R_VEC(i)];
  }
  return S;
}


void debuglgtfn(double *off,int *obsid,double *featval,int *nfobs,// parameters to specify design matrix in sparse format.
		int *Y,int *nFact,int *nObs,// nFact is number of covariates (including intercept if any)
		double *betaout,double *eta,double *beta0,double *varbeta,int *fid,double *x)
{
  int i,j,id;
  REG R;
  FIDX *S;
  S = SPSTRUCT(obsid,featval,nfobs,*nFact);

  // Initializing eta
  for(i=1;i<= *nObs;++i)eta[R_VEC(i)]=off[R_VEC(i)];
  for(j=1;j<= *nFact;++j){
    for(i=1;i<= S[R_VEC(j)].nidx;++i){
      id=S[R_VEC(j)].obsidx[R_VEC(i)];
      eta[R_VEC(id)] += (S[R_VEC(j)].fval[R_VEC(i)])*betaout[R_VEC(j)];
    }
  }
  //for(i=1;i <= *nObs;++i)printf("eta[%d]=%f\n",i,eta[R_VEC(i)]);
  R.eta = eta; R.Y=Y;R.X=S;R.betacurr=betaout;R.beta0=beta0;R.varbeta=varbeta;R.nobs=*nObs;
  R.id = *fid;
  printf("val=%f\n",lgtfn(*x,&R));


  return;
}

void dbSPSTRUCT(int *obsid,double *featval,int *nfobs,int *NF){
  int i;
  FIDX *S;
  S=SPSTRUCT(obsid,featval,nfobs,*NF);
  for(i=0;i< *NF;++i)printf("nidx=%d\n",S[i].nidx);
  //for(i=0;i < S[0].nidx;++i)printf("obsid=%d\tfval=%f\n",S[0].obsidx[i],S[0].fval[i]);
  return;
}

//off: offset
//beta0: priormean
//varbeta: vector of prior variance
//X: vectorized design matrix read columnwise.
//Y: response
// *nFact: number of covariates (including intercept)
// *nObs: number of observations
// betaout: previous beta value
//ninit : number of points in envelope
// xl: lower bound(envelope)
//xu: upper bound(envelope)

// to be fixed: eta is dynamically allocated, this should be done once outside the program and re-used for all calls to this function.
void ARSLOGISTICSPARSE(double *off,double *beta0,double *varbeta,
		       int *obsid,double *featval,int *nfobs,// parameters to specify design matrix in sparse format.
		       int *Y,int *nFact,int *nObs,// nFact is number of covariates (including intercept if any)
		       double *betaout,double *qcent,double *xcent,int *ncent,
		       int *ninit,double *xl,double *xu,double *xi,double *eta,int *neval,int *verboseC){
  // betaout can be used as initial values if required.
  // first coordinate of betaout is the intercept.
  //sample betaout
  // for(i=1;i<=*nFact+1;++i)sample betaout[R_VEC(i)]
  double xsamp[1],convex,lin,XL,XU,*XI;
  int npoint,nsamp,i,j,flag,count,id;
  REG R;
  FIDX *S;

  convex=1.0;npoint=500;nsamp=1;flag=1;count=0;
  XL = *xl; XU= *xu;
  //eta = (double *)Calloc(*nObs,double);
  XI = (double *)Calloc(*ninit,double);
  S = SPSTRUCT(obsid,featval,nfobs,*nFact);

  // Initializing eta
  for(i=1;i<= *nObs;++i)eta[R_VEC(i)]=off[R_VEC(i)];
  for(j=1;j<= *nFact;++j){
    for(i=1;i<= S[R_VEC(j)].nidx;++i){
      id=S[R_VEC(j)].obsidx[R_VEC(i)];
      eta[R_VEC(id)] += (S[R_VEC(j)].fval[R_VEC(i)])*betaout[R_VEC(j)];
    }
  }

  R.eta = eta; R.Y=Y;R.X=S;R.betacurr=betaout;R.beta0=beta0;R.varbeta=varbeta;R.nobs=*nObs;
  for(j=1;j <= *nFact;++j){
    R.id=j;
    for(i=1;i <= *ninit;++i)XI[R_VEC(i)] = xi[R_VEC((*ninit)*(j-1)+i)];

    flag=arms(XI,*ninit,&XL,&XU,lgtfn,&R,&convex,npoint,0,&betaout[R_VEC(j)],xsamp,nsamp,qcent,xcent,*ncent,&neval[R_VEC(j)]);
    if(flag > 0){printf("err=%d\t in ars",flag); exit(1);}

    for(i=1;i <= *ncent;++i)xi[R_VEC((*ncent)*(j-1)+i)]=xcent[R_VEC(i)];

    //update etas
    for(i=1;i<= S[R_VEC(j)].nidx;++i){
      id=S[R_VEC(j)].obsidx[R_VEC(i)];
      eta[R_VEC(id)] += (S[R_VEC(j)].fval[R_VEC(i)])*(xsamp[0] - betaout[R_VEC(j)]);
    }
    //for(i=1;i <= *nObs;++i)eta[R_VEC(i)] += X[R_MAT(i,j,*nObs)]*(xsamp[0]-betaout[R_VEC(j)]);
    betaout[R_VEC(j)]=xsamp[0];
    R.betacurr[R_VEC(j)] = xsamp[0];
    //if(*verboseC > 0){if( j%1000 == 0)printf("number of features processed=%d\n",j);}

  }
  Free(XI);Free(S);
  return;
}

void ARSLOGISTICSPARSEAUX(double *off,double *beta0,double *varbeta,
		       int *obsid,double *featval,int *nfobs,// parameters to specify design matrix in sparse format.
		       int *Y,int *nFact,int *nObs,// nFact is number of covariates (including intercept if any)
		       double *betaout,double *qcent,double *xcent,int *ncent,
		       int *ninit,double *xl,double *xu,double *xi,double *eta,int *neval,int *verboseC,FIDX *S,double *XI){
  // betaout can be used as initial values if required.
  // first coordinate of betaout is the intercept.
  //sample betaout
  // for(i=1;i<=*nFact+1;++i)sample betaout[R_VEC(i)]
  double xsamp[1],convex,lin,XL,XU;
  int npoint,nsamp,i,j,flag,count,id;
  REG R;


  convex=1.0;npoint=500;nsamp=1;flag=1;count=0;
  XL = *xl; XU= *xu;
  //eta = (double *)Calloc(*nObs,double);
  // XI = (double *)Calloc(*ninit,double);
  //S = SPSTRUCT(obsid,featval,nfobs,*nFact);

  // Initializing eta
  for(i=1;i<= *nObs;++i)eta[R_VEC(i)]=off[R_VEC(i)];
  for(j=1;j<= *nFact;++j){
    for(i=1;i<= S[R_VEC(j)].nidx;++i){
      id=S[R_VEC(j)].obsidx[R_VEC(i)];
      eta[R_VEC(id)] += (S[R_VEC(j)].fval[R_VEC(i)])*betaout[R_VEC(j)];
    }
  }

  R.eta = eta; R.Y=Y;R.X=S;R.betacurr=betaout;R.beta0=beta0;R.varbeta=varbeta;R.nobs=*nObs;
  for(j=1;j <= *nFact;++j){
    R.id=j;
    for(i=1;i <= *ninit;++i)XI[R_VEC(i)] = xi[R_VEC((*ninit)*(j-1)+i)];

    flag=arms(XI,*ninit,&XL,&XU,lgtfn,&R,&convex,npoint,0,&betaout[R_VEC(j)],xsamp,nsamp,qcent,xcent,*ncent,&neval[R_VEC(j)]);
    if(flag > 0){printf("err=%d\t in ars",flag); exit(1);}

    for(i=1;i <= *ncent;++i)xi[R_VEC((*ncent)*(j-1)+i)]=xcent[R_VEC(i)];

    //update etas
    for(i=1;i<= S[R_VEC(j)].nidx;++i){
      id=S[R_VEC(j)].obsidx[R_VEC(i)];
      eta[R_VEC(id)] += (S[R_VEC(j)].fval[R_VEC(i)])*(xsamp[0] - betaout[R_VEC(j)]);
    }
    //for(i=1;i <= *nObs;++i)eta[R_VEC(i)] += X[R_MAT(i,j,*nObs)]*(xsamp[0]-betaout[R_VEC(j)]);
    betaout[R_VEC(j)]=xsamp[0];
    R.betacurr[R_VEC(j)] = xsamp[0];
    //if(*verboseC > 0){if( j%1000 == 0)printf("number of features processed=%d\n",j);}

  }
  return;
}

void ARSSAMP(double *off,double *beta0,double *varbeta,int *obsid,double *featval,int *nfobs,
	     int *Y,int *nFact,int *nObs,double *betaout,double *smean,double *ssq,double *qcent,double *xcent,int *ncent,
	     int *ninit,double *xl,double *xu,double *xi,double *eta,int *neval,int *nsamp,int *verboseC){

  int i,k,flag,vninit,vnFact;
  double *xprev,diff,*XI;
  FIDX *S;

  vninit = *ninit; vnFact= *nFact;
  xprev = (double *)Calloc(vninit*vnFact,double);
  S = SPSTRUCT(obsid,featval,nfobs,vnFact);
  XI = (double *)Calloc(vninit,double);
  for(i=0;i < vninit*vnFact;++i)xprev[i] = xi[i];

  for(i=0;i < vnFact;++i){
    smean[i]=0.0;ssq[i]=0.0;
  }

  for(k=0;k < *nsamp;++k){
    ARSLOGISTICSPARSEAUX(off,beta0,varbeta,obsid,featval,nfobs,Y,nFact,nObs,betaout,qcent,xcent,ncent,ninit,xl,xu,xi,eta,neval,verboseC,S,XI);
    if(*verboseC > 0)printf("beta[0]=%f\n",betaout[0]);
    for(i=0;i < vnFact;++i){
      diff = (betaout[i] - smean[i])/(k+1.0);
      smean[i] += diff;
      //smean[i] = (k*smean[i] + betaout[i])/(k+1.0);
      diff = (betaout[i]*betaout[i] - ssq[i])/(k+1.0);
      ssq[i] += diff;
      //ssq[i] = (k*ssq[i] + betaout[i]*betaout[i])/(k+1.0);
    }
    for(i=0;i < vninit*vnFact;++i){
      if(ISNA(xi[i]) || ISNAN(xi[i]))flag=1;
    }
    if(flag){
      for(i=0;i < vninit*vnFact;++i)xi[i] = xprev[i];}
    else {
      for(i=0;i < vninit*vnFact;++i)xprev[i] = xi[i];
    }
  }
  Free(xprev);Free(XI);Free(S);
  return;
}


// matrix factorization via sgd.
// U: vec(u) where u = [u1:u2:..uM];
void SGD(int *user, int *item, double *U,double *V,int *nFact,double *Y,double *w,int *nobs,double *lrate,double *lambda,int *niter){
  double lam,lr,LR,sum,res,rescum;
  int N,nf,iter,i,j,I,J;

  lam = *lambda; LR = 1.0/(*lrate); N = *nobs; nf = *nFact;
  //printf("nf=%d\tV[125]=%f\tU[12]=%f\n",nf,V[125],U[11]);
  // get initial rmse.
  rescum=0.0;
  for(i=1;i<=N;i++){
      I = user[R_VEC(i)]; J = item[R_VEC(i)];
      sum=0.0;
      for(j=1;j<=nf;j++)sum += U[R_MAT(j,I,nf)]*V[R_MAT(j,J,nf)];
      res = Y[R_VEC(i)] - sum; rescum += w[R_VEC(i)]*res*res;
  }
  printf("initial rmse=%f\n",sqrt(rescum/N));

  for(iter=1;iter <= *niter;++iter){
    lr = 1.0/(iter - 1.0 + LR);
    rescum=0.0;
    for(i=1;i<=N;i++){
      I = user[R_VEC(i)]; J = item[R_VEC(i)];
      sum=0.0;
      for(j=1;j<=nf;j++){


	sum += U[R_MAT(j,I,nf)]*V[R_MAT(j,J,nf)];
      }

      res = Y[R_VEC(i)] - sum; rescum += w[R_VEC(i)]*res*res;
      //if(i%100000==0 || i==N)printf("I=%d,J=%d,j=%d,res=%f\tsum=%f\trescum=%f\n",I,J,j,res,sum,rescum);
      for(j=1;j<=nf;j++){
	U[R_MAT(j,I,nf)] -= 2.0*lr*(lam*U[R_MAT(j,I,nf)] - w[R_VEC(i)]*res*V[R_MAT(j,J,nf)]);
	V[R_MAT(j,J,nf)] -= 2.0*lr*(lam*V[R_MAT(j,J,nf)] - w[R_VEC(i)]*res*U[R_MAT(j,I,nf)]);
      }
      //if(i%100000==0)printf("i=%d,V[125]=%f\n",i,V[125]);
    }
    printf("iter=%d\t rmse=%f\n",iter,sqrt(rescum/N));
  }
  return;
}

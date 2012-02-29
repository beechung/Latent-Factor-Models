/*
 * Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
 * Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
 * Author: Deepak Agarwal
 */


#include <R_ext/Applic.h>
#include <R.h>
#include <Rmath.h>

#define R_VEC(i) (i-1)
#define R_MAT(i,j,nrow) (((j-1)*(nrow))+(i-1))

// structure to keep information for feature id fidx
typedef struct{
  int nidx; // number of observations: size of feature index
  //int isdense; //1 means feature is dense. When that happens, obsidx need not be recorded.
  int *obsidx; // array of size nidx containing the obs that occur with feature
  double *fval; // value of features
  //int fidx; // index of feature
}FIDX;

typedef struct{
  int id;
  int nobs;
  double *eta;
  int *Y;
  //double *X;
  FIDX *X;
  double *betacurr;
  double *beta0;
  double *varbeta;
}REG;

typedef struct{
  int id;
  int nobs;
  double *eta;
  int *Y;
  double *X;
  double *betacurr;
  double *beta0;
  double *varbeta;
  double alpha;
}REGS;

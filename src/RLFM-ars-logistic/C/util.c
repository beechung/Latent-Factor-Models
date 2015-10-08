/*
        Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
        Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Liang Zhang
*/


#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include "util.h"

void center_array(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n // array length
		  ){
  // double sum = 0.0;
  double mean = x[0];
  // get sum
  //for(int i=0; i<n; i++){sum += x[i];}
  //mean = sum / (double) n;
  //get mean ... more stable ...
  for(int i=1; i<n; i++){ 
    mean = ( (double) i / (double) (i+1) ) * mean + ( 1.0 / (double) (i+1) ) * x[i];
  }

  //  Rprintf("mean = %f\n", mean);
  // center
  for(int i=0; i<n; i++){x[i] -= mean;}
}

void center_array_2d(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim // 1 for subtracting row mean, 2 for column
){

  double sum, mean;
  int i,j;

  if(dim == 1){
    for(i = 0; i<n; i++){
      sum = 0.0;
      for(j = 0; j<m; j++){sum += x[C_MAT(i, j, n)];}
      mean = sum / (double) m;
      for(j = 0; j<m; j++){x[C_MAT(i, j, n)] -= mean;}
    }
  } else if(dim == 2){
    for(j = 0; j<m; j++){
      //sum = 0.0;
      //for(i = 0; i<n; i++){sum += x[C_MAT(i, j, n)];}
      //mean = sum / (double) n;

      mean = x[C_MAT(0, j, n)];
      for(i=1; i<n; i++){
	mean = ((double)i/(double)(i+1)) * mean + (1.0/(double) (i+1))*x[C_MAT(i,j,n)];
      }
      //      Rprintf("fac mean col %d = %f\n", j, mean);
      for(i = 0; i<n; i++){x[C_MAT(i, j, n)] -= mean;}
  }

  }
}

void center_array_online(
		  // Array to be centered
		  double *x,
		  // INPUT
		  int n, // array length
		  int * oiNum
		  ){
  // double sum = 0.0;
  double mean = 0.0;
  int nav = 0;
  // get sum
  //for(int i=0; i<n; i++){sum += x[i];}
  //mean = sum / (double) n;
  //get mean ... more stable ...
  for(int i=0; i<n; i++){ 
    if(oiNum[i]>0){
      mean = ( (double) nav / (double) (nav+1) ) * mean +	\
	( 1.0 / (double) (nav+1) ) * x[i];
      nav += 1;
    }
  }

  //  Rprintf("mean = %f\n", mean);
  // center
  for(int i=0; i<n; i++){
    if(oiNum[i]>0){
      x[i] -= mean;
    }
  }
}

void center_array_2d_online(
		     // Array to be centered
		     double *x,
		     // INPUT
		     int n, // number of rows
		     int m, // number of columns
		     int dim, // 1 for subtracting row mean, 2 for column
		     int * oiNum
){

  double sum, mean;
  int i,j;
  int nav;

  // dim == 1 is broken!!!
  if(dim == 1){
    for(i = 0; i<n; i++){
      if(oiNum[i]>0){
	sum = 0.0;
	for(j = 0; j<m; j++){sum += x[C_MAT(i, j, n)];}
	mean = sum / (double) m;
	for(j = 0; j<m; j++){x[C_MAT(i, j, n)] -= mean;}
      }
    }
  } else if(dim == 2){
    for(j = 0; j<m; j++){
      //sum = 0.0;
      //for(i = 0; i<n; i++){sum += x[C_MAT(i, j, n)];}
      //mean = sum / (double) n;
      nav = 0;
      mean = 0.0;
      for(i=0; i<n; i++){
	if(oiNum[i]>0){
	  mean = ((double)nav/(double)(nav+1)) * mean +	\
	    (1.0/(double) (nav+1))*x[C_MAT(i,j,n)];
	  nav += 1;
	}
      }
      //      Rprintf("fac mean col %d = %f\n", j, mean);
      for(i = 0; i<n; i++){
	if(oiNum[i]>0){
	  x[C_MAT(i, j, n)] -= mean;
	}
      }
    }

  }
}



// The intput/output are all R indices (start from 1, NOT 0)
void generateObsIndex2(
    //OUTPUT
    int* obsIndex, // E.g., consider user i
    int* start,    //  y[ obsIndex[ start[i]+(0:(num[i]-1)) ] ]
    int* num,      //  are the observations of user i
    //INPUT
    const int* effIndex /* user or item */,
    const int* nObs, const int* nEff,
    //OTHER
    const int* debug
){
    int *ind, i, j, cIndex, eIndex, oIndex;
    ind = (int*)Calloc(*nEff, int);

    for(i=0; i<*nObs; i++){
        eIndex = effIndex[i];
        if(*debug > 0){
            CHK_R_INDEX(eIndex, *nEff);
        }
        num[R_VEC(eIndex)]++;
    }

    start[0] = 1; ind[0] = 1;
    for(i=1; i<*nEff; i++){
        start[i] = start[i-1]+num[i-1];
        ind[i] = start[i];
    }

    for(i=0; i<*nObs; i++){
        cIndex = R_VEC(effIndex[i]);
        obsIndex[R_VEC(ind[cIndex])] = i+1;
        ind[cIndex]++;
    }

    if(*debug > 0){
        for(i=0; i<*nEff; i++){
            if(ind[i] != start[i]+num[i]) error("logical error (level 1)");
            if(*debug > 1){
                for(j=0; j<num[i]; j++){
                    oIndex = obsIndex[R_VEC(start[i]+j)];
                    if(effIndex[R_VEC(oIndex)] != i+1)
                        error("logical error (level 2)");
                }
            }
        }
    }

    Free(ind);
}

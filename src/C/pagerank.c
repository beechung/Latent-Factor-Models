/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "util.h"

// Power iteration
// lambda * ev[i] = sum_j { (w[j] * A[j,i] + (1 - w[j]) * pi[i]) * ev_prev[j] }
//                = sum_j { w[j] * A[j,i] * ev_prev[j] } + pi[i] * sum_j { (1 - w[j]) * ev_prev[j] }
//
// Initial value ev_prev[j] = in[j]
// lambda s.t. sum_i ev[i] = sum_i in[i]
//
void power_iteration(
    double *ev,        // OUTPUT eigenvector: nNodes x 1
    double *lambda,    // OUTPUT eigenvalue: 1x1
    int *nIterUsed,    // OUTPUT 1x1
    const double *in,  // INPUT  vector: nNodes x 1
    const double *pi,  // INPUT  prior vector: nNodes x 1
    const double *w,   // INPUT  teleporting probability vector: nNodes x 1
    const int *A_from, // INPUT transition probability matrix: nEdges x 1 (index start from 0)
    const int *A_to,   // INPUT transition probability matrix: nEdges x 1 (index start from 0)
    const double *A_value, // INPUT transition probability matrix: nEdges x 1
    const int *nNodes, const int *nEdges,
    const int *nIter,  // Number of iterations
    const double *eps, // Stopping criterion (see option)
    const int *option,   // -1: no normalization (run nIter iterations)
                         //  0: run nIter iterations
                         //  1: stop when nNodes * max_i abs(ev[i] - ev_out[i]) / sum(in) < eps (no more than nIter)
    const int *debug // 0: no debugging.  1: check w.  2: positivity.  3: sum-up to one.
){
    double in_sum = 0;
    for(int j=0; j<*nNodes; j++) in_sum += in[j];

    // sanity check
    if(*debug){
        double pi_sum = 0;
        for(int j=0; j<*nNodes; j++){
            if(w[j] > 1 || w[j] < 0) error("w[%d] = %.8g", j, w[j]);
            if((*debug) >= 2){
                if(pi[j] < 0) error("pi[%d] = %.8g", j, pi[j]);
                if(in[j] < 0) error("in[%d] = %.8g", j, in[j]);
            }
            if((*debug) >= 3) pi_sum += pi[j];
        }
        if((*debug) >= 3){
            if(fabs(pi_sum - 1) > 1e-8) error("sum(pi) = %.8g", pi_sum);
            if(fabs(in_sum - 1) > 1e-8) error("sum(in) = %.8g", pi_sum);
        }
        for(int k=0; k<*nEdges; k++){
            if(A_from[k] < 0 || A_from[k] >= (*nNodes)) error("A_from[%d] = %d", k, A_from[k]);
            if(A_to[k]   < 0 || A_to[k]   >= (*nNodes)) error("A_to[%d] = %d",   k, A_to[k]);
            if((*debug) >= 2){
                if(A_value[k] <= 0) error("A_value[%d] = %.8g", k, A_value[k]);
            }
        }
    }

    // initialize ev_prev
    double *ev_prev = NULL;
    if((*nIter) == 1){
        ev_prev = (double*)in;
    }else{
        ev_prev = (double*)Calloc(*nNodes, double);
        for(int i=0; i<*nNodes; i++) ev_prev[i] = in[i];
    }

	(*nIterUsed) = 0;
    // power iteration
    for(int iter=0; iter<(*nIter); iter++){
        if(iter > 0){
            // check for termination
            if((*option) > 0){
                double max_diff = 0;
                for(int i=0; i<*nNodes; i++){
                    double diff = fabs(ev[i] - ev_prev[i]);
                    if(diff > max_diff) max_diff = diff;
                }
                if((*nNodes) * max_diff / in_sum < (*eps)) break;
            }

            // initialize ev_prev
            for(int i=0; i<*nNodes; i++) ev_prev[i] = ev[i];
        }

		(*nIterUsed)++;
        // initialize ev: ev[i] = pi[i] * sum_j { (1 - w[j]) * ev_prev[j] }
        double multiplier = 0;
        for(int j=0; j<*nNodes; j++) multiplier += (1 - w[j]) * ev_prev[j];
        for(int i=0; i<*nNodes; i++) ev[i] = pi[i] * multiplier;

        // go through the edges: ev[i] += sum_j { w[j] * A[j,i] * ev_prev[j] : A[j,i] != 0 }
        for(int k=0; k<*nEdges; k++){
            int j = A_from[k];
            int i = A_to[k];
            double A_ji = A_value[k];
            ev[i] += w[j] * A_ji * ev_prev[j];
        }

        // normalization
        if((*option) != -1){
            double ev_sum = 0;
            for(int i=0; i<*nNodes; i++){
                ev_sum += ev[i];
            }
            double f = in_sum / ev_sum;
            (*lambda) = 1.0;
            if(fabs(f - 1.0) >= 1e-8){
                (*lambda) = 1.0/f;
                for(int i=0; i<*nNodes; i++) ev[i] *= f;
            }
            if((*debug) >= 3 && fabs(ev_sum - 1) > 1e-8) error("sum(ev) = %.8g", ev_sum);
        }

        // sanity check
        if((*debug) >= 2){
            for(int j=0; j<*nNodes; j++) if(ev[j] < 0) error("ev[%d] = %.8g", j, ev[j]);
        }
    }

    if((*nIter) != 1) Free(ev_prev);
}
SEXP power_iteration_Call(
  SEXP ev,        // OUTPUT eigenvector: nNodes x 1
  SEXP lambda,    // OUTPUT eigenvalue: 1x1
  SEXP nIterUsed,    // OUTPUT 1x1
  SEXP in,  // INPUT  vector: nNodes x 1
  SEXP pi,  // INPUT  prior vector: nNodes x 1
  SEXP w,   // INPUT  teleporting probability vector: nNodes x 1
  SEXP A_from, // INPUT transition probability matrix: nEdges x 1 (index start from 0)
  SEXP A_to,   // INPUT transition probability matrix: nEdges x 1 (index start from 0)
  SEXP A_value, // INPUT transition probability matrix: nEdges x 1
  SEXP nNodes, SEXP nEdges,
  SEXP nIter,  // Number of iterations
  SEXP eps, // Stopping criterion (see option)
  SEXP option,   // -1: no normalization (run nIter iterations)
                       //  0: run nIter iterations
                       //  1: stop when nNodes * max_i abs(ev[i] - ev_out[i]) / sum(in) < eps (no more than nIter)
  SEXP debug // 0: no debugging.  1: check w.  2: positivity.  3: sum-up to one.
){
  power_iteration(
    MY_REAL(ev),        // OUTPUT eigenvector: nNodes x 1
    MY_REAL(lambda),    // OUTPUT eigenvalue: 1x1
    MY_INTEGER(nIterUsed),    // OUTPUT 1x1
    MY_REAL(in),  // INPUT  vector: nNodes x 1
    MY_REAL(pi),  // INPUT  prior vector: nNodes x 1
    MY_REAL(w),   // INPUT  teleporting probability vector: nNodes x 1
    MY_INTEGER(A_from), // INPUT transition probability matrix: nEdges x 1 (index start from 0)
    MY_INTEGER(A_to),   // INPUT transition probability matrix: nEdges x 1 (index start from 0)
    MY_REAL(A_value), // INPUT transition probability matrix: nEdges x 1
    MY_INTEGER(nNodes), MY_INTEGER(nEdges),
    MY_INTEGER(nIter),  // Number of iterations
    MY_REAL(eps), // Stopping criterion (see option)
    MY_INTEGER(option),   // -1: no normalization (run nIter iterations)
                         //  0: run nIter iterations
                         //  1: stop when nNodes * max_i abs(ev[i] - ev_out[i]) / sum(in) < eps (no more than nIter)
    MY_INTEGER(debug) // 0: no debugging.  1: check w.  2: positivity.  3: sum-up to one.
  );
  return R_NilValue;
}

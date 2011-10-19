/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "util.h"
#include "multithread.hpp"

struct IntervalBasedMultiThread_Arg {
	IntervalBasedMultiThread* thrObj;
	int id;
	int begin;
	int end;
};

static void* run_IntervalBasedMultiThread_(void* t){
	IntervalBasedMultiThread_Arg *arg = (IntervalBasedMultiThread_Arg *)t;
	arg->thrObj->runThread(arg->id, arg->begin, arg->end);
    pthread_exit(t);
    return NULL;
}

void IntervalBasedMultiThread::run(int nThreads, int minNumPerThread){

	int total = getTotal();
	assert(minNumPerThread >= 1);

	while(total <= nThreads * minNumPerThread && nThreads > 1){
		nThreads--;
	}

	if(nThreads > 0) init(nThreads);

    if(nThreads == 1){

    	runThread(0, 0, total-1);
    	finalize(nThreads);

    }else if(nThreads > 1){

    	IntervalBasedMultiThread_Arg* args = new IntervalBasedMultiThread_Arg[nThreads];
        pthread_t* thread = new pthread_t[nThreads];
        pthread_attr_t attr;
        int rc;
        void *status;

        /* Initialize and set thread detached attribute */
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        for(int t=0; t<nThreads; t++) {

        	args[t].thrObj = this;
        	args[t].id     = t;
        	args[t].begin  = (int)floor( ((t+0.0)/(nThreads+0.0)) * total );
        	args[t].end    = (int)floor( ((t+1.0)/(nThreads+0.0)) * total ) - 1;

        	if(args[t].begin > args[t].end){fprintf(stderr,"ERROR: begin > end (%s at line %d)\n",__FILE__,__LINE__); exit(-1);}

            rc = pthread_create(&thread[t], &attr, run_IntervalBasedMultiThread_, (void*)(&(args[t])));
            if (rc) {
                fprintf(stderr, "ERROR: return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }

        /* Free attribute and wait for the other threads */
        pthread_attr_destroy(&attr);
        for(int t=0; t<nThreads; t++) {
            rc = pthread_join(thread[t], &status);
            if (rc) {
                fprintf(stderr, "ERROR: return code from pthread_join() is %d\n", rc);
                exit(-1);
            }
        }

        finalize(nThreads);

        delete[] args;
        delete[] thread;

    }else{
        fprintf(stderr, "ERROR: nThread is less than 1");
        exit(-1);
    }
}

/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#ifndef MULTITHREAD_HPP_
#define MULTITHREAD_HPP_

/**
 * To implement multi-threading, implement a class extending this one in the following way:
 * (1) Add data members (variables inside the class) for the input and output data
 *        Usually pointers to input & output data
 * (2) Implement a constructor that initializes the data members and a destructor to free the space (if needed)
 * (4) Implement the init(nThreads) method.
 * (4) Implement the runThread(id, begin, end) method.
 * (5) Implement the finalize(nThreads) method.
 * (6) Implement the getTotal() method, which gives the total number of elements to be processed.
 *
 * To use the class, do the following:
 * 	IntervalBasedMultiThread mt(...); // use the constructor
 *  mt.run(n,m); // you do not need to implement run()
 *
 * This is how mt.run(n,m) works intuitively.
 *  Partition [0, getTotal()-1] into n equiwidth intervals
 *  mt.init(n);
 *	for(i=0; i<n; i++){
 *		mt.runThread(i, begin_i, end_i);  // These n threads are running concurrently
 *	}
 *	finalize(n);
 *
 */
class IntervalBasedMultiThread {
public:
	// Please implement the following ----------------------
	virtual void init(int nThreads) = 0;
	virtual void runThread(int id, int begin, int end) = 0;
	virtual void finalize(int nThreads) = 0;
	virtual int  getTotal(void) = 0;
	//------------------------------------------------------

	// You do not need to implement run(...)
	virtual void run(int nThreads, int minNumPerThread);
};

#endif /* MULTITHREAD_HPP_ */

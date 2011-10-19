/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#ifndef TREE_SMOOTHING_H_
#define TREE_SMOOTHING_H_

#include "cube.hpp"

class GaussianTreeSmoother {

private:
	double rootMean;
	double *levelVar; // node at level i ~ N(parentMean, levelVar[i])
	double obsVar;    // observation ~ N(parentMean, obsVar)
	DimSchemaHierarchical *tree;

public:

	/* NOTE:
	 * (1) Each node in the hierarchy is indexed by an ID number in [0, N-1], where
	 *     N is the total number of nodes in the tree.
	 * (2) NodeID=0 must be the root.
	 */
	GaussianTreeSmoother(const int *nodeID, const int *parentID, const int nNodes, const int nEdges){
		tree = new DimSchemaHierarchical(nodeID, parentID, nNodes, nEdges, ROOT_AT_LEVEL_0);
		levelVar = new double[tree->numLevels()];
		rootMean = 0;
	}
	~GaussianTreeSmoother(void){
		delete tree;  delete[] levelVar;
	}
	const DimSchemaHierarchical* getTree(void) const { return tree; }
	const double* getLevelVar(void) const { return levelVar; }
	int numLevels(void) const { return tree->numLevels();}
	double getRootMean(void) const { return rootMean; }
	double getObsVar(void) const { return obsVar; }

	void smooth(
		// INPUT/OUTPUT:
		double mean[], //  input: mean[leaf nodes] are the sample means
		               // output: mean[i] is the posterior mean of node i
		// OUTPUT:
		double var[],  // var[i] is the posterior variance of node i
		double cov_with_parent[], // can be NULL
		// INPUT:
		const int nObs[], // nObs[i] is #observations at node i
		int nNodes,
		double rootPrior, double obsVar, const double levelVar[], int levelVar_length
	);
};

#endif // TREE_SMOOTHING_H_


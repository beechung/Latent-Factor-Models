/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include <float.h>
#include "util.h"
#include "tree-smoothing.h"

// TODO: validate this function
void GaussianTreeSmoother::smooth(
		// INPUT/OUTPUT:
		double mean[], //  input: mean[leaf nodes] are the sample means
		               // output: mean[i] is the posterior mean of node i
		// OUTPUT:
		double var[],  // var[i] is the posterior variance of node i
		double cov_with_parent[], // can be NULL
		// INPUT:
		const int nObs[], // nObs[i] is #observations at node i
		int nNodes,
		double rootPrior, double obsVar, const double level_var[], int level_var_length
){
	if(nNodes != tree->numNodes()) STOP3("nNodes=%d, tree->numNodes()=%d", nNodes, tree->numNodes());

	int nLevels = tree->numLevels();

	this->rootMean = rootMean;
	this->obsVar = obsVar;
	if(level_var_length == 1){
		for(int i=0; i<nLevels; i++) levelVar[i] = level_var[0];
	}else if(level_var_length == nLevels){
		for(int i=0; i<nLevels; i++) levelVar[i] = level_var[i];
	}else STOP_HERE("level_var_length == tree->numLevels()");

	// Shift to have mean of the root = 0
	for(int node=0; node<nNodes; node++){
		if(tree->numChildren(node) == 0){
			if(nObs[node] <= 0) STOP3("nObs[%d] = %d for a leaf", node, nObs[node]);
			mean[node] -= rootMean;
		}else{
			if(nObs[node] > 0) STOP3("nObs[%d] = %d for an internal node", node, nObs[node]);
			mean[node] = DBL_MAX;
		}
	}

	// Filtering (from leaves up)
	for(int level=nLevels-1; level>=0; level--){
		int nNodeThisLevel = tree->numNodesAtLevel(level);
		const int* nodes = tree->getNodeIDsAtLevel(level);
		double sum_var_upward = 0;  for(int i=0; i<=level; i++) sum_var_upward += levelVar[i];

		for(int k=0; k<nNodeThisLevel; k++){
			int node = nodes[k];
			int nChildren = tree->numChildren(node);

			if(nChildren == 0){
				// leaf node
				double temp = (obsVar + nObs[node]*sum_var_upward);
				mean[node] = sum_var_upward * nObs[node] * mean[node] / temp;
				var[node]  = sum_var_upward * obsVar / temp;
			}else{
				// internal node
				const int* children = tree->getChildren(node);
				int j = level+1;
				if(j >= nLevels) STOP_HERE("j >= nLevels");
				double B_j = sum_var_upward / (sum_var_upward + level_var[j]);

				double sum_of_child_adj_mean = 0;
				double sum_of_child_adj_inv_var = 0;

				for(int i=0; i<nChildren; i++){
					int child = children[i];
					if(mean[child] == DBL_MAX) STOP_HERE("Child node has not been processed");
					double mean_rc = B_j * mean[child];
					double var_rc  = B_j * (B_j * var[child] + levelVar[j]);
					sum_of_child_adj_mean += mean_rc / var_rc;
					sum_of_child_adj_inv_var += 1/var_rc - 1/sum_var_upward;
				}
				if(sum_var_upward + sum_of_child_adj_inv_var <= 0) STOP_HERE("var <= 0");
				var[node] = 1 / (sum_var_upward + sum_of_child_adj_inv_var);
				mean[node] = var[node] * sum_of_child_adj_mean;
			}
		}
	}

	if(cov_with_parent != NULL) cov_with_parent[0] = 0;

	// Smoothing (from root down)
	for(int level=1; level<nLevels; level++){
		int nNodeThisLevel = tree->numNodesAtLevel(level);
		const int* nodes = tree->getNodeIDsAtLevel(level);
		double sum_var_upward = 0;  for(int i=0; i<=level; i++) sum_var_upward += levelVar[i];
		double B_j = (sum_var_upward - levelVar[level]) / sum_var_upward;

		for(int k=0; k<nNodeThisLevel; k++){
			int node = nodes[k];
			int parent = tree->getParent(node);
			if(parent < 0) STOP_HERE("parent < 0");

			double mean_pa_r = B_j * mean[node];
			double var_pa_r  = B_j * (B_j * var[node] + levelVar[level]);
			double var_rr_TIMES_B_j_OVER_var_pa_r = var[node] * B_j / var_pa_r;

			mean[node] = mean[node] + var_rr_TIMES_B_j_OVER_var_pa_r * (mean[parent] - mean_pa_r);
			var[node]  = var[node] + SQR(var_rr_TIMES_B_j_OVER_var_pa_r) * (var[parent] - var_pa_r);

			if(cov_with_parent != NULL){
				cov_with_parent[node] = var_rr_TIMES_B_j_OVER_var_pa_r * var[parent];
			}
		}
	}

	// Shift back
	for(int node=0; node<nNodes; node++) mean[node] += rootMean;
}

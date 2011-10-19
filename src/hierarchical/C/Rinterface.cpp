/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#include "model.h"

extern "C" void fit_BetaBinomial_mean_size(
	// OUTPUT
	double *mu, double *gamma, int *code, int *nIter,
	double *muME, double *gammaME,
	// INPUT
	const int *pos,  // length x 1
	const int *size, // length x 1
	const int *length,
	const int *threshold,
	const int *nSamplesForLarge,
	const double *mu0, const double *var_mu, const double *gamma0, const double *var_gamma,
	const double *epsilon, const double *stepSize,
	const int *maxIter1, const int *nLnsrchStep1, const int *maxIter2, const int *nLnsrchStep2,
	const double *minMean, const double *maxMean, const double *minSize, const double *maxSize,
	const int *verbose, const int *debug
){
	BetaBinMEStats *me  = new BetaBinMEStats(pos, size, NULL, *length);
	(*muME)    = me->getMean();
	(*gammaME) = me->getSize();
	BetaBinStats *stats = new BetaBinStats(pos, size, *length, *threshold, *nSamplesForLarge, *mu0, *var_mu, *gamma0, *var_gamma);
	if((*verbose) >= 10) stats->print(stdout);
	double param[2];
	param[0] = (*muME);  param[1] = (*gammaME);
	(*code) = stats->getPostMeanSize(param,nIter,*minMean,*maxMean,*minSize,*maxSize,*epsilon,*stepSize,*maxIter1,*nLnsrchStep1,*maxIter2,*nLnsrchStep2,*verbose,*debug);
	mu[0] = param[0]; gamma[0] = param[1];
	delete stats;
	delete me;
}

extern "C" void fit_BetaBinomialDAG(
	// OUTPUT
	double *param, // #pNodes x #tNodes x 2
	int *status,   // #pNodes x #tNodes
	// INPUT observation data
	const int *data_pnode_id,  // data_len x 1
	const int *data_t_id,      // data_len x 1
	const int *data_nPositive, // data_len x 1
	const int *data_nTrial,    // data_len x 1
	const int *data_len,       // 1x1
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT treatment DAG
	const int *tDAG_child_id,  // num_tDAG_edges x 1
	const int *tDAG_parent_id, // num_tDAG_edges x 1
	const int *num_tNodes,     // 1x1
	const int *num_tDAG_edges, // 1x1
	// INPUT mapping from treatment IDs to nodes
	const int *tMap_t_id,      // num_tMap_edges x 1
	const int *tMap_tnode_id,  // num_tMap_edges x 1
	const int *num_t_ids,      // 1x1
	const int *num_tMap_edges, // 1x1
	// INPUT for moment estimates
	const int *ME_nRefine, // 1x1
	// INPUT for posterior mode estimates
	const int *option, //-1: use naive estimates (based on mean and var of nPositives/nTrials)
					   // 0: skip posterior mode estimation
					   // 1: use the likelihood criterion for selecting the parents & ignore out of bound cases
					   // 2: use (root_priorMean, root_priorSize) for all nodes
					   // 3: use the likelihood criterion for selecting the parents & retain out of bound cases
	const double *root_priorMean,   const double *root_priorSize,
	const double *var_logitMean, 	const double *var_logSize,
	const int *ME_threshold_numObs, const double *ME_threshold_numTrials,
	const int *nTrials_ForSmall, 	const int *nObs_ForLarge,
	const double *epsilon,  // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
						    //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
	const double *stepSize, // IN: max Newton step size
	const int *maxIter1, const int *nLnsrchStep1, const int *maxIter2, const int *nLnsrchStep2,
	const double *minMean, const double *maxMean, const double *minSize, const double *maxSize,
	const int *minNumObs, const int *minTotalTrials, const int *minTotalPositives,
	// INPUT other options
	const int *verbose, const int *debug
){
	BetaBinomialDAG dag(
			param, status, data_pnode_id, data_t_id, data_nPositive, data_nTrial, *data_len,
			pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges,
			tDAG_child_id, tDAG_parent_id, *num_tNodes, *num_tDAG_edges,
			tMap_t_id, tMap_tnode_id, *num_t_ids, *num_tMap_edges,
			*minMean,*maxMean,*minSize,*maxSize,
			*minNumObs, *minTotalTrials, *minTotalPositives,
			*verbose, ((*debug) != 0)
	);

	dag.computeMomentEstimates(*ME_nRefine, (*option)==-1);

	if((*option) <= 0){
		for(int i=0; i<(*num_pNodes)*(*num_tNodes); i++){
			if(status[i] == ME_TEMP) status[i] = ME;
		}
		return;
	}

	if((*option) == 1 || (*option) == 2 || (*option) == 3){
		BetaBinPriorSetter *priorSetter = NULL;

		if((*option) == 1)      priorSetter = new BetaBinLikelihood_PriorSetter(true);
		else if((*option) == 2) priorSetter = new Fixed_PriorSetter();
		else if((*option) == 3) priorSetter = new BetaBinLikelihood_PriorSetter(false);
		else STOP_HERE("");

		dag.computePosteriorMode(
			*root_priorMean, *root_priorSize,
			*var_logitMean, *var_logSize, *ME_threshold_numObs, *ME_threshold_numTrials,
			*nTrials_ForSmall, *nObs_ForLarge, priorSetter, *epsilon, *stepSize,
			*maxIter1, *nLnsrchStep1, *maxIter2, *nLnsrchStep2
		);

		if(priorSetter != NULL) delete priorSetter;
	}else STOP2("Unknow option: %d", *option);
}

extern "C" void loglik_BetaBinomialDAG(
	// OUTPUT
	double *loglik,      // 1x1
	double *node_loglik, // #pNodes x #tNodes (or NULL)
	// INPUT model
	double *param, // #pNodes x #tNodes x 2
	int *status,   // #pNodes x #tNodes
	// INPUT observation data
	const int *data_pnode_id,  // data_len x 1
	const int *data_t_id,      // data_len x 1
	const int *data_nPositive, // data_len x 1
	const int *data_nTrial,    // data_len x 1
	const int *data_len,       // 1x1
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT treatment DAG
	const int *tDAG_child_id,  // num_tDAG_edges x 1
	const int *tDAG_parent_id, // num_tDAG_edges x 1
	const int *num_tNodes,     // 1x1
	const int *num_tDAG_edges, // 1x1
	// INPUT mapping from treatment IDs to nodes
	const int *tMap_t_id,      // num_tMap_edges x 1
	const int *tMap_tnode_id,  // num_tMap_edges x 1
	const int *num_t_ids,      // 1x1
	const int *num_tMap_edges, // 1x1
	// OPTIONS
	const int *option, //   0: leaves only,  1: the whole DAG,  2: frontiers only
	const int *output_nodeLoglik, // 0: No,  1: Yes
	const int *nTrials_ForSmall, 	const int *nObs_ForLarge,
	const double *minMean, const double *maxMean, const double *minSize, const double *maxSize,
	const int *avgloglik, // whether to compute the average log likelihood
	// INPUT other options
	const int *verbose, const int *debug
){
	BetaBinomialDAG dag(
			param, status, NULL, NULL, NULL, NULL, 0,
			pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges,
			tDAG_child_id, tDAG_parent_id, *num_tNodes, *num_tDAG_edges,
			tMap_t_id, tMap_tnode_id, *num_t_ids, *num_tMap_edges,
			*minMean,*maxMean,*minSize,*maxSize,
			0,0,0,
			*verbose, ((*debug) != 0), false
	);

	if((*output_nodeLoglik) == 0) node_loglik = NULL;
	else if((*output_nodeLoglik) != 1) STOP_HERE("");

	(*loglik) = dag.testLogLikelihood(
			node_loglik, data_pnode_id, data_t_id, data_nPositive, data_nTrial, *data_len,
			*nTrials_ForSmall, *nObs_ForLarge, *option, (*avgloglik)!=0
	);
}


extern "C" void confIntCoverage_BetaBinomialDAG(
	// OUTPUT
	double *coverage,      // 1x1
	double *node_nObsInCI, // #pNodes x #tNodes (optional)
	int *node_nObs,        // #pNodes x #tNodes (optional)
	int *nFailed,
	// INPUT model
	double *param,     // #pNodes x #tNodes x 2
	int *status,       // #pNodes x #tNodes
	const double *prob,
	const int *select, // #pNodes x #tNodes (optional)
	// INPUT observation data
	const int *data_pnode_id,  // data_len x 1
	const int *data_t_id,      // data_len x 1
	const int *data_nPositive, // data_len x 1
	const int *data_nTrial,    // data_len x 1
	const int *data_len,       // 1x1
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT treatment DAG
	const int *tDAG_child_id,  // num_tDAG_edges x 1
	const int *tDAG_parent_id, // num_tDAG_edges x 1
	const int *num_tNodes,     // 1x1
	const int *num_tDAG_edges, // 1x1
	// INPUT mapping from treatment IDs to nodes
	const int *tMap_t_id,      // num_tMap_edges x 1
	const int *tMap_tnode_id,  // num_tMap_edges x 1
	const int *num_t_ids,      // 1x1
	const int *num_tMap_edges, // 1x1
	// OPTIONS
	const int *useSelect,        // 0: select all,  1: use select
	const int *output_nodeStats, // 0: No,  1: Yes
	const double *minMean, const double *maxMean, const double *minSize, const double *maxSize,
	// INPUT other options
	const int *verbose, const int *debug
){
	BetaBinomialDAG dag(
			param, status, NULL, NULL, NULL, NULL, 0,
			pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges,
			tDAG_child_id, tDAG_parent_id, *num_tNodes, *num_tDAG_edges,
			tMap_t_id, tMap_tnode_id, *num_t_ids, *num_tMap_edges,
			*minMean,*maxMean,*minSize,*maxSize,
			0,0,0,
			*verbose, ((*debug) != 0), false
	);

	if((*output_nodeStats) == 0){
		node_nObs = NULL;
		node_nObsInCI = NULL;
	}
	else if((*output_nodeStats) != 1) STOP_HERE("");

	if((*useSelect) == 0){
		select = NULL;
	}
	else if((*useSelect) != 1) STOP_HERE("");

	(*coverage) = dag.confIntervalCoverage(
			node_nObsInCI, node_nObs, *prob, select,
			data_pnode_id, data_t_id, data_nPositive, data_nTrial, *data_len,
			nFailed
	);
}


/**
 * Input must be ordered by inData_t_id
 */
extern "C" void cubeAggregate_pDAG(
	// OUTPUT
	int *outData_pnode_id,  // outData_len x 1
	int *outData_t_id,      // outData_len x 1
	int *outData_nPositive, // outData_len x 1
	int *outData_nTrial,    // outData_len x 1
	// INPUT
	const int *outData_len,      // 1x1
	// INPUT observation data
	const int *inData_pnode_id,  // data_len x 1
	const int *inData_t_id,      // data_len x 1
	const int *inData_nPositive, // data_len x 1
	const int *inData_nTrial,    // data_len x 1
	const int *inData_len,       // 1x1
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT other options
	const int *verbose, const int *debug
){
	SimpleDagIndex *pDagIndex = new SimpleDagIndex(pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges);
	int *nodes      = new int[*num_pNodes];
	int *nPositives = new int[*num_pNodes];
	int *nTrials    = new int[*num_pNodes];

	int k = 0;
	int prev_t_id = -INT_MAX;

	for(int m=0; m<(*inData_len); m++){
		int pNodeID = inData_pnode_id[m];
		int t_id    = inData_t_id[m];
		int nPos    = inData_nPositive[m];
		int nTri    = inData_nTrial[m];

		if(t_id > prev_t_id){
			if(prev_t_id != -INT_MAX){
				for(int i=0; i<(*num_pNodes); i++) if(nTrials[i] > 0){
					if(k >= (*outData_len)) STOP_HERE("index out of bound");
					outData_pnode_id[k]  = i;
					outData_t_id[k]      = prev_t_id;
					outData_nPositive[k] = nPositives[i];
					outData_nTrial[k]    = nTrials[i];
					k++;
				}
			}
			for(int i=0; i<(*num_pNodes); i++){
				nPositives[i] = 0;
				nTrials[i] = 0;
			}
		}else if(t_id < prev_t_id) STOP_HERE("Input data is not sorted");

		int num = pDagIndex->getReachableNodes(nodes, pNodeID, *num_pNodes, (*debug)!=0);
		for(int i=0; i<num; i++){
			nPositives[nodes[i]] += nPos;
			nTrials[nodes[i]]    += nTri;
		}

		prev_t_id = t_id;
	}

	if(prev_t_id != -INT_MAX){
		for(int i=0; i<(*num_pNodes); i++) if(nTrials[i] > 0){
			if(k >= (*outData_len)) STOP_HERE("index out of bound");
			outData_pnode_id[k]  = i;
			outData_t_id[k]      = prev_t_id;
			outData_nPositive[k] = nPositives[i];
			outData_nTrial[k]    = nTrials[i];
			k++;
		}
	}

	if(k != (*outData_len)) STOP_HERE("k != (*outData_len)");

	delete pDagIndex;
	delete[] nodes;
	delete[] nPositives;
	delete[] nTrials;
}

/**
 * Input must be ordered by inData_t_id
 */
extern "C" void cubeAggregate_pDAG_get_nObs(
	// OUTPUT
	int *nObs,      // 1x1
	// INPUT observation data
	const int *inData_pnode_id,  // data_len x 1
	const int *inData_t_id,      // data_len x 1
	const int *inData_len,       // 1x1
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT other options
	const int *verbose, const int *debug
){
	SimpleDagIndex *pDagIndex = new SimpleDagIndex(pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges);
	int *nodes    = new int[*num_pNodes];
	bool *reached = new bool[*num_pNodes];

	int k = 0;
	int prev_t_id = -INT_MAX;

	for(int m=0; m<(*inData_len); m++){
		int pNodeID = inData_pnode_id[m];
		int t_id    = inData_t_id[m];

		if(t_id > prev_t_id){
			if(prev_t_id != -INT_MAX){
				for(int i=0; i<(*num_pNodes); i++) if(reached[i]) k++;
			}
			for(int i=0; i<(*num_pNodes); i++) reached[i] = false;
		}else if(t_id < prev_t_id) STOP_HERE("Input data is not sorted");

		int num = pDagIndex->getReachableNodes(nodes, pNodeID, *num_pNodes, (*debug)!=0);
		for(int i=0; i<num; i++){
			reached[nodes[i]] = true;
		}

		prev_t_id = t_id;
	}

	if(prev_t_id != -INT_MAX){
		for(int i=0; i<(*num_pNodes); i++) if(reached[i]) k++;
	}

	(*nObs) = k;

	delete pDagIndex;
	delete[] nodes;
	delete[] reached;
}


// If the input nFrontiers=0, this function only counts the number of frontier nodes
extern "C" void get_frontierIndices(
	// OUTPUT
	int *indices,    // nFrontiers x 2
	int *nFrontiers, // 1x1
	// INPUT
	const int *status,     // #pNodes x #tNodes
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT treatment DAG
	const int *tDAG_child_id,  // num_tDAG_edges x 1
	const int *tDAG_parent_id, // num_tDAG_edges x 1
	const int *num_tNodes,     // 1x1
	const int *num_tDAG_edges, // 1x1
	// INPUT other options
	const int *verbose, const int *debug
){
	SimpleDagIndex *pDagChildIndex = new SimpleDagIndex(pDAG_parent_id, pDAG_child_id, *num_pNodes, *num_pDAG_edges);
	SimpleDagIndex *tDagChildIndex = new SimpleDagIndex(tDAG_parent_id, tDAG_child_id, *num_tNodes, *num_tDAG_edges);

	const int input_nFrontiers = *nFrontiers;
	const int nTNodes = *num_tNodes;
	const int nPNodes = *num_pNodes;
	int k=0;

	for(int pnode=0; pnode<nPNodes; pnode++){
		int nPChildren = pDagChildIndex->numNextNodes(pnode);
		const int *pChildren = pDagChildIndex->getNextNodes(pnode);

		for(int tnode=0; tnode<nTNodes; tnode++){

			if(status[C_MAT(pnode,tnode,nPNodes)] == LOW_SUPPORT) continue;

			int nTChildren = tDagChildIndex->numNextNodes(tnode);
			const int *tChildren = tDagChildIndex->getNextNodes(tnode);

			bool hasLowSuppChild = false;
			for(int i=0; (i<nPChildren) && (!hasLowSuppChild); i++){
				if(status[C_MAT(pChildren[i],tnode,nPNodes)] == LOW_SUPPORT) hasLowSuppChild=true;
			}
			for(int i=0; (i<nTChildren) && (!hasLowSuppChild); i++){
				if(status[C_MAT(pnode,tChildren[i],nPNodes)] == LOW_SUPPORT) hasLowSuppChild=true;
			}

			if((nPChildren == 0 && nTChildren == 0) || hasLowSuppChild){
				if(input_nFrontiers != 0){
					if(k >= input_nFrontiers) STOP_HERE("index out of bound");
					indices[C_MAT(k,0,input_nFrontiers)] = pnode;
					indices[C_MAT(k,1,input_nFrontiers)] = tnode;
				}
				k++;
			}
		}
	}
	(*nFrontiers) = k;
	if(input_nFrontiers != 0 && input_nFrontiers != k) STOP_HERE("");

	delete pDagChildIndex;
	delete tDagChildIndex;
}

// TODO: Validate this
extern "C" void fit_BetaBinomialCube(
	// OUTPUT
	double *paramMean, // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	double *paramSize, // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	int *status,       // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	// INPUT observation data
	const int *data_pnode_id,  // data_len x 1
	const int *data_tnode_id,  // data_len x nTDim
	const int *data_nPositive, // data_len x 1
	const int *data_nTrial,    // data_len x 1
	const int *data_len,       // 1x1
	const int *nTDim,
	// INPUT population DAG
	const int *pDAG_child_id,  // num_pDAG_edges x 1
	const int *pDAG_parent_id, // num_pDAG_edges x 1
	const int *num_pNodes,     // 1x1
	const int *num_pDAG_edges, // 1x1
	// INPUT treatment Cube
	const int *tcube_dim,       // size: tcube_length
	const int *tcube_node_id,   // size: tcube_length
	const int *tcube_parent_id, // size: tcube_length
	const int *tcube_length,
	// INPUT for posterior mode estimates
	const int *option, // 0: skip posterior mode estimation
					   // 1: use the likelihood criterion for selecting the parents
	const double *root_priorMean,   const double *root_priorSize,
	const double *var_logitMean, 	const double *var_logSize,
	const int *ME_threshold_numObs, const double *ME_threshold_numTrials,
	const int *nTrials_ForSmall, 	const int *nObs_ForLarge,
	const double *epsilon,  // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
						    //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
	const double *stepSize, // IN: max Newton step size
	const int *maxIter1, const int *nLnsrchStep1, const int *maxIter2, const int *nLnsrchStep2,
	const double *minMean, const double *maxMean, const double *minSize, const double *maxSize,
	// INPUT other options
	const int *verbose, const int *debug
){
	DimHierarchyArray tdim(tcube_dim, tcube_node_id, tcube_parent_id, *nTDim, *tcube_length);

	BetaBinomialCube cube(
		paramMean, paramSize, status, data_pnode_id, data_tnode_id, data_nPositive, data_nTrial,
		*data_len, *nTDim,
		pDAG_child_id, pDAG_parent_id, *num_pNodes, *num_pDAG_edges,
		tdim.hierarchy,
		*minMean,*maxMean,*minSize,*maxSize,
		*verbose, ((*debug) != 0)
	);

	cube.computeMomentEstimates();

	if((*option) == 0) return;

	if((*option) == 1){
		BetaBinLikelihood_PriorSetter priorSetter(true);
		cube.computePosteriorMode(
			*root_priorMean, *root_priorSize,
			*var_logitMean, *var_logSize, *ME_threshold_numObs, *ME_threshold_numTrials,
			*nTrials_ForSmall, *nObs_ForLarge, &priorSetter, *epsilon, *stepSize,
			*maxIter1, *nLnsrchStep1, *maxIter2, *nLnsrchStep2
		);
	}else STOP2("Unknow option: %d", *option);
}

extern "C" void get_numDagEdges_from_cube(
	// OUTPUT
	int *nDagNodes, int *nDagEdges,
	// INPUT
	const int *cube_dim,       // size: cube_length
	const int *cube_node_id,   // size: cube_length
	const int *cube_parent_id, // size: cube_length
	const int *cube_length,
	const int *nDim
){
	DimHierarchyArray dim(cube_dim, cube_node_id, cube_parent_id, *nDim, *cube_length);
	(*nDagNodes) = 1;
	for(int i=0; i<*nDim; i++) (*nDagNodes) *= dim.hierarchy[i]->numNodes();
	int* nLevels = new int[*nDim];
	for(int i=0; i<*nDim; i++) nLevels[i] = dim.hierarchy[i]->numLevels();
	LatticeEnumeration lat(nLevels, *nDim);
	// note: leaf level (0,...,0)
	//       root level (nLevels[0]-1, ..., nLevels[nDim-1]-1)
	(*nDagEdges) = 0;
	for(int i=0; i<lat.length; i++){
		int nNodes = 1;
		int nParents = 0;
		for(int d=0; d<*nDim; d++){
			if(lat.level[i][d] < nLevels[d]-1) nParents++;
			nNodes *= dim.hierarchy[d]->numNodesAtLevel(lat.level[i][d]);
		}
		(*nDagEdges) += nNodes * nParents;
	}
	delete[] nLevels;
}

extern "C" void gen_DAG_from_cube(
	// OUTPUT
	int *nodeID, // nDagNodes x 1
	int *dimVec, // nDagNodes x nDim
	int *dag_child,  // nDagEdges x 1  (dag_child, dag_parent) specifies an edge from a
	int *dag_parent, // nDagEdges x 1  child node to its parent
	// INPUT
	const int *cube_dim,       // size: cube_length
	const int *cube_node_id,   // size: cube_length
	const int *cube_parent_id, // size: cube_length
	const int *cube_length,
	const int *nDim, const int *nDagNodes, const int *nDagEdges
){
	DimHierarchyArray dim(cube_dim, cube_node_id, cube_parent_id, *nDim, *cube_length);
	BottomUpCubeProcessor cube(dim.hierarchy, *nDim);
	if((*nDagNodes) != cube.NumSubsets) STOP_HERE("(*nDagNodes) != cube.NumSubsets");
	int k=0;
	int *dimvec = new int[*nDim];
	for(int index=0; index < cube.NumSubsets; index++){
		cube.getDimVIDs(index, dimvec);
		nodeID[index] = index;
		for(int d=0; d<*nDim; d++) dimVec[C_MAT(index,d,cube.NumSubsets)] = dimvec[d];
		for(int d=0; d<*nDim; d++){
			int pa = dim.hierarchy[d]->getParent(dimvec[d]);
			if(pa >= 0){
				int temp = dimvec[d];
				dimvec[d] = pa;
				int parentIndex = cube.getSubsetIndex(dimvec);
				if(k >= (*nDagEdges)) STOP_HERE("k >= (*nDagEdges)");

				dag_child[k] = index;
				dag_parent[k] = parentIndex;
				k++;

				dimvec[d] = temp;
			}
		}
	}
	if(k != (*nDagEdges)) STOP_HERE("k != (*nDagEdges)");
	delete[] dimvec;
}

extern "C" void fit_BetaBinomial_mean_size_old(
	// OUTPUT
	double *mu, double *gamma, int *code, int *nIter,
	// INPUT
	const int *pos,  // length x 1
	const int *size, // length x 1
	const int *length,
	const int *threshold,
	const int *nSamplesForLarge,
	const double *mu0, const double *var_mu, const double *gamma0, const double *var_gamma,
	const double *epsilon, const double *stepSize, const int *maxIter, const int *linesearch,
	const int *transform, const int *verbose, const int *debug
){
	BetaBinStats *stats = new BetaBinStats(pos, size, *length, *threshold, *nSamplesForLarge, *mu0, *var_mu, *gamma0, *var_gamma);
	if((*verbose) >= 10) stats->print(stdout);
	double param[2];
	if((*transform)!=0){
		(*code) = postMode_BetaBinomial_mean_size(param, nIter, stats, *epsilon, *stepSize, *maxIter, *linesearch, *verbose, *debug);
	}else{
		(*code) = postMode_BetaBinomial_mean_size_old(param, nIter, stats, *epsilon, *stepSize, *maxIter, (*linesearch)!=0, *verbose, *debug);
	}
	mu[0] = param[0]; gamma[0] = param[1];
	delete stats;
}

extern "C" void do_inc(double *x){
	(*x)++;
}

extern "C" void assign_to(double *x, double *y){
	(*y) = (*x);
}

/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include "model.h"
#include <stdlib.h>
#include <float.h>

#include "util.h"

//#ifdef USE_R
#include <R.h>
//#else
//#define MATHLIB_STANDALONE
//#endif
#include <Rmath.h>
#include <R_ext/Applic.h>

BetaBinomialCube::BetaBinomialCube(
	double *paramMean, // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	double *paramSize, // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	int *status,       // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	const int *data_pnode_id, // data_len x 1
	const int *data_tnode_id, // data_len x nTDim
	const int *data_nPositive,
	const int *data_nTrial,
	int data_len, int nTDim,
	const int *pDAG_child_id, const int *pDAG_parent_id, int num_pNodes, int num_pDAG_edges,
	DimSchemaHierarchical *tCubeDim[], // nTDim hierarchies
	double minMean, double maxMean,
	double minSize, double maxSize,
	int verbose, bool debug
){
	MyTimer timer;
	if(verbose >= 1) printf("START BetaBinomialCube constructor\n");

	this->verbose = verbose;
	this->debug   = debug;

	pDagIndex = new SimpleDagIndex(pDAG_child_id, pDAG_parent_id, num_pNodes, num_pDAG_edges);
	cubeProcessor = new BottomUpCubeProcessor(tCubeDim, nTDim);

	data_length = data_len;
	data_pNodeID = const_cast<int*>(data_pnode_id);
	data_tNodeID = new R_2DArray<int>(const_cast<int*>(data_tnode_id), data_len, nTDim);
	data_nPositives = const_cast<int*>(data_nPositive);
	data_nTrials = const_cast<int*>(data_nTrial);

	nPNodes     = num_pNodes;
	this->nTDim = nTDim;

	nTrialsForSmall = 0;
	nObsForLarge = 100000;

	root_priorMean = -1;
	root_priorSize = -1;
	var_logitMean = -1;
	var_logSize   = -1;
	ME_threshold_nObs = -1;
	ME_threshold_nTrials = -1;

	this->minMean = minMean;
	this->maxMean = maxMean;
	this->minSize = minSize;
	this->maxSize = maxSize;

	this->paramMean = new R_MultiDimArray<double>(paramMean, cubeProcessor->dimSize, nTDim);
	this->paramSize = new R_MultiDimArray<double>(paramSize, cubeProcessor->dimSize, nTDim);
	this->status = new R_MultiDimArray<int>(status, cubeProcessor->dimSize, nTDim);

	for(int i=0; i<this->paramMean->length; i++) paramMean[i] = -1;
	for(int i=0; i<this->paramSize->length; i++) paramSize[i] = -1;
	for(int i=0; i<this->status->length; i++) status[i] = INVALID;

	// check ID ordering
	for(int i=0; i<nPNodes; i++){
		int nParents = pDagIndex->numNextNodes(i);
		if(nParents > 0){
			const int *parent = pDagIndex->getNextNodes(i);
			for(int j=0; j<nParents; j++) if(parent[j] >= i || parent[j] < 0) STOP3("parent[%d] = %d", i, parent[j]);
		}
	}
	for(int d=0; d<nTDim; d++){
		DimSchemaHierarchical *hier = tCubeDim[d];
		int num = hier->numNodes();
		for(int i=1; i<num; i++) if(hier->getParent(i) >= i) STOP3("parent[%d] = %d", i, hier->getParent(i));
	}

	data_startIndex = new int[nPNodes+1];

	int prev_pNodeID = -1;
	for(int k=0; k<data_length; k++){
		int pNodeID = data_pNodeID[k];
		if(pNodeID > prev_pNodeID){
			if(pNodeID < 0 || pNodeID >= nPNodes) STOP_HERE("pNodeID < 0 || pNodeID >= nPNodes");
			for(int id=prev_pNodeID+1; id<=pNodeID; id++) data_startIndex[id] = k;
		}else if(pNodeID < prev_pNodeID) STOP3("Input observation data is not ordered: %d vs. %d", pNodeID, prev_pNodeID);
		prev_pNodeID = pNodeID;
	}
	for(int id=prev_pNodeID+1; id<=nPNodes; id++) data_startIndex[id] = data_length;

	temp = new int[nTDim+1];
	if(verbose >= 1) printf("END   BetaBinomialCube constructor (used %d sec)\n", timer.wallTime());
	MEStats = NULL;
	PMStats = NULL;
}

BetaBinomialCube::~BetaBinomialCube(void){
	delete pDagIndex;  delete cubeProcessor; delete data_tNodeID;
	delete paramMean;  delete paramSize;     delete status;
	delete[] data_startIndex;
	delete[] temp;
}

void BetaBinomialCube::computeMomentEstimates(void)
{
	MyTimer timer;
	if(verbose >= 1) printf("START computeMomentEstimates()\n");
	int nTNodes = cubeProcessor->NumSubsets;
	BetaBinMEStats **stats_space = new BetaBinMEStats*[nTNodes];
	MEStats = new R_MultiDimArray<BetaBinMEStats*>(stats_space, cubeProcessor->dimSize, cubeProcessor->nDim);
	int *dim = new int[cubeProcessor->nDim];

	action = ME;

	for(pNodeID=0; pNodeID < nPNodes; pNodeID++){

		MyTimer timer2;
		if(verbose >= 2) printf("  pNodeID=%d -------------------------------------------\n", pNodeID);
		for(int i=0; i<nTNodes; i++) stats_space[i] = NULL;

		cubeProcessor->compute(this);

		int nSuccess = 0; int nSkip = 0; int nFail = 0;
		for(int tnodeIndex=0; tnodeIndex<nTNodes; tnodeIndex++){
			cubeProcessor->getDimVIDs(tnodeIndex, dim);
			BetaBinMEStats *stats = MEStats->get(dim);
			if(stats != NULL && stats->num > 2 && stats->sum_p > 0){
				double mu = stats->getMean();
				double gamma = stats->getSize();
				if(mu > 0 && mu < 1 && gamma > 0){
					setMean(pNodeID, dim, mu);
					setSize(pNodeID, dim, gamma);
					setStatus(pNodeID, dim, ME_TEMP);
					nSuccess++;
				}else nFail++;
			}else nSkip++;

			if(verbose >= 5) printf("    (pnode:%d, tnode:%d) => (mean:%.10g, size:%.10g)\n", pNodeID, tnodeIndex, getMean(pNodeID,dim), getSize(pNodeID,dim));
		}

		for(int i=0; i<nTNodes; i++) if(stats_space[i] != NULL) delete stats_space[i];
		if(verbose >=2) printf("  #Succeeded: %d   #Skipped: %d   #Failed: %d\n", nSuccess, nSkip, nFail);
		if(verbose >=2) printf("  Time used: %d sec\n", timer2.wallTime());
	}

	delete MEStats;  MEStats = NULL;
	delete[] stats_space;
	delete[] dim;

	if(verbose >= 1) printf("END   computeMomentEstimates (used %d sec)\n", timer.wallTime());
}

void BetaBinomialCube::computePosteriorMode(
	double root_priorMean, double root_priorSize,
	double varLogitMean, double varLogSize, int ME_threshold_numObs, double ME_threshold_numTrials,
	int nTrials_ForSmall, int nObs_ForLarge, BetaBinPriorSetter *priorSetter,
	double epsilon,    // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
					   //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
	double stepSize,   // IN: max Newton step size
	int maxIter1, int nLnsrchStep1, int maxIter2, int nLnsrchStep2
){
	MyTimer timer;
	if(verbose >= 1) printf("START computePosteriorMode(var_logitMean=%.8g, var_logSize=%.8g)\n", varLogitMean, varLogSize);

	this->root_priorMean = root_priorMean;  this->root_priorSize = root_priorSize;
	var_logitMean = varLogitMean; var_logSize = varLogSize;
	ME_threshold_nObs = ME_threshold_numObs; ME_threshold_nTrials = ME_threshold_numTrials;
	nTrialsForSmall = nTrials_ForSmall;  nObsForLarge = nObs_ForLarge;

	action = PM;

	int nTNodes = cubeProcessor->NumSubsets;
	BetaBinStats **stats_space  = new BetaBinStats*[nTNodes];
	BetaBinSizeInfo *info_space = new BetaBinSizeInfo[nTNodes];
	PMStats = new R_MultiDimArray<BetaBinStats*>(stats_space, cubeProcessor->dimSize, cubeProcessor->nDim);
	sizeInfo = new R_MultiDimArray<BetaBinSizeInfo>(info_space, cubeProcessor->dimSize, cubeProcessor->nDim);
	int *dim = new int[cubeProcessor->nDim];

	for(pNodeID=0; pNodeID < nPNodes; pNodeID++){

		MyTimer timer2;
		if(verbose >= 2) printf("  pNodeID=%d -------------------------------------------\n", pNodeID);
		for(int i=0; i<nTNodes; i++){
			info_space[i].nLargeObs = 0;  info_space[i].nObs = 0;  info_space[i].obsCountNeg_length = 0;
			info_space[i].obsCountPos_length = 0; info_space[i].obsCount_length = 0; info_space[i].totalTrials = 0;
			stats_space[i] = NULL;
		}

		// Compute the statistics
		cubeProcessor->compute(this);

		// Find the posterior modes
		int nFail = 0; int nSuccess = 0; int nSkip = 0;
		double ans[2];
		for(int tnodeIndex=0; tnodeIndex<nTNodes; tnodeIndex++){

			cubeProcessor->getDimVIDs(tnodeIndex, dim);
			BetaBinStats *stats = PMStats->get(dim);

			if(verbose >= 6) printf("    processing (pnode:%d, tnode:%d):\n", pNodeID, tnodeIndex);

			if(stats != NULL){

				stats->finalize();
				int code = priorSetter->set(this, pNodeID, dim, stats, verbose);
				if(code != 0){
					if(verbose >= 10){
						printf("FAILED!! code=%d\n", code);
						stats->print(stdout);
					}
					setStatus(pNodeID, dim, PM_FAILED); nFail++;  continue;
				}

				if(verbose >= 6) printf("      mu0=%.10g,  gamma0=%.10g\n", stats->mu0, stats->gamma0);

				int prevStatus = getStatus(pNodeID, dim);
				if(prevStatus == ME_TEMP || prevStatus == ME){
					ans[0] = getMean(pNodeID,dim);  ans[1] = getSize(pNodeID,dim);
				}else{
					ans[0] = 0;  ans[1] = 0;
				}

				int nIterUsed;
				code = stats->getPostMeanSize(ans,&nIterUsed,minMean,maxMean,minSize,maxSize,epsilon,stepSize,maxIter1,nLnsrchStep1,maxIter2,nLnsrchStep2,verbose-10,debug);
				if(code != 0){
					if(verbose >= 10){
						printf("FAILED!! code=%d, niter=%d\n",code, nIterUsed);
						stats->print(stdout);
					}
					setStatus(pNodeID, dim, PM_FAILED); nFail++; continue;
				}

				setMean(pNodeID, dim, ans[0]);
				setSize(pNodeID, dim, ans[1]);
				setStatus(pNodeID, dim, PM);
				nSuccess ++;

			}else nSkip++;

			if(verbose >= 5) printf("    (pnode:%d, tnode:%d) => (mean:%.10g, size:%.10g)\n", pNodeID, tnodeIndex, getMean(pNodeID,dim), getSize(pNodeID,dim));
		}

		for(int i=0; i<nTNodes; i++) if(stats_space[i] != NULL) delete stats_space[i];

		if(verbose >=2) printf("  #Succeeded: %d   #Failed: %d   #Skipped: %d\n", nSuccess, nFail, nSkip);
		if(verbose >=2) printf("  Time used: %d sec\n", timer2.wallTime());
	}

	delete PMStats;  PMStats  = NULL;
	delete sizeInfo; sizeInfo = NULL;
	delete[] stats_space;
	delete[] info_space;
	delete[] dim;

	if(verbose >= 1) printf("END   computePosteriorMode (used %d sec)\n", timer.wallTime());
}

void BetaBinomialCube::initialize_ME(BottomUpCubeProcessor *processor){
	int *tdim = new int[nTDim];
	for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
		if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
		data_tNodeID->getRow(k, tdim);
		int nPositives = data_nPositives[k];
		int nTrials    = data_nTrials[k];
		BetaBinMEStats *stats = MEStats->get(tdim);
		if(stats == NULL){
			stats = new BetaBinMEStats();  MEStats->set(tdim, stats);
		}
		stats->add(nPositives, nTrials, 1);
	}
	delete[] tdim;
}
void BetaBinomialCube::processOneBaseSubset_ME(int *tnode, int ndim, BottomUpCubeProcessor *processor){ }
void BetaBinomialCube::processOneNonbaseSubset_ME(int *tnode, int ndim, BottomUpCubeProcessor *processor){
	BetaBinMEStats *stats = MEStats->get(tnode);
	if(stats != NULL) STOP_HERE("stats != NULL");
	stats = new BetaBinMEStats();  MEStats->set(tnode, stats);
	int d = processor->getDimWithFewestChildren(tnode);
	int nChildren    = processor->CubeDimension[d]->numChildren(tnode[d]);
	const int *child = processor->CubeDimension[d]->getChildren(tnode[d]);
	int *child_node = new int[nTDim]; for(int i=0; i<nTDim; i++) child_node[i] = tnode[i];
	for(int i=0; i<nChildren; i++){
		child_node[d] = child[i];
		BetaBinMEStats *child_stats = MEStats->get(child_node);
		if(child_stats != NULL) stats->add(child_stats);
	}
	delete[] child_node;
}
void BetaBinomialCube::finalizeCube_ME(BottomUpCubeProcessor *processor){ }

void BetaBinomialCube::initialize_PM(BottomUpCubeProcessor *processor){
	// Compute some basic stats
	int *tdim = new int[nTDim];
	for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
		if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
		data_tNodeID->getRow(k, tdim);
		int nPositives = data_nPositives[k];
		int nTrials    = data_nTrials[k];
		if(nTrials < 0 || nPositives < 0 || nPositives > nTrials) STOP_HERE("input data error");
		BetaBinSizeInfo info = sizeInfo->get(tdim);
		info.nObs++;
		info.totalTrials += nTrials;
		if(nTrials <= nTrialsForSmall){
			if(nTrials > info.obsCount_length)       info.obsCount_length    = nTrials;
			if(nPositives > info.obsCountPos_length) info.obsCountPos_length = nPositives;
			if((nTrials-nPositives) > info.obsCountNeg_length) info.obsCountNeg_length = nTrials-nPositives;
		}else{
			info.nLargeObs++;
		}
		sizeInfo->set(tdim, info);
	}

	// Decide which nodes need to compute posterior modes
	for(int t=0; t<PMStats->length; t++){
		BetaBinSizeInfo info = sizeInfo->getBy1DIndex(t);
		BetaBinStats *stats   = PMStats->getBy1DIndex(t);
		if(stats != NULL) STOP_HERE("stats != NULL");
		if(info.nObs >= ME_threshold_nObs && ((double)info.totalTrials)/info.nObs >= ME_threshold_nTrials){
			// do nothing
		}else{
			stats = new BetaBinStats(nTrialsForSmall,nObsForLarge,info.nLargeObs,info.obsCount_length,info.obsCountPos_length,info.obsCountNeg_length,var_logitMean,var_logSize);
			PMStats->setBy1DIndex(t, stats);
		}
	}

	// Update the stats for computing the posterior mode
	for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
		if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
		data_tNodeID->getRow(k, tdim);
		int nPositives = data_nPositives[k];
		int nTrials    = data_nTrials[k];
		BetaBinStats *stats = PMStats->get(tdim);
		if(stats != NULL){
			stats->add(nPositives, nTrials);
		}else{
			if(getStatus(pNodeID, tdim) == ME_TEMP) setStatus(pNodeID, tdim, ME);
		}
	}
	delete[] tdim;
}
void BetaBinomialCube::processOneBaseSubset_PM(int *tnode, int ndim, BottomUpCubeProcessor *processor){ }
void BetaBinomialCube::processOneNonbaseSubset_PM(int *tnode, int ndim, BottomUpCubeProcessor *processor){

	BetaBinStats *stats = PMStats->get(tnode);
	if(stats != NULL) STOP_HERE("stats != NULL");
	BetaBinSizeInfo info = sizeInfo->get(tnode);
	if(info.nObs != 0) STOP_HERE("info.nObs != 0");

	int d = processor->getDimWithFewestChildren(tnode);
	int nChildren    = processor->CubeDimension[d]->numChildren(tnode[d]);
	const int *child = processor->CubeDimension[d]->getChildren(tnode[d]);
	int *child_node = new int[nTDim]; for(int i=0; i<nTDim; i++) child_node[i] = tnode[i];

	for(int i=0; i<nChildren; i++){
		child_node[d] = child[i];
		BetaBinSizeInfo child_info = sizeInfo->get(child_node);
		info.nLargeObs += child_info.nLargeObs;
		info.nObs += child_info.nObs;
		info.obsCountNeg_length = MAX(info.obsCountNeg_length, child_info.obsCountNeg_length);
		info.obsCountPos_length = MAX(info.obsCountPos_length, child_info.obsCountPos_length);
		info.obsCount_length = MAX(info.obsCount_length, child_info.obsCount_length);
		info.totalTrials += child_info.totalTrials;
	}
	sizeInfo->set(tnode, info);

	if(info.nObs >= ME_threshold_nObs && ((double)info.totalTrials)/info.nObs >= ME_threshold_nTrials){
		if(getStatus(pNodeID, tnode) == ME_TEMP) setStatus(pNodeID, tnode, ME);
		delete[] child_node;
		return;
	}

	BetaBinStats **child_stats = new BetaBinStats*[nChildren];
	for(int i=0; i<nChildren; i++){
		child_node[d] = child[i];
		child_stats[i] = PMStats->get(child_node);
		if(child_stats[i] == NULL){
			delete[] child_node;
			delete[] child_stats;
			return;
		}
	}

	stats = new BetaBinStats(nTrialsForSmall,nObsForLarge,info.nLargeObs,info.obsCount_length,info.obsCountPos_length,info.obsCountNeg_length,var_logitMean,var_logSize);
	stats->set(child_stats, nChildren);
	PMStats->set(tnode, stats);
	delete[] child_node;
	delete[] child_stats;
}
void BetaBinomialCube::finalizeCube_PM(BottomUpCubeProcessor *processor){ }


BetaBinomialDAG::BetaBinomialDAG(
	double *param, // (output)  #pNodes * #tNodes * 2
	int *status,   // (output)  #pNodes * #tNodes
	const int *data_pnode_id, const int *data_t_id, const int *data_nPositive, const int *data_nTrial, int data_len,
	const int *pDAG_child_id, const int *pDAG_parent_id, int num_pNodes, int num_pDAG_edges,
	const int *tDAG_child_id, const int *tDAG_parent_id, int num_tNodes, int num_tDAG_edges,
	const int *tMap_t_id, const int *tMap_tnode_id, int num_t_ids, int num_tMap_edges,
	double minMean, double maxMean, double minSize, double maxSize,
	int minNumObs, int minTotalTrials, int minTotalPositives,
	int verbose, bool debug, bool resetParamAndStatus
){
	this->verbose = verbose;
	this->debug   = debug;

	MyTimer timer;
	if(verbose >= 1) printf("START BetaBinomialDAG constructor\n");

	pDagIndex = new SimpleDagIndex(pDAG_child_id, pDAG_parent_id, num_pNodes, num_pDAG_edges);
	tDagIndex = new SimpleDagIndex(tDAG_child_id, tDAG_parent_id, num_tNodes, num_tDAG_edges);
	tMapping  = new SimpleDagIndex(tMap_t_id, tMap_tnode_id, num_t_ids, num_tMap_edges);

	data_length = data_len;
	data_pNodeID = const_cast<int*>(data_pnode_id);
	data_tID = const_cast<int*>(data_t_id);
	data_nPositives = const_cast<int*>(data_nPositive);
	data_nTrials = const_cast<int*>(data_nTrial);

	nPNodes = num_pNodes;
	nTNodes = num_tNodes;

	nTrialsForSmall = 0;
	nObsForLarge = 100000;

	var_logitMean = -1;
	var_logSize   = -1;

	this->minMean = minMean;
	this->maxMean = maxMean;
	this->minSize = minSize;
	this->maxSize = maxSize;

	this->minNumObs = minNumObs;
	this->minTotalTrials = minTotalTrials;
	this->minTotalPositives = minTotalPositives;

	this->param  = param;
	this->status = status;
	if(resetParamAndStatus){
		for(int i=0; i<nPNodes*nTNodes*2; i++) param[i] = -1;
		for(int i=0; i<nPNodes*nTNodes; i++) status[i] = INVALID;
	}

	// check ID ordering
	for(int i=0; i<nPNodes; i++){
		int nParents = pDagIndex->numNextNodes(i);
		if(nParents > 0){
			const int *parent = pDagIndex->getNextNodes(i);
			for(int j=0; j<nParents; j++) if(parent[j] >= i || parent[j] < 0) STOP3("parent[%d] = %d", i, parent[j]);
		}
	}
	for(int i=0; i<nTNodes; i++){
		int nParents = tDagIndex->numNextNodes(i);
		if(nParents > 0){
			const int *parent = tDagIndex->getNextNodes(i);
			for(int j=0; j<nParents; j++) if(parent[j] >= i || parent[j] < 0) STOP3("parent[%d] = %d", i, parent[j]);
		}
	}

	data_startIndex = new int[nPNodes+1];

	int prev_pNodeID = -1;
	for(int k=0; k<data_length; k++){
		int pNodeID = data_pNodeID[k];
		if(pNodeID > prev_pNodeID){
			if(pNodeID < 0 || pNodeID >= nPNodes) STOP_HERE("pNodeID < 0 || pNodeID >= nPNodes");
			for(int id=prev_pNodeID+1; id<=pNodeID; id++) data_startIndex[id] = k;
		}else if(pNodeID < prev_pNodeID) STOP3("Input observation data is not ordered: %d vs. %d", pNodeID, prev_pNodeID);
		prev_pNodeID = pNodeID;
	}
	for(int id=prev_pNodeID+1; id<=nPNodes; id++) data_startIndex[id] = data_length;

	if(verbose >= 1) printf("END   BetaBinomialDAG constructor (used %d sec)\n", timer.wallTime());
}

BetaBinomialDAG::~BetaBinomialDAG(void){
	delete pDagIndex;
	delete tDagIndex;
	delete tMapping;
	delete[] data_startIndex;
}

void BetaBinomialDAG::computeMomentEstimates(int nRefine, bool naive)
{
	MyTimer timer;
	if(verbose >= 1) printf("START computeMomentEstimates(%d)\n", nRefine);

	if(naive && nRefine != 0) STOP_HERE("naive && nRefine != 0");

	BetaBinMEStats **stats = new BetaBinMEStats*[nTNodes];
	int *tNodes = new int[nTNodes];

	for(int pNodeID=0; pNodeID < nPNodes; pNodeID++){

		MyTimer timer2;
		if(verbose >= 2) printf("  pNodeID=%d -------------------------------------------\n", pNodeID);

		for(int r=0; r<=nRefine; r++){

			for(int i=0; i<nTNodes; i++) stats[i] = new BetaBinMEStats();
			for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
				if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
				int tID 	   = data_tID[k];
				int nPositives = data_nPositives[k];
				int nTrials    = data_nTrials[k];
				int  num   = tMapping->numNextNodes(tID);
				const int* nodes = tMapping->getNextNodes(tID);
				num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);
				for(int j=0; j<num; j++){
					int tNodeID = tNodes[j];
					double w = (r == 0 ? 1 : nTrials/(1.0 + (nTrials-1)/(getSize(pNodeID,tNodeID)+1)));
					stats[tNodeID]->add(nPositives, nTrials, w);
				}
			}

			int nSuccess = 0; int nSkip = 0; int nFail = 0;
			for(int tnode=0; tnode<nTNodes; tnode++){

				if(stats[tnode]->num < minNumObs || stats[tnode]->totalPostives < minTotalPositives || stats[tnode]->totalTrials < minTotalTrials){
					setStatus(pNodeID, tnode, LOW_SUPPORT);
					nSkip++;
				}else if(stats[tnode]->num > 2 && stats[tnode]->sum_p > 0){

					double mu = stats[tnode]->getMean();
					double gamma = stats[tnode]->getSize(naive);

					if(mu > 0 && mu < 1 && gamma > 0){
						mu = FORCE_RANGE(mu, minMean, maxMean);
						gamma = FORCE_RANGE(gamma, minSize, maxSize);
						setMean(pNodeID, tnode, mu);
						setSize(pNodeID, tnode, gamma);
						setStatus(pNodeID, tnode, ME_TEMP);
						nSuccess++;
					}else{
						setStatus(pNodeID, tnode, INVALID);
						nFail++;
					}
				}else{
					setStatus(pNodeID, tnode, INVALID);
					nSkip++;
				}

				if(verbose >= 5) printf("    (pnode:%d, tnode:%d) => (mean:%.10g, size:%.10g)\n", pNodeID, tnode, getMean(pNodeID,tnode), getSize(pNodeID,tnode));
			}

			for(int i=0; i<nTNodes; i++) delete stats[i];
			if(verbose >=2) printf("  Round %d => #Succeeded: %d   #Skipped: %d   #Failed: %d\n", r, nSuccess, nSkip, nFail);
		}

		if(verbose >=2) printf("  Time used: %d sec\n", timer2.wallTime());
	}

	delete[] stats;
	delete[] tNodes;

	if(verbose >= 1) printf("END   computeMomentEstimates (used %d sec)\n", timer.wallTime());
}

void BetaBinomialDAG::computePosteriorMode(
	double root_priorMean, double root_priorSize,
	double varLogitMean, double varLogSize, int ME_threshold_numObs, double ME_threshold_numTrials,
	int nTrials_ForSmall, int nObs_ForLarge, BetaBinPriorSetter *priorSetter,
	double epsilon,    // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
					   //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
	double stepSize,   // IN: max Newton step size
	int maxIter1, int nLnsrchStep1, int maxIter2, int nLnsrchStep2
){
	MyTimer timer;
	if(verbose >= 1) printf("START computePosteriorMode(var_logitMean=%.8g, var_logSize=%.8g)\n", varLogitMean, varLogSize);

	this->root_priorMean = root_priorMean;  this->root_priorSize = root_priorSize;
	var_logitMean = varLogitMean; var_logSize = varLogSize;
	ME_threshold_nObs = ME_threshold_numObs; ME_threshold_nTrials = ME_threshold_numTrials;
	nTrialsForSmall = nTrials_ForSmall;  nObsForLarge = nObs_ForLarge;

	BetaBinStats **stats = new BetaBinStats*[nTNodes];
	int *tNodes = new int[nTNodes];
	int *nLargeObs = new int[nTNodes];
	int *obsCount_length = new int[nTNodes];
	int *obsCountPos_length = new int[nTNodes];
	int *obsCountNeg_length = new int[nTNodes];
	int *nObs = new int[nTNodes];
	int *totalTrials = new int[nTNodes];
	int *totalPostives = new int[nTNodes];
	double temp[2];

	for(int pNodeID=0; pNodeID < nPNodes; pNodeID++){

		MyTimer timer2;
		if(verbose >= 2) printf("  pNodeID=%d -------------------------------------------\n", pNodeID);

		for(int i=0; i<nTNodes; i++) nLargeObs[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCount_length[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCountPos_length[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCountNeg_length[i] = 0;
		for(int i=0; i<nTNodes; i++) nObs[i] = 0;
		for(int i=0; i<nTNodes; i++) totalTrials[i] = 0;
		for(int i=0; i<nTNodes; i++) totalPostives[i] = 0;

		// Compute some basic stats
		for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
			if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
			int tID 	   = data_tID[k];
			int nPositives = data_nPositives[k];
			int nTrials    = data_nTrials[k];
			int  num   = tMapping->numNextNodes(tID);
			const int* nodes = tMapping->getNextNodes(tID);
			num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);
			for(int j=0; j<num; j++){
				int tNodeID = tNodes[j];
				if(nTrials < 0 || nPositives < 0 || nPositives > nTrials) STOP_HERE("input data error");
				nObs[tNodeID]++;
				totalTrials[tNodeID] += nTrials;
				totalPostives[tNodeID] += nPositives;
				if(nTrials <= nTrialsForSmall){
					if(nTrials > obsCount_length[tNodeID])    obsCount_length[tNodeID]    = nTrials;
					if(nPositives  > obsCountPos_length[tNodeID]) obsCountPos_length[tNodeID] = nPositives;
					if((nTrials-nPositives) > obsCountNeg_length[tNodeID]) obsCountNeg_length[tNodeID] = nTrials-nPositives;
				}else{
					nLargeObs[tNodeID]++;
				}
			}
		}

		// Decide which nodes need to compute posterior modes
		for(int t=0; t<nTNodes; t++){
			stats[t] = NULL;
			temp[0] = getMean(pNodeID,t);
			temp[1] = getSize(pNodeID,t);
			if(nObs[t] < minNumObs || totalPostives[t] < minTotalPositives || totalTrials[t] < minTotalTrials){
				if(getStatus(pNodeID,t) != LOW_SUPPORT) STOP_HERE("Inconsistent low support");
			}else if(IS_PARAM_IN_RANGE(temp, minMean, maxMean, minSize, maxSize) &&
					 nObs[t] >= ME_threshold_nObs && ((double)totalTrials[t])/nObs[t] >= ME_threshold_numTrials &&
					 getStatus(pNodeID,t) == ME_TEMP){
				setStatus(pNodeID, t, ME);
			}else{
				stats[t] = new BetaBinStats(nTrialsForSmall,nObsForLarge,nLargeObs[t],obsCount_length[t],obsCountPos_length[t],obsCountNeg_length[t],var_logitMean,var_logSize);
			}
		}

		// Update the stats for computing the posterior mode
		for(int k=data_startIndex[pNodeID]; k<data_startIndex[pNodeID+1]; k++){
			if(data_pNodeID[k] != pNodeID) STOP_HERE("data_pNodeID[k] != pNodeID");
			int tID 	   = data_tID[k];
			int nPositives = data_nPositives[k];
			int nTrials    = data_nTrials[k];
			int  num   = tMapping->numNextNodes(tID);
			const int* nodes = tMapping->getNextNodes(tID);
			num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);
			for(int j=0; j<num; j++){
				int tNodeID = tNodes[j];
				if(stats[tNodeID] != NULL) stats[tNodeID]->add(nPositives, nTrials);
			}
		}

		// Find the posterior modes
		int nFail = 0; int nSuccess = 0; int nSkip = 0; bool setPriorFailed;
		double ans[2];
		for(int tnode=0; tnode<nTNodes; tnode++){

			if(verbose >= 6) printf("    processing (pnode:%d, tnode:%d):\n", pNodeID, tnode);

			if(stats[tnode] != NULL){

				stats[tnode]->finalize();
				setPriorFailed = false;
				int code = priorSetter->set(this, pNodeID, tnode, stats[tnode], verbose);
				if(code != 0){
					if(verbose >= 10){
						printf("FAILED!! code=%d\n", code);
						stats[tnode]->print(stdout);
					}
					setPriorFailed = true;
				}

				if(verbose >= 6) printf("      mu0=%.10g,  gamma0=%.10g\n", stats[tnode]->mu0, stats[tnode]->gamma0);

				if(getStatus(pNodeID, tnode) == ME_TEMP || getStatus(pNodeID, tnode) == ME){
					ans[0] = getMean(pNodeID,tnode);  ans[1] = getSize(pNodeID, tnode);
				}else{
					ans[0] = 0;  ans[1] = 0;
				}

				int nIterUsed;
				code = stats[tnode]->getPostMeanSize(ans,&nIterUsed,minMean,maxMean,minSize,maxSize,epsilon,stepSize,maxIter1,nLnsrchStep1,maxIter2,nLnsrchStep2,verbose-10,debug);
				if(code != 0){
					if(verbose >= 10){
						printf("FAILED!! code=%d, niter=%d\n",code, nIterUsed);
						stats[tnode]->print(stdout);
					}
					setStatus(pNodeID, tnode, PM_FAILED); nFail++; continue;
				}

				ans[0] = FORCE_RANGE(ans[0], minMean, maxMean);
				ans[1] = FORCE_RANGE(ans[1], minSize, maxSize);

				setMean(pNodeID, tnode, ans[0]);
				setSize(pNodeID, tnode, ans[1]);
				if(setPriorFailed){
					setStatus(pNodeID, tnode, PM_FAILED);
					nFail++;
				}else{
					setStatus(pNodeID, tnode, PM);
					nSuccess ++;
				}
			}else nSkip++;

			if(verbose >= 5) printf("    (pnode:%d, tnode:%d) => (mean:%.10g, size:%.10g)\n", pNodeID, tnode, getMean(pNodeID,tnode), getSize(pNodeID,tnode));
		}

		for(int i=0; i<nTNodes; i++) if(stats[i] != NULL) delete stats[i];

		if(verbose >=2) printf("  #Succeeded: %d   #Failed: %d   #Skipped: %d\n", nSuccess, nFail, nSkip);
		if(verbose >=2) printf("  Time used: %d sec\n", timer2.wallTime());
	}

	delete[] stats;
	delete[] tNodes;
	delete[] nLargeObs;
	delete[] obsCount_length;
	delete[] obsCountPos_length;
	delete[] obsCountNeg_length;
	delete[] nObs;
	delete[] totalTrials;
	delete[] totalPostives;

	if(verbose >= 1) printf("END   computePosteriorMode (used %d sec)\n", timer.wallTime());
}

double BetaBinomialDAG::testLogLikelihood(
	double *node_loglik, // (output) #pNodes * #tNodes: LogLikelihood per observation of each node
	                     //          set node_loglik=NULL to disable the output
						 //          0 means no data in the node
	const int *testdata_pNodeID, const int *testdata_tID, const int *testdata_nPositives, const int *testdata_nTrials, int testdata_len,
	int nTrials_ForSmall, int nObs_ForLarge,
	int option, // 0: leaves only,  1: the whole DAG,  2: frontiers only
	bool averageLoglik
){
	MyTimer timer;
	if(verbose >= 1) printf("START testLogLikelihood(option=%d)\n", option);

	BetaBinStats **stats = new BetaBinStats*[nTNodes];
	int *tNodes = new int[nTNodes];
	int *nLargeObs = new int[nTNodes];
	int *obsCount_length = new int[nTNodes];
	int *obsCountPos_length = new int[nTNodes];
	int *obsCountNeg_length = new int[nTNodes];
	int *nObs = new int[nTNodes];
	int *totalTrials = new int[nTNodes];
	bool *computeLoglik = new bool[nTNodes];

	int *testdata_startIndex = new int[nPNodes+1];

	int prev_pNodeID = -1;
	for(int k=0; k<testdata_len; k++){
		int pNodeID = testdata_pNodeID[k];
		if(pNodeID > prev_pNodeID){
			if(pNodeID < 0 || pNodeID >= nPNodes) STOP_HERE("pNodeID < 0 || pNodeID >= nPNodes");
			for(int id=prev_pNodeID+1; id<=pNodeID; id++) testdata_startIndex[id] = k;
		}else if(pNodeID < prev_pNodeID) STOP3("Input observation data is not ordered: %d vs. %d", pNodeID, prev_pNodeID);
		prev_pNodeID = pNodeID;
	}
	for(int id=prev_pNodeID+1; id<=nPNodes; id++) testdata_startIndex[id] = testdata_len;


	int nNodes = 0;
	double loglik = 0;
	int numFailed = 0;

	for(int pNodeID=0; pNodeID < nPNodes; pNodeID++){

		MyTimer timer2;
		if(verbose >= 2) printf("  pNodeID=%d -------------------------------------------\n", pNodeID);

		for(int i=0; i<nTNodes; i++) nLargeObs[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCount_length[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCountPos_length[i] = 0;
		for(int i=0; i<nTNodes; i++) obsCountNeg_length[i] = 0;
		for(int i=0; i<nTNodes; i++) nObs[i] = 0;
		for(int i=0; i<nTNodes; i++) totalTrials[i] = 0;
		for(int i=0; i<nTNodes; i++) computeLoglik[i] = false;

		if(node_loglik != NULL) for(int i=0; i<nTNodes; i++) node_loglik[C_MAT(pNodeID,i,nPNodes)] = 0;

		// Compute some basic stats
		for(int k=testdata_startIndex[pNodeID]; k<testdata_startIndex[pNodeID+1]; k++){
			if(testdata_pNodeID[k] != pNodeID) STOP_HERE("testdata_pNodeID[k] != pNodeID");
			int tID 	   = testdata_tID[k];
			int nPositives = testdata_nPositives[k];
			int nTrials    = testdata_nTrials[k];
			int  num   = tMapping->numNextNodes(tID);
			const int* nodes = tMapping->getNextNodes(tID);

			for(int i=0; i<num; i++) computeLoglik[nodes[i]] = true;

			if(option == 1 || option == 2){
				num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);
			}else if(option == 0){
				for(int i=0; i<num; i++) tNodes[i] = nodes[i];
			}else STOP_HERE("");

			for(int j=0; j<num; j++){
				int tNodeID = tNodes[j];
				if(nTrials < 0 || nPositives < 0 || nPositives > nTrials) STOP_HERE("input data error");
				nObs[tNodeID]++;
				totalTrials[tNodeID] += nTrials;
				if(nTrials <= nTrials_ForSmall){
					if(nTrials > obsCount_length[tNodeID])        obsCount_length[tNodeID]    = nTrials;
					if(nPositives  > obsCountPos_length[tNodeID]) obsCountPos_length[tNodeID] = nPositives;
					if((nTrials-nPositives) > obsCountNeg_length[tNodeID]) obsCountNeg_length[tNodeID] = nTrials-nPositives;
				}else{
					nLargeObs[tNodeID]++;
				}
			}
		}

		if(option == 1){
			for(int t=0; t<nTNodes; t++) computeLoglik[t] = true;
		}else if(option == 2){
			for(int t=0; t<nTNodes; t++){
				if(getStatus(pNodeID,t) == LOW_SUPPORT){
					int num = tDagIndex->numNextNodes(t);
					const int* nodes = tDagIndex->getNextNodes(t);
					for(int i=0; i<num; i++) computeLoglik[nodes[i]] = true;
				}
			}
		}

		// Decide which nodes need to compute loglikelihood
		for(int t=0; t<nTNodes; t++){
			stats[t] = NULL;
			if(nObs[t] > 0 && getStatus(pNodeID,t) > 0 && computeLoglik[t]){
				stats[t] = new BetaBinStats(nTrials_ForSmall,nObs_ForLarge,nLargeObs[t],obsCount_length[t],obsCountPos_length[t],obsCountNeg_length[t],var_logitMean,var_logSize, false);
			}
		}

		// Update the stats for computing loglikelihood
		for(int k=testdata_startIndex[pNodeID]; k<testdata_startIndex[pNodeID+1]; k++){
			if(testdata_pNodeID[k] != pNodeID) STOP_HERE("testdata_pNodeID[k] != pNodeID");
			int tID 	   = testdata_tID[k];
			int nPositives = testdata_nPositives[k];
			int nTrials    = testdata_nTrials[k];
			int  num   = tMapping->numNextNodes(tID);
			const int* nodes = tMapping->getNextNodes(tID);

			if(option == 1){
				num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);
			}else if(option == 0){
				for(int i=0; i<num; i++) tNodes[i] = nodes[i];
			}else STOP_HERE("");

			for(int j=0; j<num; j++){
				int tNodeID = tNodes[j];
				if(stats[tNodeID] != NULL) stats[tNodeID]->add(nPositives, nTrials);
			}
		}

		// Compute loglikelihood
		double ans[2];
		for(int tnode=0; tnode<nTNodes; tnode++){
			int status = getStatus(pNodeID, tnode);
			if(status < 0) numFailed++;
			if(stats[tnode] != NULL && status > 0){
				stats[tnode]->finalize();

				ans[0] = FORCE_RANGE(getMean(pNodeID, tnode), minMean, maxMean);
				ans[1] = FORCE_RANGE(getSize(pNodeID, tnode), minSize, maxSize);

				double v = stats[tnode]->logLikelihood(ans, false, averageLoglik);
				if(node_loglik != NULL) node_loglik[C_MAT(pNodeID,tnode,nPNodes)] = v;
				loglik += v;
				nNodes++;
			}
		}

		for(int i=0; i<nTNodes; i++) if(stats[i] != NULL) delete stats[i];

		if(verbose >=2) printf("  Time used: %d sec\n", timer2.wallTime());
	}

	if(numFailed > 0) printf("  %d nodes do not have validate estimates\n", numFailed);

	delete[] stats;
	delete[] tNodes;
	delete[] nLargeObs;
	delete[] obsCount_length;
	delete[] obsCountPos_length;
	delete[] obsCountNeg_length;
	delete[] nObs;
	delete[] totalTrials;
	delete[] computeLoglik;
	delete[] testdata_startIndex;

	if(verbose >= 1) printf("END   testLogLikelihood (used %d sec)\n", timer.wallTime());

	if(averageLoglik) return loglik/nNodes;
	return loglik;
}

double BetaBinomialDAG::confIntervalCoverage(
	double *nObsInCI,  // (output) #pNodes x #tNodes (NULL means disable)
	int *nObs,         // (output) #pNodes x #tNodes (NULL means disable; 0 means no output for the node)
	double prob, // probability: (prob*100)% confidence interval will be constructed
	const int *select, // #pNodes x #tNodes (NULL means select all; non-zero means selected)
	const int *testdata_pNodeID, const int *testdata_tID, const int *testdata_nPositives, const int *testdata_nTrials, int testdata_len,
	int *nFailed
){
	MyTimer timer;
	if(verbose >= 1) printf("START confIntervalCoverage(prob=%f)\n", prob);

	int *testdata_startIndex = new int[nPNodes+1];
	int *tNodes = new int[nTNodes];
	double pLower = (1-prob)/2;
	double pUpper = 1 - pLower;

	int prev_pNodeID = -1;
	for(int k=0; k<testdata_len; k++){
		int pNodeID = testdata_pNodeID[k];
		if(pNodeID > prev_pNodeID){
			if(pNodeID < 0 || pNodeID >= nPNodes) STOP_HERE("pNodeID < 0 || pNodeID >= nPNodes");
			for(int id=prev_pNodeID+1; id<=pNodeID; id++) testdata_startIndex[id] = k;
		}else if(pNodeID < prev_pNodeID) STOP3("Input observation data is not ordered: %d vs. %d", pNodeID, prev_pNodeID);
		prev_pNodeID = pNodeID;
	}
	for(int id=prev_pNodeID+1; id<=nPNodes; id++) testdata_startIndex[id] = testdata_len;

	double totalObs=0, totalObsInCI=0;
	int numFailed = 0;
	if(nObsInCI != NULL) for(int i=0; i<nPNodes*nTNodes; i++) nObsInCI[i] = 0;
	if(nObs != NULL) for(int i=0; i<nPNodes*nTNodes; i++) nObs[i] = 0;

	bool *failed = new bool[nTNodes];

	for(int pNodeID=0; pNodeID < nPNodes; pNodeID++){
		for(int i=0; i<nTNodes; i++) failed[i] = false;
		for(int k=testdata_startIndex[pNodeID]; k<testdata_startIndex[pNodeID+1]; k++){
			if(testdata_pNodeID[k] != pNodeID) STOP_HERE("testdata_pNodeID[k] != pNodeID");
			int tID 	   = testdata_tID[k];
			int nPositives = testdata_nPositives[k];
			int nTrials    = testdata_nTrials[k];
			int  num   = tMapping->numNextNodes(tID);
			const int* nodes = tMapping->getNextNodes(tID);

			num = tDagIndex->getReachableNodes(tNodes, nodes, num, nTNodes, debug);

			for(int j=0; j<num; j++){
				int tNodeID = tNodes[j];
				if(select != NULL && select[C_MAT(pNodeID,tNodeID,nPNodes)] == 0) continue;
				if(getStatus(pNodeID,tNodeID) <= 0){
					if(getStatus(pNodeID,tNodeID) < 0) failed[tNodeID] = true;
					continue;
				}
				double mean = FORCE_RANGE(getMean(pNodeID, tNodeID), minMean, maxMean);
				double size = FORCE_RANGE(getSize(pNodeID, tNodeID), minSize, maxSize);
				double prev_p = -1;
				double p = pbetabinom_MeanSize(nPositives, nTrials, mean, size, &prev_p);
				if(p < -1e-8 || p > 1+1e-8 || prev_p < -1e-8 || prev_p > 1+1e-8){
					printf("WARNING: p=%.16g and prev_p=%.16g! pNodeID=%d, tNodeID=%d, status=%d, mean=%.16g, size=%.16g, nPos=%d, nTrial=%d\n",
							p, prev_p, pNodeID, tNodeID, getStatus(pNodeID,tNodeID), mean, size, nPositives, nTrials);
					failed[tNodeID] = true;
					continue;
				}

				totalObs++;
				if(nObs != NULL) nObs[C_MAT(pNodeID,tNodeID,nPNodes)]++;

				if(prev_p >= pLower && p <= pUpper){
					totalObsInCI++;
					if(nObsInCI != NULL) nObsInCI[C_MAT(pNodeID,tNodeID,nPNodes)]++;
				}else if(p > pLower && prev_p < pLower){
					double frac = (p-pLower)/(p-prev_p);
					totalObsInCI += frac;
					if(nObsInCI != NULL) nObsInCI[C_MAT(pNodeID,tNodeID,nPNodes)] += frac;
				}else if(p > pUpper && prev_p < pUpper){
					double frac = (pUpper-prev_p)/(p-prev_p);
					totalObsInCI += frac;
					if(nObsInCI != NULL) nObsInCI[C_MAT(pNodeID,tNodeID,nPNodes)] += frac;
				}
			}
		}
		for(int i=0; i<nTNodes; i++) if(failed[i]) numFailed++;
	}
	(*nFailed) = numFailed;

	delete[] testdata_startIndex;
	delete[] tNodes;
	delete[] failed;

	if(verbose >= 1) printf("END   confIntervalCoverage (used %d sec)\n", timer.wallTime());
	return totalObsInCI/totalObs;
}

int Fixed_PriorSetter::set(BetaBinomialDAG *dag, int pnode, int tnode, BetaBinStats *stats, int verbose){
	stats->mu0 = dag->root_priorMean;  stats->gamma0 = dag->root_priorSize;
	return 0;
}
int Fixed_PriorSetter::set(BetaBinomialCube *cube, int pnode, int tvector[], BetaBinStats *stats, int verbose){
	stats->mu0 = cube->root_priorMean;  stats->gamma0 = cube->root_priorSize;
	return 0;
}

int BetaBinLikelihood_PriorSetter::set(BetaBinomialDAG *dag, int pnode, int tnode, BetaBinStats *stats, int verbose){

	stats->mu0 = dag->root_priorMean;  stats->gamma0 = dag->root_priorSize;
	if(pnode == 0 && tnode == 0) return 0;

	double best_logLik = -DBL_MAX;
	double ans[2];

	// rollup along the treatment
	int num = dag->tDagIndex->numNextNodes(tnode);
	if(num == 0 && tnode != 0) STOP2("tNode %d has no parents", tnode);
	const int *tParents = dag->tDagIndex->getNextNodes(tnode);
	for(int i=0; i<num; i++){

		if(dag->getStatus(pnode, tParents[i]) <= 0) continue;

		ans[0] = dag->getMean(pnode, tParents[i]);
		ans[1] = dag->getSize(pnode, tParents[i]);

		if(ignoreOutOfBound && !IS_PARAM_IN_RANGE(ans,dag->minMean,dag->maxMean,dag->minSize,dag->maxSize)) continue;

		double logLik = logLik_BetaBinomial_mean_size(stats, ans, false);

		if(logLik > best_logLik){
			stats->mu0 = ans[0];  stats->gamma0 = ans[1];
			best_logLik = logLik;
		}

		if(verbose >= 10) printf("      (pnode:%d, tnode:%d, tparent:%d) => status=%d, loglik=%f (mu=%f, gamma=%f)\n", pnode, tnode, tParents[i], dag->getStatus(pnode, tParents[i]), logLik, ans[0], ans[1]);
	}

	// rollup along the population
	num = dag->pDagIndex->numNextNodes(pnode);
	if(num == 0 && pnode != 0) STOP2("pNode %d has no parents", pnode);
	const int *pParents = dag->pDagIndex->getNextNodes(pnode);
	for(int i=0; i<num; i++){

		if(dag->getStatus(pParents[i], tnode) <= 0) continue;

		ans[0] = dag->getMean(pParents[i], tnode);
		ans[1] = dag->getSize(pParents[i], tnode);

		if(ignoreOutOfBound && !IS_PARAM_IN_RANGE(ans,dag->minMean,dag->maxMean,dag->minSize,dag->maxSize)) continue;

		double logLik = logLik_BetaBinomial_mean_size(stats, ans, false);

		if(logLik > best_logLik){
			stats->mu0 = ans[0];  stats->gamma0 = ans[1];
			best_logLik = logLik;
		}

		if(verbose >= 10) printf("      (pnode:%d, tnode:%d, pparent:%d) => status=%d, loglik=%f (mu=%f, gamma=%f)\n", pnode, tnode, pParents[i], dag->getStatus(pParents[i], tnode), logLik, ans[0], ans[1]);
	}
	if(verbose >= 20) stats->print(stdout);

	if(best_logLik == -DBL_MAX) return -1;
	return 0;
}

int BetaBinLikelihood_PriorSetter::set(BetaBinomialCube *cube, int pnode, int tvector[], BetaBinStats *stats, int verbose){

	bool isAllZero = true;
	for(int i=0; i<cube->nTDim; i++) if(tvector[i] != 0) isAllZero = false;

	stats->mu0 = cube->root_priorMean;  stats->gamma0 = cube->root_priorSize;
	if(pnode == 0 && isAllZero) return 0;

	double best_logLik = -DBL_MAX;
	double ans[2];

	// rollup along the treatment
	int *tParent = new int[cube->nTDim]; for(int i=0; i<cube->nTDim; i++) tParent[i] = tvector[i];
	for(int i=0; i<cube->nTDim; i++){

		if(tvector[i] == 0) continue;
		tParent[i] = cube->cubeProcessor->CubeDimension[i]->getParent(tvector[i]);

		if(cube->getStatus(pnode, tParent) <= 0) continue;

		ans[0] = cube->getMean(pnode, tParent);
		ans[1] = cube->getSize(pnode, tParent);
		double logLik = logLik_BetaBinomial_mean_size(stats, ans, false);

		if(logLik > best_logLik){
			stats->mu0 = ans[0];  stats->gamma0 = ans[1];
			best_logLik = logLik;
		}

		if(verbose >= 10){
			printf("      (pnode:%d, tnode:(", pnode);
			print_row_vector(tvector,cube->nTDim,stdout,",");
			printf("), tparent:(");
			print_row_vector(tParent,cube->nTDim,stdout,",");
			printf(")) => status=%d, loglik=%f (mu=%f, gamma=%f)\n", cube->getStatus(pnode, tParent), logLik, ans[0], ans[1]);
		}
		tParent[i] = tvector[i];
	}

	// rollup along the population
	int num = cube->pDagIndex->numNextNodes(pnode);
	if(num == 0 && pnode != 0) STOP2("pNode %d has no parents", pnode);
	const int *pParents = cube->pDagIndex->getNextNodes(pnode);
	for(int i=0; i<num; i++){

		if(cube->getStatus(pParents[i], tvector) <= 0) continue;

		ans[0] = cube->getMean(pParents[i], tvector);
		ans[1] = cube->getSize(pParents[i], tvector);
		double logLik = logLik_BetaBinomial_mean_size(stats, ans, false);

		if(logLik > best_logLik){
			stats->mu0 = ans[0];  stats->gamma0 = ans[1];
			best_logLik = logLik;
		}

		if(verbose >= 10){
			printf("      (pnode:%d, tnode:(", pnode);
			print_row_vector(tvector,cube->nTDim,stdout,",");
			printf("), pparent:%d) => status=%d, loglik=%f (mu=%f, gamma=%f)\n", pParents[i], cube->getStatus(pParents[i], tvector), logLik, ans[0], ans[1]);
		}
	}
	if(verbose >= 20) stats->print(stdout);

	if(best_logLik == -DBL_MAX) return -1;

	delete[] tParent;
	return 0;
}


double pbetabinom_MeanSize(int nPositives, int nTrials, double mean, double size, double *prevProb){
	const double alpha = mean * size, beta = size - alpha;
	if(alpha <= 0 || beta <= 0 || size <= 0 || isnan(alpha) || isnan(beta)) STOP4("alpha=%.10g, beta=%.10g, size=%.10g", alpha, beta, size);

	double logP = lgammafn(beta+nTrials) + lgammafn(size) - lgammafn(beta) - lgammafn(size+nTrials);

	double sum = 1.0;
	double Tr  = 1.0;
	double bn  = -beta - nTrials;

	for(int i=0; i<nPositives; i++){
		double r  =  i;
		double rp = (i+1);
		Tr = Tr * ((r+alpha)*(r-nTrials))/(rp*(bn+rp));
		sum = sum+Tr;
		if(prevProb != NULL && i == (nPositives-2)) (*prevProb) = sum;
	}

	if(prevProb != NULL){
		if(nPositives == 0)      (*prevProb) = 0;
		else if(nPositives == 1) (*prevProb) = exp(logP);
		else                     (*prevProb) = exp(logP + log(*prevProb));
	}

	logP = logP + log(sum);

	return exp(logP);
}


double logLik_BetaBinomial_mean_size(const BetaBinStats* stats, const double *param, bool usePrior, bool ignoreConst){
	const double mu = param[0], gamma = param[1];
	const double alpha = mu * gamma, beta = gamma - alpha;
	if(alpha <= 0 || gamma <= 0 || isnan(alpha) || isnan(gamma)) STOP3("alpha=%.10g, gamma=%.10g", alpha, gamma);
	if(!stats->ready) STOP_HERE("!stats->ready");
	double ans = 0;
	// prior
	if(usePrior){
		if(stats->mu0 <= 0 || stats->mu0 >= 1 || stats->gamma0 <= 0 || isnan(stats->mu0) || isnan(stats->gamma0)) STOP3("mu0=%.10g, gamma0=%.10g", stats->mu0, stats->gamma0);
		double err = LOGIT(mu) - LOGIT(stats->mu0);
		ans -= err*err / (2*stats->var_mu);
		err = log(gamma) - log(stats->gamma0);
		ans -= err*err / (2*stats->var_gamma);
	}
	// small observations
	double small_part = 0;
	for(int j=0; j<stats->obsCountPos_length; j++) small_part += stats->obsCountPos[j] * log(alpha + j);
	for(int j=0; j<stats->obsCountNeg_length; j++) small_part += stats->obsCountNeg[j] * log(beta  + j);
	for(int j=0; j<stats->obsCount_length; j++)    small_part -= stats->obsCount[j]    * log(gamma + j);
	int *x=stats->largeObsPos, *n=stats->largeObsSize;
	// large observations
	double large_part = 0;
	for(int i=0; i<stats->nSampledLargeObs; i++){
		large_part += lgammafn(x[i]+alpha) + lgammafn(n[i]-x[i]+beta) + lgammafn(gamma)
				      - lgammafn(alpha) - lgammafn(beta) - lgammafn(n[i]+gamma);
	}
	// combine the two parts
	ans += small_part + large_part; // TODO: Weigh the two parts correctly

	if(!ignoreConst){
		// do not ignore the constants
		if(stats->sum_logConst_forSmall == -DBL_MAX || stats->sum_logConst_forLarge == -DBL_MAX || usePrior) STOP_HERE("Please initialize BetaBinStats with ignoreConst=false and do not use prior");
		ans += stats->sum_logConst_forSmall + stats->sum_logConst_forLarge;  // TODO: Weigh the two parts correctly
	}

	return ans;
}

void d_LogLik_BetaBinomial_mean_size(double* gr, const BetaBinStats* stats, const double *param, int verbose){
	const double mu = param[0], gamma = param[1];
	const double alpha = mu * gamma, beta = gamma - alpha;
	if(alpha <= 0 || gamma <= 0 || isnan(alpha) || isnan(gamma)) STOP3("alpha=%.10g, gamma=%.10g", alpha, gamma);
	if(stats->mu0 <= 0 || stats->mu0 >= 1 || stats->gamma0 <= 0 || isnan(stats->mu0) || isnan(stats->gamma0)) STOP3("mu0=%.10g, gamma0=%.10g", stats->mu0, stats->gamma0);
	if(!stats->ready) STOP_HERE("!stats->ready");
	// prior
	double mu_prior_part = - (LOGIT(mu) - LOGIT(stats->mu0)) * ((1/mu) - (1/(1-mu))) / stats->var_mu;
	double gamma_prior_part = - (log(gamma) - log(stats->gamma0)) / (gamma * stats->var_gamma);
	// small observations
	double mu_small_part = 0, gamma_small_part = 0;
	int *m_plus = stats->obsCountPos, *m_neg = stats->obsCountNeg, *m = stats->obsCount;
	for(int j=0; j<stats->obsCountPos_length; j++){
		mu_small_part    += m_plus[j] / (mu    + (j/gamma));
		gamma_small_part += m_plus[j] / (gamma + (j/mu));
	}
	for(int j=0; j<stats->obsCountNeg_length; j++){
		mu_small_part    -= m_neg[j] / (1 - mu + (j/gamma));
		gamma_small_part += m_neg[j] / (gamma  + (j/(1-mu)));
	}
	for(int j=0; j<stats->obsCount_length; j++){
		gamma_small_part -= m[j] / (gamma+j);
	}
	int *x=stats->largeObsPos, *n=stats->largeObsSize;
	// large observations
	double mu_large_part = 0, gamma_large_part = 0;
	for(int i=0; i<stats->nSampledLargeObs; i++){
		double Psi_x_alpha = digamma(x[i]+alpha),  Psi_n_x_beta = digamma(n[i]-x[i]+beta),
			   Psi_alpha   = digamma(alpha),       Psi_beta     = digamma(beta);
		mu_large_part    += gamma * (Psi_x_alpha - Psi_n_x_beta - Psi_alpha + Psi_beta);
		gamma_large_part += Psi_x_alpha*mu + Psi_n_x_beta*(1-mu) + digamma(gamma)
		    			  - Psi_alpha * mu - Psi_beta*(1-mu)     - digamma(n[i]+gamma);
	}
	// combine the parts
	gr[0] = mu_prior_part + mu_small_part + mu_large_part; // TODO: Weigh the two parts correctly
	gr[1] = gamma_prior_part + gamma_small_part + gamma_large_part; // TODO: Weigh the two parts correctly
	if(verbose >= 10){
		printf("      gradient[0] = %.16g + %.16g + %.16g\n", mu_prior_part, mu_small_part, mu_large_part);
		printf("      gradient[1] = %.16g + %.16g + %.16g\n", gamma_prior_part, gamma_small_part, gamma_large_part);
	}
}

void d2_LogLik_BetaBinomial_mean_size(double* H, const BetaBinStats* stats, const double *param, int verbose){
	const double mu = param[0], gamma = param[1];
	const double alpha = mu * gamma, beta = gamma - alpha;
	if(alpha <= 0 || gamma <= 0 || isnan(alpha) || isnan(gamma)) STOP3("alpha=%.10g, gamma=%.10g", alpha, gamma);
	if(stats->mu0 <= 0 || stats->mu0 >= 1 || stats->gamma0 <= 0 || isnan(stats->mu0) || isnan(stats->gamma0)) STOP3("mu0=%.10g, gamma0=%.10g", stats->mu0, stats->gamma0);
	if(!stats->ready) STOP_HERE("!stats->ready");
	// prior
	double mu_prior_part = - (LOGIT(mu) - LOGIT(stats->mu0)) * (-1/SQR(mu) + 1/SQR(1-mu)) / stats->var_mu
						   - SQR(1/mu + 1/(1-mu)) / stats->var_mu;
	double gamma_prior_part = - (1 - log(gamma) + log(stats->gamma0)) / (SQR(gamma) * stats->var_gamma);
	// small observations
	double mu_small_part = 0, gamma_small_part = 0, muga_small_part = 0;
	int *m_plus = stats->obsCountPos, *m_neg = stats->obsCountNeg, *m = stats->obsCount;
	for(int j=0; j<stats->obsCountPos_length; j++){
		mu_small_part    -= m_plus[j] / SQR(mu    + (j/gamma));
		gamma_small_part -= m_plus[j] / SQR(gamma + (j/mu));
		muga_small_part  += m_plus[j] * j / SQR(alpha + j);
	}
	for(int j=0; j<stats->obsCountNeg_length; j++){
		mu_small_part    -= m_neg[j] / SQR(1 - mu + (j/gamma));
		gamma_small_part -= m_neg[j] / SQR(gamma  + (j/(1-mu)));
		muga_small_part  -= m_neg[j] * j / SQR(beta + j);
	}
	for(int j=0; j<stats->obsCount_length; j++){
		gamma_small_part += m[j] / SQR(gamma+j);
	}
	int *x=stats->largeObsPos, *n=stats->largeObsSize;
	// large observations
	double mu_large_part = 0, gamma_large_part = 0, muga_large_part = 0;
	for(int i=0; i<stats->nSampledLargeObs; i++){
		double dPsi_x_alpha = trigamma(x[i]+alpha),  dPsi_n_x_beta = trigamma(n[i]-x[i]+beta),
			   dPsi_alpha   = trigamma(alpha),       dPsi_beta     = trigamma(beta);
		mu_large_part    += SQR(gamma) * (dPsi_x_alpha + dPsi_n_x_beta - dPsi_alpha - dPsi_beta);
		gamma_large_part += dPsi_x_alpha*SQR(mu) + dPsi_n_x_beta*SQR(1-mu) + trigamma(gamma)
				          - dPsi_alpha * SQR(mu) - dPsi_beta*SQR(1-mu)     - trigamma(n[i]+gamma);
		muga_large_part  += digamma(x[i]+alpha) - digamma(n[i]-x[i]+beta) - digamma(alpha) + digamma(beta)
		 			       + dPsi_x_alpha*alpha - dPsi_n_x_beta*beta - dPsi_alpha*alpha + dPsi_beta*beta;
	}
	// combine the parts
	H[0] = mu_prior_part + mu_small_part + mu_large_part; // TODO: Weigh the two parts correctly
	H[1] = muga_small_part + muga_large_part; // TODO: Weigh the two parts correctly
	H[2] = H[1];
	H[3] = gamma_prior_part + gamma_small_part + gamma_large_part; // TODO: Weigh the two parts correctly

	if(verbose >= 10){
		printf("      hessian[0] = %.16g + %.16g + %.16g\n", mu_prior_part, mu_small_part, mu_large_part);
		printf("      hessian[1] = %.16g + %.16g\n", muga_small_part, muga_large_part);
		printf("      hessian[3] = %.16g + %.16g + %.16g\n", gamma_prior_part, gamma_small_part, gamma_large_part);
	}
}


double logLik_BetaBinomial_logitMean_logSize(const BetaBinStats* stats, const double *param, bool usePrior){
	const double logit_mu = param[0], log_gamma = param[1], mu = INV_LOGIT(param[0]), gamma = exp(param[1]);
	const double alpha = mu * gamma, beta = gamma - alpha;
	if(alpha <= 0 || gamma <= 0 || isnan(alpha) || isnan(gamma)) STOP3("alpha=%.10g, gamma=%.10g", alpha, gamma);
	if(!stats->ready) STOP_HERE("!stats->ready");
	double ans = 0;
	// prior
	if(usePrior){
		if(stats->mu0 <= 0 || stats->mu0 >= 1 || stats->gamma0 <= 0 || isnan(stats->mu0) || isnan(stats->gamma0)) STOP3("mu0=%.10g, gamma0=%.10g", stats->mu0, stats->gamma0);
		double err = logit_mu - LOGIT(stats->mu0);
		ans -= err*err / (2*stats->var_mu);
		err = log_gamma - log(stats->gamma0);
		ans -= err*err / (2*stats->var_gamma);
	}
	// small observations
	double small_part = 0;
	for(int j=0; j<stats->obsCountPos_length; j++) small_part += stats->obsCountPos[j] * log(alpha + j);
	for(int j=0; j<stats->obsCountNeg_length; j++) small_part += stats->obsCountNeg[j] * log(beta  + j);
	for(int j=0; j<stats->obsCount_length; j++)    small_part -= stats->obsCount[j]    * log(gamma + j);
	int *x=stats->largeObsPos, *n=stats->largeObsSize;
	// large observations
	double large_part = 0;
	for(int i=0; i<stats->nSampledLargeObs; i++){
		large_part += lgammafn(x[i]+alpha) + lgammafn(n[i]-x[i]+beta) + lgammafn(gamma)
				      - lgammafn(alpha) - lgammafn(beta) - lgammafn(n[i]+gamma);
	}
	// combine the two parts
	ans += small_part + large_part; // TODO: Weigh the two parts correctly
	return ans;
}


void d_d2_LogLik_BetaBinomial_logitMean_logSize(double* gr, double* H, const BetaBinStats* stats, const double *param, int verbose){
	const double logit_mu = param[0], log_gamma = param[1], mu = INV_LOGIT(param[0]), gamma = exp(param[1]);
	const double alpha = mu * gamma, beta = gamma - alpha;
	const int *m_plus = stats->obsCountPos, *m_neg = stats->obsCountNeg, *m = stats->obsCount;
	const int *x=stats->largeObsPos, *n=stats->largeObsSize;

	if(alpha <= 0 || gamma <= 0 || isnan(alpha) || isnan(gamma)) STOP3("alpha=%.10g, gamma=%.10g", alpha, gamma);
	if(stats->mu0 <= 0 || stats->mu0 >= 1 || stats->gamma0 <= 0 || isnan(stats->mu0) || isnan(stats->gamma0)) STOP3("mu0=%.10g, gamma0=%.10g", stats->mu0, stats->gamma0);
	if(!stats->ready) STOP_HERE("!stats->ready");

	// Prior
	double gr_mu_prior_part    = - (logit_mu - LOGIT(stats->mu0)) / stats->var_mu;
	double gr_gamma_prior_part = - (log_gamma - log(stats->gamma0)) / stats->var_gamma;

	double H_mu_prior_part    = -1 / stats->var_mu;
	double H_gamma_prior_part = -1 / stats->var_gamma;

	// Small observations
	double gr_mu_small_part = 0, gr_gamma_small_part = 0;
	double H_mu_small_part = 0,  H_gamma_small_part = 0,  H_muga_small_part = 0;
	for(int j=0; j<stats->obsCountPos_length; j++){
		gr_mu_small_part    += m_plus[j] / (mu    + (j/gamma));
		gr_gamma_small_part += m_plus[j] / (gamma + (j/mu));

		H_mu_small_part    -= m_plus[j] / SQR(mu    + (j/gamma));
		H_gamma_small_part -= m_plus[j] / SQR(gamma + (j/mu));
		H_muga_small_part  += m_plus[j] * j / SQR(alpha + j);
	}
	for(int j=0; j<stats->obsCountNeg_length; j++){
		gr_mu_small_part    -= m_neg[j] / (1 - mu + (j/gamma));
		gr_gamma_small_part += m_neg[j] / (gamma  + (j/(1-mu)));

		H_mu_small_part    -= m_neg[j] / SQR(1 - mu + (j/gamma));
		H_gamma_small_part -= m_neg[j] / SQR(gamma  + (j/(1-mu)));
		H_muga_small_part  -= m_neg[j] * j / SQR(beta + j);
	}
	for(int j=0; j<stats->obsCount_length; j++){
		gr_gamma_small_part -= m[j] / (gamma+j);

		H_gamma_small_part += m[j] / SQR(gamma+j);
	}

	// Large observations
	double gr_mu_large_part = 0, gr_gamma_large_part = 0;
	double H_mu_large_part = 0, H_gamma_large_part = 0, H_muga_large_part = 0;
	for(int i=0; i<stats->nSampledLargeObs; i++){
		double Psi_x_alpha = digamma(x[i]+alpha),    Psi_n_x_beta = digamma(n[i]-x[i]+beta),
			   Psi_alpha   = digamma(alpha),         Psi_beta     = digamma(beta),
			   dPsi_x_alpha = trigamma(x[i]+alpha),  dPsi_n_x_beta = trigamma(n[i]-x[i]+beta),
			   dPsi_alpha   = trigamma(alpha),       dPsi_beta     = trigamma(beta);

		gr_mu_large_part    += gamma * (Psi_x_alpha - Psi_n_x_beta - Psi_alpha + Psi_beta);
		gr_gamma_large_part += Psi_x_alpha*mu + Psi_n_x_beta*(1-mu) + digamma(gamma)
						     - Psi_alpha * mu - Psi_beta*(1-mu)     - digamma(n[i]+gamma);

		H_mu_large_part    += SQR(gamma) * (dPsi_x_alpha + dPsi_n_x_beta - dPsi_alpha - dPsi_beta);
		H_gamma_large_part += dPsi_x_alpha*SQR(mu) + dPsi_n_x_beta*SQR(1-mu) + trigamma(gamma)
						    - dPsi_alpha * SQR(mu) - dPsi_beta*SQR(1-mu)     - trigamma(n[i]+gamma);
		H_muga_large_part  += Psi_x_alpha - Psi_n_x_beta - Psi_alpha + Psi_beta
						    + dPsi_x_alpha*alpha - dPsi_n_x_beta*beta - dPsi_alpha*alpha + dPsi_beta*beta;
	}

	double gr_mu_obs    = gr_mu_small_part + gr_mu_large_part;        // TODO: Weigh the two parts correctly
	double gr_gamma_obs = gr_gamma_small_part + gr_gamma_large_part;  // TODO: Weigh the two parts correctly
	double H_mu_obs     = H_mu_small_part + H_mu_large_part;          // TODO: Weigh the two parts correctly
	double H_gamma_obs  = H_gamma_small_part + H_gamma_large_part;    // TODO: Weigh the two parts correctly
	double H_muga_obs   = H_muga_small_part + H_muga_large_part;      // TODO: Weigh the two parts correctly

	// Convert to the derivative w.r.t. logit(mu) and log(gamma)
	double gr_logitMu_obs_part  = mu * (1-mu) * gr_mu_obs;
	double gr_logGamma_obs_part = gamma * gr_gamma_obs;

	double H_logitMu_obs_part  = mu * (1-mu) * (H_mu_obs*mu*(1-mu) + gr_mu_obs*(1 - 2*mu));
	double H_logGamma_obs_part = gamma * (H_gamma_obs * gamma + gr_gamma_obs);
	double H_logitMu_logGamma  = mu * (1-mu) * gamma * H_muga_obs;

	// combine the parts
	gr[0] = gr_mu_prior_part    + gr_logitMu_obs_part;
	gr[1] = gr_gamma_prior_part + gr_logGamma_obs_part;

	H[0] = H_mu_prior_part + H_logitMu_obs_part;
	H[1] = H_logitMu_logGamma;
	H[2] = H[1];
	H[3] = H_gamma_prior_part + H_logGamma_obs_part;

	if(verbose >= 10){
		printf("      gradient[0] = %.16g + %.16g\n", gr_mu_prior_part, gr_logitMu_obs_part);
		printf("      gradient[1] = %.16g + %.16g\n", gr_gamma_prior_part, gr_logGamma_obs_part);
		printf("       hessian[0] = %.16g + %.16g\n", H_mu_prior_part, H_logitMu_obs_part);
		printf("       hessian[1] = %.16g\n", H_logitMu_logGamma);
		printf("       hessian[3] = %.16g + %.16g\n", H_gamma_prior_part, H_logGamma_obs_part);
	}
}

#define MIN_MU 	  1e-20
#define MIN_GAMMA 1e-20

typedef struct{
	BetaBinStats* stats;
	double param[2];
	double dir[2];
} BetaBinLineSearchStruct;

double negLogLik_BetaBinomial_logitMean_logSize(double step, void* obj){
	BetaBinLineSearchStruct* data = (BetaBinLineSearchStruct*)obj;
	double x[2];
	x[0] = data->param[0] - step * data->dir[0];
	x[1] = data->param[1] - step * data->dir[1];
	return - logLik_BetaBinomial_logitMean_logSize(data->stats, x, true);
}

void put_in_range(double *step, double param, double delta, double min, double max){
	if(param - (*step)*delta > max){
		(*step) = (param - max)/delta;
	}else if(param - (*step)*delta < min){
		(*step) = (param - min)/delta;
	}
}


// if param[0] == 0 then use mu0
// if param[1] == 0 then use gamma0
int postMode_BetaBinomial_mean_size(
	double *param, int *nIter, const BetaBinStats *stats, double epsilon, double stepSize,
	int maxIter, int nLnsrchSteps,
	double minMean, double maxMean, double minSize, double maxSize,
	int verbose, int debug
){
	double gr[2]; // gradient
	double H[4];  // hessian
	if(nIter != NULL) (*nIter) = 0;
	if(stepSize <= 0 || stepSize > 1) STOP_HERE("stepSize <= 0 || stepSize > 1.  Set to 1.");
	if(stats->mu0 <= 0 || stats->mu0 >= 1) STOP_HERE("stats->mu0 <= 0 || stats->mu0 >= 1");
	if(stats->gamma0 <= 0) STOP_HERE("stats->gamma0 <= 0");

	double minLogitMean=LOGIT(minMean);
	double maxLogitMean=LOGIT(maxMean);
	double minLogSize=log(minSize);
	double maxLogSize=log(maxSize);

	double mu    = (param[0] == 0 ? stats->mu0    : param[0]),
		   gamma = (param[1] == 0 ? stats->gamma0 : param[1]);

	if(verbose > 0){
		printf("Find the posterior mode of Beta Binomial (logit(mu), log(gamma))\n");
		printf("  Init mu = %.16g, gamma = %.16g\n", mu, gamma);
	}

	if(mu <= 0 || mu >= 1) STOP_HERE("mu <= 0 || mu >= 1");
	if(gamma <= 0) STOP_HERE("gamma <= 0");

	param[0] = LOGIT(mu);	param[0] = FORCE_RANGE(param[0], minLogitMean, maxLogitMean);
	param[1] = log(gamma);  param[1] = FORCE_RANGE(param[1], minLogSize, maxLogSize);

	double best_loglik = -DBL_MAX, best_logitMu = LOGIT(stats->mu0), best_logGamma = log(stats->gamma0);

	// Newton's method
	double prev_mu=0, prev_cv=0;
	BetaBinLineSearchStruct lnsrch_obj;

	for(int iter=0; iter<maxIter; iter++){

		if(verbose > 1){
			double loglik2 = logLik_BetaBinomial_logitMean_logSize(stats, param, true);
			printf("  Iteration %d: mu=%.16g, gamma=%.16g, loglik=%.16g\n", iter, INV_LOGIT(param[0]), exp(param[1]), loglik2);
		}

		d_d2_LogLik_BetaBinomial_logitMean_logSize(gr, H, stats, param, verbose);

		if(verbose > 2){
			printf("    gr: %.16g\t%.16g\n", gr[0], gr[1]);
			printf("     H: %.16g\t%.16g\n", H[0], H[1]);
			printf("        %.16g\t%.16g\n", H[2], H[3]);
		}

		if(nIter != NULL) (*nIter) ++;

		double detH = H[0]*H[3] - H[1]*H[2];
		if(detH < 1e-50 && detH > -1e-50){
			if(verbose > 0) printf("  FAILED!!  det(H) = %.16g\n", detH);
			param[0] = INV_LOGIT(best_logitMu); param[1] = exp(best_logGamma);
			return 1;
		}

		double delta0 = ( H[3]*gr[0] - H[2]*gr[1])/detH;
		double delta1 = (-H[1]*gr[0] + H[0]*gr[1])/detH;

		if(verbose > 2){
			printf("    delta: %.16g\t%.16g\n", delta0, delta1);
		}

		if(nLnsrchSteps > 0){

			lnsrch_obj.stats = const_cast<BetaBinStats*>(stats);
			lnsrch_obj.param[0] = param[0]; lnsrch_obj.param[1] = param[1];
			lnsrch_obj.dir[0]   = delta0;   lnsrch_obj.dir[1]   = delta1;

			// set a proper range
			double min_step = -stepSize;
			double max_step = stepSize;
			put_in_range(&min_step, param[0], delta0, minLogitMean, maxLogitMean);
			put_in_range(&max_step, param[0], delta0, minLogitMean, maxLogitMean);
			put_in_range(&min_step, param[1], delta1, minLogSize, maxLogSize);
			put_in_range(&max_step, param[1], delta1, minLogSize, maxLogSize);
			if(min_step >= max_step){
				param[0] = INV_LOGIT(best_logitMu); param[1] = exp(best_logGamma);
				if(verbose > 0) printf("  FAILED in linesearch!!  min_step (%.16g) >= max_step (%.16g)\n", min_step, max_step);
				return 3;
			}

			double step = Brent_fmin(min_step, max_step, negLogLik_BetaBinomial_logitMean_logSize, &lnsrch_obj, 1e-6 /*, nLnsrchSteps*/);

			if(verbose > 1){
				printf("    line search: step=%.16g\n", step);
			}

			param[0] -= step * delta0;
			param[1] -= step * delta1;
			param[0] = FORCE_RANGE(param[0], minLogitMean, maxLogitMean);
			param[1] = FORCE_RANGE(param[1], minLogSize, maxLogSize);

			double loglik = logLik_BetaBinomial_logitMean_logSize(stats, param, true);

			if(loglik > best_loglik){
				best_logitMu  = param[0];
				best_logGamma = param[1];
				best_loglik = loglik;
			}

		}else{
			param[0] -= stepSize * delta0;
			param[1] -= stepSize * delta1;
			param[0] = FORCE_RANGE(param[0], minLogitMean, maxLogitMean);
			param[1] = FORCE_RANGE(param[1], minLogSize, maxLogSize);
		}

		mu    = INV_LOGIT(param[0]);
		gamma = exp(param[1]);

		if(param[0] < -500 || param[0] > 500 || param[1] < -500 || param[1] > 500){
			param[0] = INV_LOGIT(best_logitMu); param[1] = exp(best_logGamma);
			if(verbose > 0) printf("  FAILED in linesearch!!  Out of range\n");
			return 3;
		}

		// Stopping criterion
		double cv = sqrt(VAR_BETABIN_MEAN_SIZE(mu,gamma)) / mu;
		if(fabs(mu - prev_mu) <= epsilon && fabs(cv - prev_cv) <= epsilon){
			param[0] = mu; param[1] = gamma;
			if(verbose > 2) printf("  Stop:  mu=%.16g  prev_mu=%.16g  cv=%.16g  prev_cv=%.16g\n", mu, prev_mu, cv, prev_cv);
			if(verbose > 0) printf("  OUTPUT:  mean=%.16g  size=%.16g\n", param[0], param[1]);
			return 0;
		}
		prev_mu = mu; prev_cv = cv;
	}
	param[0] = INV_LOGIT(best_logitMu); param[1] = exp(best_logGamma);
	if(verbose > 0) printf("  FAILED!!  Maximum number of iteration reached\n");
	return 2;
}

int postMode_BetaBinomial_mean_size_old(
	double *param, int *nIter, const BetaBinStats *stats, double epsilon, double stepSize, int maxIter, bool halfStepsizing,
	int verbose, int debug
){
	double gr[2]; // gradient
	double H[4];  // hessian
	double temp[2];
	if(nIter != NULL) (*nIter) = 0;
	if(stepSize <= 0 || stepSize > 1) STOP_HERE("stepSize <= 0 || stepSize > 1.  Set to 1.");

	if(verbose > 0){
		printf("Find the posterior mode of Beta Binomial\n");
	}

	param[0] = stats->mu0; param[1] = stats->gamma0;
	if(param[0] < MIN_MU) param[0] = MIN_MU;
	if(param[0] > 1 - MIN_MU) param[0] = 1 - MIN_MU;
	if(param[1] < MIN_GAMMA) param[1] = MIN_GAMMA;
	// Newton's method
	double prev_mu=0, prev_cv=0, loglik=-DBL_MAX;
	if(halfStepsizing){
		loglik = logLik_BetaBinomial_mean_size(stats, param, true);
	}
	int iter = 0;
	for(iter=0; iter<maxIter; iter++){

		if(verbose > 1){
			double loglik2 = logLik_BetaBinomial_mean_size(stats, param, true);
			printf("  Iteration %d: mu=%.16g, gamma=%.16g, loglik=%.16g\n", iter, param[0], param[1], loglik2);
		}

		d_LogLik_BetaBinomial_mean_size(gr, stats, param, verbose);
		d2_LogLik_BetaBinomial_mean_size(H, stats, param, verbose);

		if(verbose > 2){
			printf("    gr: %.16g\t%.16g\n", gr[0], gr[1]);
			printf("     H: %.16g\t%.16g\n", H[0], H[1]);
			printf("        %.16g\t%.16g\n", H[2], H[3]);
		}

		if(nIter != NULL) (*nIter) ++;

		double detH = H[0]*H[3] - H[1]*H[2];
		if(detH == 0) return 1;

		double delta0 = ( H[3]*gr[0] - H[2]*gr[1])/detH;
		double delta1 = (-H[1]*gr[0] + H[0]*gr[1])/detH;

		if(verbose > 2){
			printf("    delta: %.16g\t%.16g\n", delta0, delta1);
		}

		if(halfStepsizing){
			double newloglik = -DBL_MAX, step = stepSize;
			for(int k=0; k<10; k++){
				temp[0] = param[0] - step * delta0;
				temp[1] = param[1] - step * delta1;

				if(temp[0] < MIN_MU) temp[0] = MIN_MU;
				if(temp[0] > 1 - MIN_MU) temp[0] = 1 - MIN_MU;
				if(temp[1] < MIN_GAMMA) temp[1] = MIN_GAMMA;

				newloglik = logLik_BetaBinomial_mean_size(stats, temp, true);

				if(verbose >= 5){
					printf("    stepSize=%f:  mu=%.16g, gamma=%.16g, loglik=%.16g\n", step, temp[0], temp[1], newloglik);
				}

				if(newloglik > loglik) break;
				else                   step /= 2;
			}
			param[0] = temp[0];
			param[1] = temp[1];
			loglik = newloglik;
		}else{
			param[0] -= stepSize * delta0;
			param[1] -= stepSize * delta1;

			if(param[0] < MIN_MU) param[0] = MIN_MU;
			if(param[0] > 1 - MIN_MU) param[0] = 1 - MIN_MU;
			if(param[1] < MIN_GAMMA) param[1] = MIN_GAMMA;
		}

		// Stopping criterion
		double cv = sqrt(VAR_BETABIN_MEAN_SIZE(param[0],param[1])) / param[0];
		if(fabs(param[0] - prev_mu) <= epsilon && fabs(cv - prev_cv) <= epsilon) return 0;
		prev_mu = param[0]; prev_cv = cv;
	}
	return 2;
}

//-----------------------------------------------------------------------------
//      Data Structures
//-----------------------------------------------------------------------------

BetaBinStats::BetaBinStats(
	const int *pos, const int *size, const int length, const int threshold, const int nSamplesForLarge,
	const double mu0, const double var_mu, const double gamma0, const double var_gamma
){
	this->threshold = threshold;  this->mu0 = mu0;  this->var_mu = var_mu;
	this->gamma0 = gamma0; this->var_gamma = var_gamma;
	sum_logConst_forSmall = -DBL_MAX;
	sum_logConst_forLarge = -DBL_MAX;
	nSmallObs = 0;
	nLargeObs = 0;
	obsCount_length = 0;
	obsCountPos_length = 0;
	obsCountNeg_length = 0;
	this->nSamplesForLarge = nSamplesForLarge;
	for(int i=0; i<length; i++){
		if(size[i] < 0 || pos[i] < 0 || pos[i] > size[i]) STOP_HERE("input data error");
		if(size[i] <= threshold){
			nSmallObs++;
			if(size[i] > obsCount_length)    obsCount_length    = size[i];
			if(pos[i]  > obsCountPos_length) obsCountPos_length = pos[i];
			if((size[i]-pos[i]) > obsCountNeg_length) obsCountNeg_length = size[i]-pos[i];
		}else{
			nLargeObs++;
		}
	}
	nSampledLargeObs = nSamplesForLarge >= nLargeObs ? nLargeObs : nSamplesForLarge;
	obsCount = new int[obsCount_length];		for(int i=0; i<obsCount_length; i++) obsCount[i] = 0;
	obsCountPos = new int[obsCountPos_length];	for(int i=0; i<obsCountPos_length; i++) obsCountPos[i] = 0;
	obsCountNeg = new int[obsCountNeg_length]; 	for(int i=0; i<obsCountNeg_length; i++) obsCountNeg[i] = 0;
	largeObsPos  = new int[nSampledLargeObs];
	largeObsSize = new int[nSampledLargeObs];

	int k=0;
	for(int i=0; i<length; i++){
		if(size[i] <= threshold){
			if(size[i] > 0) obsCount[size[i]-1]++;
			if(pos[i]  > 0) obsCountPos[pos[i]-1]++;
			if(size[i]-pos[i] > 0) obsCountNeg[size[i]-pos[i]-1]++;
		}else{
			// Reservoir Sampling
			if(k < nSampledLargeObs){
				largeObsPos[k] = pos[i];
				largeObsSize[k] = size[i];
			}else{
				// TODO: Use a better random # generator
				// When k > RAND_MAX, the result won't be random enough
				int r = (int)floor(k * ((double)rand()) / ((double)RAND_MAX));
				if(r < nSampledLargeObs){
					largeObsPos[r] = pos[i];
					largeObsSize[r] = size[i];
				}
			}
			k++;
		}
	}
	for(int j=obsCount_length-2; j>=0; j--) obsCount[j] += obsCount[j+1];
	for(int j=obsCountPos_length-2; j>=0; j--) obsCountPos[j] += obsCountPos[j+1];
	for(int j=obsCountNeg_length-2; j>=0; j--) obsCountNeg[j] += obsCountNeg[j+1];

	ready = true;
}

BetaBinStats::BetaBinStats(
	const int threshold, const int nSamplesForLarge,
	const int nLargeObs, const int obsCount_length, const int obsCountPos_length, const int obsCountNeg_length,
	const double var_mu, const double var_gamma, const bool ignoreConst
){
	this->threshold = threshold;  this->mu0 = -1;  this->var_mu = var_mu;
	this->gamma0 = -1; this->var_gamma = var_gamma;
	sum_logConst_forSmall = -DBL_MAX;
	sum_logConst_forLarge = -DBL_MAX;

	if(!ignoreConst){
		sum_logConst_forSmall = 0;
		sum_logConst_forLarge = 0;
	}

	nSmallObs = 0;
	this->nLargeObs = 0;
	this->nSamplesForLarge = nSamplesForLarge;

	this->obsCount_length = obsCount_length;
	this->obsCountPos_length = obsCountPos_length;
	this->obsCountNeg_length = obsCountNeg_length;

	nSampledLargeObs = nSamplesForLarge >= nLargeObs ? nLargeObs : nSamplesForLarge;
	obsCount = new int[obsCount_length];		for(int i=0; i<obsCount_length; i++) obsCount[i] = 0;
	obsCountPos = new int[obsCountPos_length];	for(int i=0; i<obsCountPos_length; i++) obsCountPos[i] = 0;
	obsCountNeg = new int[obsCountNeg_length]; 	for(int i=0; i<obsCountNeg_length; i++) obsCountNeg[i] = 0;
	largeObsPos  = new int[nSampledLargeObs];
	largeObsSize = new int[nSampledLargeObs];

	ready = false;
	k = 0;
}

void BetaBinStats::add(int nPositives, int nTrials){
	if(ready) STOP_HERE("it's already finalized");
	if(nTrials < 0 || nPositives < 0 || nPositives > nTrials) STOP_HERE("input data error");
	if(nTrials <= threshold){
		if(sum_logConst_forSmall != -DBL_MAX) sum_logConst_forSmall += lchoose(nTrials, nPositives);
		nSmallObs++;
		if(nTrials > obsCount_length) STOP_HERE("nTrials > obsCount_length");
		if(nPositives > obsCountPos_length) STOP_HERE(nPositives > obsCountPos_length);
		if(nTrials-nPositives > obsCountNeg_length) STOP_HERE("nTrials-nPositives > obsCountNeg_length");
		if(nTrials > 0) obsCount[nTrials-1]++;
		if(nPositives  > 0) obsCountPos[nPositives-1]++;
		if(nTrials-nPositives > 0) obsCountNeg[nTrials-nPositives-1]++;
	}else{
		nLargeObs++;
		// Reservoir Sampling
		if(k < nSampledLargeObs){
			if(sum_logConst_forLarge != -DBL_MAX) sum_logConst_forLarge += lchoose(nTrials, nPositives);
			largeObsPos[k] = nPositives;
			largeObsSize[k] = nTrials;
		}else{
			// TODO: Use a better random # generator
			// When k > RAND_MAX, the result won't be random enough
			int r = (int)floor(k * ((double)rand()) / ((double)RAND_MAX));
			if(r < nSampledLargeObs){
				if(sum_logConst_forLarge != -DBL_MAX){
					sum_logConst_forLarge -= lchoose(largeObsSize[r], largeObsPos[r]);
					sum_logConst_forLarge += lchoose(nTrials, nPositives);
				}
				largeObsPos[r] = nPositives;
				largeObsSize[r] = nTrials;
			}
		}
		k++;
	}
}

void BetaBinStats::set(BetaBinStats *input[], int num){
	int total = 0;
	double *startIndex = new double[num+1];
	int *pick = new int[num];
	// handle small obs
	for(int j=0; j<obsCount_length; j++)    obsCount[j]    = 0;
	for(int j=0; j<obsCountPos_length; j++) obsCountPos[j] = 0;
	for(int j=0; j<obsCountNeg_length; j++) obsCountNeg[j] = 0;
	nSmallObs = 0; nLargeObs = 0;
	if(sum_logConst_forSmall != -DBL_MAX) sum_logConst_forSmall = 0;
	if(sum_logConst_forLarge != -DBL_MAX) sum_logConst_forLarge = 0;
	for(int i=0; i<num; i++){
		if(input[i] == NULL) STOP_HERE("input[i] == NULL");
		if(input[i]->ready) STOP_HERE("input[i]->ready");
		if(input[i]->obsCount_length > obsCount_length) STOP_HERE("input[i]->obsCount_length > obsCount_length");
		if(input[i]->obsCountPos_length > obsCountPos_length) STOP_HERE("input[i]->obsCountPos_length > obsCountPos_length");
		if(input[i]->obsCountNeg_length > obsCountNeg_length) STOP_HERE("input[i]->obsCountNeg_length > obsCountNeg_length");
		if(input[i]->threshold != threshold) STOP_HERE("input[i]->threshold != threshold");

		if(sum_logConst_forSmall != -DBL_MAX){
			if(input[i]->sum_logConst_forSmall == -DBL_MAX) STOP_HERE("");
			sum_logConst_forSmall += input[i]->sum_logConst_forSmall;
		}

		nSmallObs += input[i]->nSmallObs;
		for(int j=0; j<input[i]->obsCount_length; j++)    obsCount[j]    += input[i]->obsCount[j];
		for(int j=0; j<input[i]->obsCountPos_length; j++) obsCountPos[j] += input[i]->obsCountPos[j];
		for(int j=0; j<input[i]->obsCountNeg_length; j++) obsCountNeg[j] += input[i]->obsCountNeg[j];

		total += input[i]->nSampledLargeObs;
		nLargeObs += input[i]->nLargeObs;

		pick[i] = i;
		startIndex[i] = input[i]->nLargeObs;
	}
	for(int i=0; i<num; i++) startIndex[i] *= ((double)nSampledLargeObs)/nLargeObs;
	// handle large obs
	if(nLargeObs == nSampledLargeObs){
		k = 0;
		for(int i=0; i<num; i++){
			for(int j=0; j<input[i]->nSampledLargeObs; j++){
				largeObsPos[k]  = input[i]->largeObsPos[j];
				largeObsSize[k] = input[i]->largeObsSize[j];
				k++;
			}
		}
	}else if(nLargeObs < nSampledLargeObs){
		STOP_HERE("nLargeObs < nSampledLargeObs");
	}else{
		// randomize the order
		for(int i=0; i<num-1; i++){
			int r = i + rand() % (num-i);
			if(i != r){
				swap(pick[i], pick[r]);
				swap(startIndex[i], startIndex[r]);
			}
		}
		for(int i=1; i<num; i++) startIndex[i] += startIndex[i-1];
		for(int i=0; i<num; i++) startIndex[i] = floor(startIndex[i]);
		startIndex[num] = nSampledLargeObs;

		k=0;
		for(int i=0; i<num; i++){
			int nPoints = (int)startIndex[i+1] - k;
			if(nPoints >= input[i]->nSampledLargeObs){
				if(nPoints > input[i]->nSampledLargeObs){
					printf("Warning: nPoints = %d, input[%d]->nSampledLargeObs = %d\n", nPoints, i, input[i]->nSampledLargeObs);
				}
				for(int j=0; j<input[i]->nSampledLargeObs; j++){

					if(sum_logConst_forLarge != -DBL_MAX) sum_logConst_forLarge += lchoose(input[i]->largeObsSize[j], input[i]->largeObsPos[j]);

					largeObsPos[k]  = input[i]->largeObsPos[j];
					largeObsSize[k] = input[i]->largeObsSize[j];
					k++;
				}
			}else{
				for(int nSampled = 0; nSampled < nPoints; nSampled++){
					int r = nSampled + rand() % (input[i]->nSampledLargeObs - nSampled);

					if(sum_logConst_forLarge != -DBL_MAX) sum_logConst_forLarge += lchoose(input[i]->largeObsSize[r], input[i]->largeObsPos[r]);

					largeObsPos[k]  = input[i]->largeObsPos[r];
					largeObsSize[k] = input[i]->largeObsSize[r];
					if(r != nSampled){
						swap(input[i]->largeObsPos[r], input[i]->largeObsPos[nSampled]);
						swap(input[i]->largeObsSize[r], input[i]->largeObsSize[nSampled]);
					}
					k++;
				}
			}
		}
		if(k < nSampledLargeObs){
			printf("Warning: k = %d, nSampledLargeObs = %d\n", k, nSampledLargeObs);
			nSampledLargeObs = k;
		}
	}

	delete[] startIndex;
	delete[] pick;

	/*
	int *x = new int[total];
	int *n = new int[total];
	double *w = new double[total];
	k = 0;
	for(int i=0; i<num; i++){
		for(int j=0; j<input[i]->nSampledLargeObs; j++){
			x[k] = input[i]->largeObsPos[j];
			n[k] = input[i]->largeObsSize[j];
			w[k] = ((double)input[i]->nLargeObs) / input[i]->nSampledLargeObs;
			k++;
		}
	}
	delete[] x;
	delete[] n;
	delete[] w;
	*/
}

void BetaBinStats::finalize(void){
	if(ready) STOP_HERE("it's already finalized");

	// print(stdout);

	if(k < nSampledLargeObs) STOP3("k < nSampledLargeObs: k=%d, nSampledLargeObs=%d", k, nSampledLargeObs);
	if(obsCount_length > 0 && obsCount[obsCount_length-1] <= 0) STOP_HERE("obsCount[obsCount_length-1] <= 0");
	if(obsCountPos_length > 0 && obsCountPos[obsCountPos_length-1] <= 0) STOP_HERE("obsCount[obsCount_length-1] <= 0");
	if(obsCountNeg_length > 0 && obsCountNeg[obsCountNeg_length-1] <= 0) STOP_HERE("obsCount[obsCount_length-1] <= 0");
	if(nLargeObs > nSamplesForLarge && nSampledLargeObs != nSamplesForLarge) STOP_HERE("nLargeObs > nSamplesForLarge && nSampledLargeObs != nSamplesForLarge");
	if(nLargeObs < nSamplesForLarge && nSampledLargeObs != nLargeObs) STOP_HERE("nLargeObs < nSamplesForLarge && nSampledLargeObs != nLargeObs");

	for(int j=obsCount_length-2; j>=0; j--) obsCount[j] += obsCount[j+1];
	for(int j=obsCountPos_length-2; j>=0; j--) obsCountPos[j] += obsCountPos[j+1];
	for(int j=obsCountNeg_length-2; j>=0; j--) obsCountNeg[j] += obsCountNeg[j+1];

	ready = true;
}


BetaBinStats::~BetaBinStats(void){
	if(obsCount != NULL) delete[] obsCount;
	if(obsCountPos != NULL) delete[] obsCountPos;
	if(obsCountNeg != NULL) delete[] obsCountNeg;
	if(largeObsPos != NULL) delete[] largeObsPos;
	if(largeObsSize != NULL) delete[] largeObsSize;
}

void BetaBinStats::print(FILE * out){
	fprintf(out, "  BetaBinStats:\n");
	fprintf(out, "     nSmallObs=%d, nLargeObs=%d, threshold=%d, nSampledLargeObs=%d\n",nSmallObs,nLargeObs,threshold,nSampledLargeObs);
	fprintf(out, "     mu0=%f, var_mu=%f, gamma0=%f, var_gamma=%f\n",mu0,var_mu,gamma0,var_gamma);
	fprintf(out, "     --------------------------------------------------------------\n");
	fprintf(out, "     Count\tPos\tNeg\n");
	for(int i=0; i<obsCount_length; i++){
		int p=0, n=0;
		if(i<obsCountPos_length) p = obsCountPos[i];
		if(i<obsCountNeg_length) n = obsCountNeg[i];
		fprintf(out, "     %d\t%d\t%d\n", obsCount[i], p, n);
	}
	fprintf(out, "     --------------------------------------------------------------\n");
	fprintf(out, "     Size\tPos\n");
	for(int i=0; i<nSampledLargeObs; i++){
		fprintf(out, "     %d\t%d\n", largeObsSize[i], largeObsPos[i]);
	}
	fprintf(out, "  End\n");
}

double BetaBinStats::logLikelihood(const double *param, const bool ignoreConst, const bool avgloglik){
	double ans = logLik_BetaBinomial_mean_size(this, param, false, ignoreConst);
	ans /= (nSmallObs + nSampledLargeObs); // TODO: weigh the two parts correctly
	if(avgloglik) return ans;
	return ans * (nSmallObs + nLargeObs);
}

#define MIN_MEAN 1e-10
#define MAX_MEAN 1-1e-10
#define MIN_SIZE 1e-10
#define MAX_SIZE 1e10

/**
 * Return 0 if successful
 *        1 if det(H) = 0
 *        2 if it reaches maxIter
 * 1. Try initial param with the specified maxIter1 and nLnsrchStep1
 * 2. Try (mu0, gamma0) with the specified maxIter1 and nLnsrchStep1
 * 3. Try initial param with the specified maxIter2 and nLnsrchStep2
 * 4. Try (mu0, gamma0) with the specified maxIter2 and nLnsrchStep2
 */
int BetaBinStats::getPostMeanSize(
	double *param,  // IN/OUT: {mean, size}
					// IN: initial value (0 means using mu0 and gamma0)
	int *nIter,     // OUT: number of iterations taken
	double minMean, double maxMean,
	double minSize, double maxSize,
	double epsilon, // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
					//             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
	double stepSize,// IN: max Newton step size
	int maxIter1, int nLnsrchStep1,
	int maxIter2, int nLnsrchStep2, int verbose, int debug
){
	if(!ready) STOP_HERE("Not ready!!");

	int code, niter, code1=1, code2=1;
	double init_param[2];

	(*nIter) = 0;
	init_param[0] = param[0];  init_param[1] = param[1];

	double extended_minMean = MIN(MIN_MEAN,minMean),
		   extended_maxMean = MAX(MAX_MEAN,maxMean),
		   extended_minSize = MIN(MIN_SIZE,minSize),
		   extended_maxSize = MAX(MAX_SIZE,maxSize);

	// (1) Using init param, maxIter1, nLnsrchStep1, extended range
	code1 = 1;
	if(init_param[0] != 0 || init_param[1] != 0){
		param[0] = init_param[0]; param[1] = init_param[1];
		code1 = postMode_BetaBinomial_mean_size(
			param, &niter, this, epsilon, stepSize, maxIter1, nLnsrchStep1,
			extended_minMean, extended_maxMean, extended_minSize, extended_maxSize, verbose, debug
		);

		if(nIter) (*nIter) += niter;
		if(code1 == 0 && IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
	}

	// (2) Using mu0, gamma0, maxIter1, nLnsrchStep1, extended range
	param[0] = 0; param[1] = 0;
	code2 = postMode_BetaBinomial_mean_size(
		param, &niter, this, epsilon, stepSize, maxIter1, nLnsrchStep1,
		extended_minMean, extended_maxMean, extended_minSize, extended_maxSize, verbose, debug
	);
	if(nIter) (*nIter) += niter;
	if(code2 == 0 && IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;

	// (3) Using init param, maxIter1, nLnsrchStep2, extended range
	code = 1;
	if(code1 != 0 && code2 != 0 && (init_param[0] != 0 || init_param[1] != 0)){
		param[0] = init_param[0]; param[1] = init_param[1];
		code = postMode_BetaBinomial_mean_size(
			param, &niter, this, epsilon, stepSize, maxIter1, nLnsrchStep2,
			extended_minMean, extended_maxMean, extended_minSize, extended_maxSize, verbose, debug
		);
		if(nIter) (*nIter) += niter;
		if(code == 0 && IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
	}

	// (4) Using mu0, gamma0, maxIter1, nLnsrchStep2, extended range
	if(code1 != 0 && code2 != 0 && code != 0){
		param[0] = 0; param[1] = 0;
		code = postMode_BetaBinomial_mean_size(
			param, &niter, this, epsilon, stepSize, maxIter1, nLnsrchStep2,
			extended_minMean, extended_maxMean, extended_minSize, extended_maxSize, verbose, debug
		);
		if(nIter) (*nIter) += niter;
		if(code == 0 && IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
	}

	double best_param[2] = {-1, -1};
	double best_loglik = -DBL_MAX;

	// (5) Using init param, maxIter2, nLnsrchStep2, limited range
	if(init_param[0] != 0 || init_param[1] != 0){
		param[0] = init_param[0]; param[1] = init_param[1];
		code = postMode_BetaBinomial_mean_size(
			param, &niter, this, epsilon, stepSize, maxIter2, nLnsrchStep2,
			minMean, maxMean, minSize, maxSize, verbose, debug
		);
		if(nIter) (*nIter) += niter;
		if(code == 0){
			if(IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
			best_loglik = logLik_BetaBinomial_mean_size(this,param,true);
			best_param[0] = param[0]; best_param[1] = param[1];
		}
	}
	// (6) Using mu0, gamma0, maxIter2, nLnsrchStep2, limited range
	param[0] = 0; param[1] = 0;
	code = postMode_BetaBinomial_mean_size(
		param, &niter, this, epsilon, stepSize, maxIter2, nLnsrchStep2,
		minMean, maxMean, minSize, maxSize, verbose, debug
	);
	if(nIter) (*nIter) += niter;
	if(code == 0){
		if(IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
		double loglik = logLik_BetaBinomial_mean_size(this,param,true);
		if(loglik > best_loglik){
			best_loglik = loglik;
			best_param[0] = param[0]; best_param[1] = param[1];
		}
	}

	if(best_loglik > -DBL_MAX){
		param[0] = best_param[0];  param[1] = best_param[1];
		return 0;
	}

	/* backup
	// Using initial parameter setting with the specified range constraints
	param[0] = init_param[0]; param[1] = init_param[1];
	code = postMode_BetaBinomial_mean_size(
		param, &niter, this, epsilon, stepSize,
		maxIter1, nLnsrchStep1, minMean,maxMean,minSize,maxSize, verbose, debug
	);
	if(nIter) (*nIter) += niter;
	if(code == 0){
		if(IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
		best_loglik = logLik_BetaBinomial_mean_size(this,param,true);
		best_param[0] = param[0]; best_param[1] = param[1];
	}

	param[0] = 0; param[1] = 0;
	code = postMode_BetaBinomial_mean_size(
		param, &niter, this, epsilon, stepSize,
		maxIter1, nLnsrchStep1, minMean,maxMean,minSize,maxSize, verbose, debug
	);
	if(nIter) (*nIter) += niter;
	if(code == 0){
		if(IS_PARAM_IN_RANGE(param,minMean,maxMean,minSize,maxSize)) return 0;
		double loglik = logLik_BetaBinomial_mean_size(this,param,true);
		if(loglik > best_loglik){
			best_loglik = loglik;
			best_param[0] = param[0]; best_param[1] = param[1];
		}
	}
	if(best_loglik > -DBL_MAX){
		param[0] = best_param[0];  param[1] = best_param[1];
		return 0;
	}
	*/
	return code;
}


BetaBinMEStats::BetaBinMEStats(const int *pos, const int *size, const double *weight, const int length){
	num = 0; totalPostives = 0; totalTrials = 0;
	total = 0; sum_p = 0; sum_p2 = 0; sum_w_by_n = 0; sum_w2_by_n = 0; sum_w2 = 0;
	for(int i=0; i<length; i++){
		if(size[i] < 0) STOP_HERE("size[i] < 0");
		if(pos[i] > size[i]) STOP_HERE("pos[i] > size[i]");
		if(size[i] == 0) continue;
		double w_i = (weight == NULL ? 1 : weight[i]);
		double p_i = ((double)pos[i]) / ((double)size[i]);
		num ++;
		totalPostives += pos[i];
		totalTrials   += size[i];
		total  += w_i;
		sum_w2 += w_i * w_i;
		sum_p  += w_i * p_i;
		sum_p2 += w_i * p_i * p_i;
		sum_w_by_n  += w_i / size[i];
		sum_w2_by_n += w_i * w_i / size[i];
	}
}

BetaBinMEStats::BetaBinMEStats(const BetaBinMEStats *stats, const int length){
	num = 0; totalPostives = 0; totalTrials = 0;
	total = 0; sum_p = 0; sum_p2 = 0; sum_w_by_n = 0; sum_w2_by_n = 0; sum_w2 = 0;
	for(int i=0; i<length; i++){
		num    += stats[i].num;
		totalPostives += stats[i].totalPostives;
		totalTrials   += stats[i].totalTrials;
		total  += stats[i].total;
		sum_w2 += stats[i].sum_w2;
		sum_p  += stats[i].sum_p;
		sum_p2 += stats[i].sum_p2;
		sum_w_by_n  += stats[i].sum_w_by_n;
		sum_w2_by_n += stats[i].sum_w2_by_n;
	}
}

BetaBinMEStats::BetaBinMEStats(void){
	num = 0; totalPostives = 0; totalTrials = 0;
	total = 0; sum_p = 0; sum_p2 = 0; sum_w_by_n = 0; sum_w2_by_n = 0; sum_w2 = 0;
}

void BetaBinMEStats::add(int nPositives, int nTrials, double weight){
	if(nTrials < 0) STOP_HERE("nTrials < 0");
	if(nPositives > nTrials) STOP_HERE("nPositives > nTrials");
	if(nTrials == 0) return;
	double w_i = weight;
	double p_i = ((double)nPositives) / ((double)nTrials);
	num ++;
	totalPostives += nPositives;
	totalTrials   += nTrials;
	total  += w_i;
	sum_w2 += w_i * w_i;
	sum_p  += w_i * p_i;
	sum_p2 += w_i * p_i * p_i;
	sum_w_by_n  += w_i / nTrials;
	sum_w2_by_n += w_i * w_i / nTrials;
}

void BetaBinMEStats::add(const BetaBinMEStats *stats){
	num += stats->num;
	totalPostives += stats->totalPostives;
	totalTrials   += stats->totalTrials;
	total  += stats->total;
	sum_w2 += stats->sum_w2;
	sum_p  += stats->sum_p;
	sum_p2 += stats->sum_p2;
	sum_w_by_n  += stats->sum_w_by_n;
	sum_w2_by_n += stats->sum_w2_by_n;
}

double BetaBinMEStats::getMean(void){
	return sum_p / total;
}

double BetaBinMEStats::getSize(bool naive){
	double p = sum_p / total;
	double pq = p * (1 - p);
	double S = sum_p2 - total * p * p;
	if(naive){
		double var = S/(num-1);
		if(var >= pq) return 1e-20;
		if(var == 0)  return 1e20;
		return pq/var - 1;
	}
	S = (num-1) * S / num;
	double r = (S - pq*(sum_w_by_n - sum_w2_by_n/total)) / (pq*(total - sum_w2/total - sum_w_by_n + sum_w2_by_n/total));
	if(r < 1e-20) r = 1e-20;
	if(r > 1 - 1e-20) r = 1 - 1e-20;
	return 1/r - 1;
}


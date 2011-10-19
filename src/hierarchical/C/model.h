/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#ifndef HIERARCHICAL_MODEL_H_
#define HIERARCHICAL_MODEL_H_

#include "stdio.h"
#include "util.h"
#include "dag.hpp"
#include "cube.hpp"

#define VAR_BETABIN_MEAN_SIZE(mean, size) (mean*(1.0-mean)/(size+1.0))
#define IS_PARAM_IN_RANGE(param, minMean, maxMean, minSize, maxSize) (param[0] > minMean && param[0] < maxMean && param[1] > minSize && param[1] < maxSize)

/**
 * x_i ~ Binomial(p_i, size=n_i)
 * p_i ~ Beta(mu, gamma)
 * logit(mu)  ~ N(logit(mu0),  var_mu)
 * log(gamma) ~ N(log(gamma0), var_gamma)
 */
class BetaBinStats{

public:
	bool ready;

	int nSmallObs;
	int threshold;   // Small vs. Large

	int* obsCount;    // obsCount[j] = #obs with n_i - 1       >= j and n_i <= threshold
	int* obsCountPos; // obsCount[j] = #obs with x_i - 1       >= j and n_i <= threshold
	int* obsCountNeg; // obsCount[j] = #obs with n_i - x_i - 1 >= j and n_i <= threshold

	int obsCount_length;
	int obsCountPos_length;
	int obsCountNeg_length;

	int nLargeObs;
	int nSampledLargeObs;
	int *largeObsSize; // largeObsSize[i] = n_i; array of length nSampledLargeObs
	int *largeObsPos;  // largeObsPos[i]  = x_i; array of length nSampledLargeObs

	double sum_logConst_forSmall; // -DBL_MAX means ignore the constants
	double sum_logConst_forLarge; // -DBL_MAX means ignore the constants

	double mu0, var_mu, gamma0, var_gamma; // prior
	// NOTE: var_mu = var(logit(mu));  var_gamma = var(log(gamma))

	int nSamplesForLarge;

	BetaBinStats(
		const int *nPositives, const int *nTrials, const int length, const int threshold, const int nSamplesForLarge,
		const double mu0, const double var_mu, const double gamma0, const double var_gamma
	);
	BetaBinStats(
		const int threshold, const int nSamplesForLarge,
		const int nLargeObs, const int obsCount_length, const int obsCountPos_length, const int obsCountNeg_length,
		const double var_mu, const double var_gamma, const bool ignoreConst=true
	);
	void add(int nPositives, int nTrials);
	void set(BetaBinStats *input[], int num);
	void finalize(void);

	~BetaBinStats(void);

	// return log likelihood
	double logLikelihood(const double *param, const bool ignoreConst, const bool avgloglik); // param = {mean, size}

	/**
	 * Return 0 if successful
	 *        1 if det(H) = 0
	 *        2 if it reaches maxIter
	 *        3 if line search fails
	 * 1. Try initial param with the specified maxIter1 and nLnsrchStep1
	 * 2. Try (mu0, gamma0) with the specified maxIter1 and nLnsrchStep1
	 * 3. Try initial param with the specified maxIter2 and nLnsrchStep2
	 * 4. Try (mu0, gamma0) with the specified maxIter2 and nLnsrchStep2
	 */
	int getPostMeanSize(
		double *param,  // IN/OUT: {mean, size}
						// IN: initial value (0 means using mu0 and gamma0)
		int *nIter,     // OUT: number of iterations taken (nIter==NULL if don't care)
		double minMean=1e-10, double maxMean=1-1e-10,
		double minSize=1e-10, double maxSize=1e10,
		double epsilon=1e-10, // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
							  //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
		double stepSize=1,    // IN: max Newton step size
		int maxIter1=20,  int nLnsrchStep1=0,
		int maxIter2=100, int nLnsrchStep2=10, int verbose=0, int debug=0
	);

	void print(FILE* out);

private:
	int k;
};

/**
 * Moment estimator for the Beta-Binomial model
 */
class BetaBinMEStats;
class BetaBinMEStats {
public:
	double num;           // number of observations
	double totalPostives; // sum_i n_i p_i
	double totalTrials;   // sum_i n_i

	double total; 		// w = sum_i w_i
	double sum_p; 		// sum_i w_i p_i
	double sum_p2;		// sum_i w_i p_i^2
	double sum_w_by_n;	// sum_i w_i / n_i
	double sum_w2_by_n; // sum_i w_i^2 / n_i
	double sum_w2;      // sum_i w_i^2

	// weight can be null
	BetaBinMEStats(const int *nPositives, const int *nTrials, const double *weight, const int length);
	BetaBinMEStats(const BetaBinMEStats *stats, const int length);

	BetaBinMEStats(void);
	void add(int nPositives, int nTrials, double weight=1);
	void add(const BetaBinMEStats *stats);

	double getMean(void);
	double getSize(bool naive=false);
};

/**
 * return Pr(x <= nPositives), where x ~ BetaBinomial(mean, size, nTrials)
 * (*prevProb) = Pr(x-1 <= nPositive)
 */
double pbetabinom_MeanSize(int nPositives, int nTrials, double mean, double size, double *prevProb=NULL);

/**
 *                  [0]                    [1]                    [2]                    [3]
 * param = {    mean=mu,            size=gamma}
 * gr    = {(d LL/d mu),         (d LL/d gamma)}
 * H     = {(d^2 LL/d mu^2), (d^2 LL/d mu d gamma), (d^2 LL/d mu d gamma), (d^2 LL/d gamma^2)}
 */
double logLik_BetaBinomial_mean_size(const BetaBinStats* stats, const double *param, bool usePrior, bool ignoreConst=true);
void d_LogLik_BetaBinomial_mean_size(double* gr, const BetaBinStats* stats, const double *param, int verbose);
void d2_LogLik_BetaBinomial_mean_size(double* H, const BetaBinStats* stats, const double *param, int verbose);

/**
 * Return 0 if successful
 *        1 if det(H) = 0
 *        2 if it reaches maxIter
 *        3 if line search fails
 * Stop when |mean_current - mean_previous| <= epsilon and
 *           |sd_current/mean_current - sd_prev/mean_current| <= epsilon
 * Set stepSize = 1 (the regular Newton's method) if you do not know how to set it.
 *
 */
int postMode_BetaBinomial_mean_size(
	double *param, int *nIter, const BetaBinStats *stats, double epsilon, double stepSize,
	int maxIter, int nLnsrchSteps=10,
	double minMean=1e-10, double maxMean=1-1e-10, double minSize=1e-10, double maxSize=1e10,
	int verbose=0, int debug=0
);
int postMode_BetaBinomial_mean_size_old(
	double *param, int *nIter, const BetaBinStats *stats, double epsilon=1e-10, double stepSize=1,
	int maxIter=100, bool linesearch=true, int verbose=0, int debug=0
);

enum BetaBinomialStatus {
	PM_FAILED=-3, ME_TEMP=-2, INVALID=-1, LOW_SUPPORT=0, ME=1, PM=2
};

class BetaBinomialDAG;
class BetaBinomialCube;

class BetaBinPriorSetter {
public:
	virtual ~BetaBinPriorSetter(void){}
	/**
	 * This function set the prior stats->mu0 and stats->gamma0
	 * return status code: 0 -> success
	 */
	virtual int set(BetaBinomialDAG *dag, int pnode, int tnode, BetaBinStats *stats, int verbose){ STOP_HERE("Please implement a BetaBinPriorSetter"); return -1; }
	virtual int set(BetaBinomialCube *cube, int pnode, int tvector[], BetaBinStats *stats, int verbose){ STOP_HERE("Please implement a BetaBinPriorSetter"); return -1; }
};

class BetaBinLikelihood_PriorSetter : public BetaBinPriorSetter {
	bool ignoreOutOfBound;
public:
	BetaBinLikelihood_PriorSetter(bool ignoreOutOfBound){ this->ignoreOutOfBound = ignoreOutOfBound; }
	int set(BetaBinomialDAG *dag, int pnode, int tnode, BetaBinStats *stats, int verbose);
	int set(BetaBinomialCube *cube, int pnode, int tvector[], BetaBinStats *stats, int verbose);
};

class Fixed_PriorSetter : public BetaBinPriorSetter {
public:
	int set(BetaBinomialDAG *dag, int pnode, int tnode, BetaBinStats *stats, int verbose);
	int set(BetaBinomialCube *cube, int pnode, int tvector[], BetaBinStats *stats, int verbose);
};

/**
 *	NODE ID ORDERING:
 *	  * Root should have ID = 0.
 *	  * Population nodes and treatment nodes are labeled separately; i.e., both start from 0.
 *	  * For any given node, every one of its ancestors must have ID < the ID of the given node.
 *
 *	INPUT OBSERVATION DATA ORDERING:
 *	  * Observation data table: data_(pNodeID, tID, nPositives, nTrials)
 *	  * This table must be ordered by pNodeID
 */
class BetaBinomialDAG {

public:

	bool debug;
	int verbose;

	SimpleDagIndex *pDagIndex; // population DAG index
	SimpleDagIndex *tDagIndex; // treatment DAG index
	SimpleDagIndex *tMapping;  // treatment mapping (from ID to node)

	int data_length;
	int *data_pNodeID;    // population node ID    	(length: data_length)
	int *data_tID;        // treatmentID			(length: data_length)
	int *data_nPositives; // number of positives	(length: data_length)
	int *data_nTrials;    // number of trials		(length: data_length)

	int *data_startIndex; // data_pNodeID[data_startIndex[pnode], ..., data_startIndex[pnode+1]-1] = pnode

	int nPNodes; // number of population nodes
	int nTNodes; // number of treatment nodes

	int nTrialsForSmall;   // the threshold for the statistics for posterior mode estimation
	int nObsForLarge;      // number of observations in the reservoir sample for posterior mode estimation

	double root_priorMean;
	double root_priorSize;
	double var_logitMean;
	double var_logSize;

	double minMean, maxMean, minSize, maxSize;

	int ME_threshold_nObs;		 // Use moment estimates for notes with
	double ME_threshold_nTrials; // #observations >= ME_threshold_nObs and avg(#trial) >= ME_threshold_nTrials

	int minNumObs;			// Do not estimate the parameters for a node if the node has too little support
	int minTotalTrials;		//   #Treatments < minNumObs OR Total # of Trials    < minTotalTrials
	int minTotalPositives;  //                           OR Total # of Positives < minTotalPositives

	virtual inline double getMean(int pnode, int tnode){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); return param[C_3DA(pnode,tnode,0,nPNodes,nTNodes)]; }
	virtual inline double getSize(int pnode, int tnode){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); return param[C_3DA(pnode,tnode,1,nPNodes,nTNodes)]; }
	virtual inline int getStatus(int pnode, int tnode){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); return status[C_MAT(pnode,tnode,nPNodes)]; }
	virtual inline void setMean(int pnode, int tnode, double value){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); param[C_3DA(pnode,tnode,0,nPNodes,nTNodes)] = value; }
	virtual inline void setSize(int pnode, int tnode, double value){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); param[C_3DA(pnode,tnode,1,nPNodes,nTNodes)] = value; }
	virtual inline void setStatus(int pnode, int tnode, BetaBinomialStatus value){ CHK_C_INDEX(pnode, nPNodes); CHK_C_INDEX(tnode, nTNodes); status[C_MAT(pnode,tnode,nPNodes)] = value; }

	/**
	 * IMPORTANT: param and status must have been allocated
	 */
	BetaBinomialDAG(
		double *param, // (output)  #pNodes * #tNodes * 2
		int *status,   // (output)  #pNodes * #tNodes
		// Observation data table: data_(pnode_id, t_id, nPositive, nTrial)
		const int *data_pnode_id, const int *data_t_id, const int *data_nPositive, const int *data_nTrial, int data_len,
		// pDAG: Graph of population nodes
		const int *pDAG_child_id, const int *pDAG_parent_id, int num_pNodes, int num_pDAG_edges,
		// tDAG: Graph of treatment nodes
		const int *tDAG_child_id, const int *tDAG_parent_id, int num_tNodes, int num_tDAG_edges,
		// tMap: Mapping from treatment IDs to treatment leaf nodes
		const int *tMap_t_id, const int *tMap_tnode_id, int num_t_ids, int num_tMap_edges,
		// options
		double minMean=1e-10, double maxMean=1-1e-10,
		double minSize=1e-10, double maxSize=1e10,
		int minNumObs=3, int minTotalTrials=10, int minTotalPositives=0,
		int verbose=0, bool debug=false, bool resetParamAndStatus=true
	);

	virtual ~BetaBinomialDAG(void);

	/**
	 * if(naive): estimate parameters based on mean and var of nPositive/nTrials
	 */
	void computeMomentEstimates(int nRefine=0, bool naive=false);

	void computePosteriorMode(
		double root_priorMean, double root_priorSize,
		double varLogitMean, double varLogSize, int ME_threshold_numObs, double ME_threshold_numTrials,
		int nTrials_ForSmall, int nObs_ForLarge, BetaBinPriorSetter *priorSetter,
		double epsilon=1e-10, // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
							  //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
		double stepSize=1,    // IN: max Newton step size
		int maxIter1=20,  int nLnsrchStep1=0, int maxIter2=100, int nLnsrchStep2=10
	);

	// return test-set log-likelihood
	// 	  averageLoglik=true: average log likelihood per observation
	//   averageLoglik=false: sum of log likelihood
	double testLogLikelihood(
		double *node_loglik, // (output) #pNodes * #tNodes: LogLikelihood per observation of each node
		                     //          set node_loglik=NULL to disable the output
							 //          0 means no data in the node
		const int *data_pnode_id, const int *data_t_id, const int *data_nPositive, const int *data_nTrial, int data_len,
		int nTrials_ForSmall, int nObs_ForLarge,
		int option, // 0: leaves only,  1: the whole DAG
		bool averageLoglik
	);

	double confIntervalCoverage(
		double *nObsInCI,  // (output) #pNodes x #tNodes (NULL means disable)
		int *nObs,         // (output) #pNodes x #tNodes (NULL means disable; 0 means no output for the node)
		double prob, // probability: (prob*100)% confidence interval will be constructed
		const int *select, // #pNodes x #tNodes (NULL means select all; non-zero means selected)
		const int *data_pnode_id, const int *data_t_id, const int *data_nPositive, const int *data_nTrial, int data_len,
		int *nFailed
	);

private:
	double *param; // nPNodes x nTNodes x 2
	int *status;   // nPNodes x nTNodes
};

typedef struct {
	int nLargeObs;
	int obsCount_length;
	int obsCountPos_length;
	int obsCountNeg_length;
	int nObs;
	int totalTrials;
} BetaBinSizeInfo;

/**
 *  TODO: Validate BetaBinomialCube
 *  TODO: handle minMean, maxMean, minSize, maxSize
 *  IMPORTANT NOTE: BetaBinomialCube has not yet been validated at all!! It may just crash
 *
 *	NODE ID ORDERING:
 *	  * Root of each hierarchy should have ID = 0.
 *	  * Population node: DAG
 *	  * Treatment entry: vector of nodes; each dimension has a tree structure
 *	  * For any given node, every one of its ancestors must have ID < the ID of the given node.
 *
 *	INPUT OBSERVATION DATA ORDERING:
 *	  * Observation data table: data_(pNodeID, tNodeID[0], ..., tNodeID[nTDim-1], nPositives, nTrials)
 *	  * This table must be ordered by pNodeID
 */
class BetaBinomialCube : public IBottomUpCubeComputable {

public:

	bool debug;
	int verbose;

	SimpleDagIndex *pDagIndex; // population DAG index
	BottomUpCubeProcessor *cubeProcessor;

	int data_length;
	int *data_pNodeID;    // population node ID    	(length: data_length)
	R_2DArray<int> *data_tNodeID; // data_length x nTDim
	int *data_nPositives; // number of positives	(length: data_length)
	int *data_nTrials;    // number of trials		(length: data_length)

	int *data_startIndex; // data_pNodeID[data_startIndex[pnode], ..., data_startIndex[pnode+1]-1] = pnode

	int nPNodes; // number of population nodes
	int nTDim;   // number of treatment dimensions

	int nTrialsForSmall;   // the threshold for the statistics for posterior mode estimation
	int nObsForLarge;      // number of observations in the reservoir sample for posterior mode estimation

	double root_priorMean;
	double root_priorSize;
	double var_logitMean;
	double var_logSize;
	double minMean, maxMean, minSize, maxSize;

	int ME_threshold_nObs;		 // Use moment estimates for notes with
	double ME_threshold_nTrials; // #observations >= ME_threshold_nObs and avg(#trial) >= ME_threshold_nTrials

	virtual inline double getMean(int pnode, int tnode[]){ return paramMean->get(getEntry(pnode,tnode)); }
	virtual inline double getSize(int pnode, int tnode[]){ return paramSize->get(getEntry(pnode,tnode)); }
	virtual inline int getStatus(int pnode, int tnode[]){  return status->get(getEntry(pnode,tnode)); }
	virtual inline void setMean(int pnode, int tnode[], double value){ paramMean->set(getEntry(pnode,tnode), value); }
	virtual inline void setSize(int pnode, int tnode[], double value){ paramSize->set(getEntry(pnode,tnode), value); }
	virtual inline void setStatus(int pnode, int tnode[], BetaBinomialStatus value){ status->set(getEntry(pnode,tnode), value); }

	/**
	 * IMPORTANT: paramMean, paramSize and status must have been allocated
	 */
	BetaBinomialCube(
		double *paramMean, // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
		double *paramSize, // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
		int *status,       // (output)  nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
		const int *data_pnode_id, // data_len x 1
		const int *data_tnode_id, // data_len x nTDim
		const int *data_nPositive,
		const int *data_nTrial,
		int data_len, int nTDim,
		const int *pDAG_child_id, const int *pDAG_parent_id, int num_pNodes, int num_pDAG_edges,
		DimSchemaHierarchical *tCubeDim[], // nTDim hierarchies for the treatments
		double minMean=1e-10, double maxMean=1-1e-10,
		double minSize=1e-10, double maxSize=1e10,
		int verbose=0, bool debug=false
	);

	virtual ~BetaBinomialCube(void);

	void computeMomentEstimates(void);

	void computePosteriorMode(
		double root_priorMean, double root_priorSize,
		double varLogitMean, double varLogSize, int ME_threshold_numObs, double ME_threshold_numTrials,
		int nTrials_ForSmall, int nObs_ForLarge, BetaBinPriorSetter *priorSetter,
		double epsilon=1e-10, // IN: stop if |mu(t) - mu(t-1)| <= epsilon and
							  //             |cv(t) - cv(t-1)| <= epsilon (cv: coeff of variation)
		double stepSize=1,    // IN: max Newton step size
		int maxIter1=20,  int nLnsrchStep1=0, int maxIter2=100, int nLnsrchStep2=10
	);

    void initialize(BottomUpCubeProcessor *processor){
    	if(action == ME) initialize_ME(processor); else if(action == PM) initialize_PM(processor); else STOP_HERE("unknown action");
    }
    void processOneBaseSubset(int *tnode, int ndim, BottomUpCubeProcessor *processor){
    	if(action == ME) processOneBaseSubset_ME(tnode, ndim, processor); else if(action == PM) processOneBaseSubset_PM(tnode, ndim, processor); else STOP_HERE("unknown action");
    }
    void processOneNonbaseSubset(int *tnode, int ndim, BottomUpCubeProcessor *processor){
    	if(action == ME) processOneNonbaseSubset_ME(tnode, ndim, processor); else if(action == PM) processOneNonbaseSubset_PM(tnode, ndim, processor); else STOP_HERE("unknown action");
    }
    void finalizeCube(BottomUpCubeProcessor *processor){
    	if(action == ME) finalizeCube_ME(processor); else if(action == PM) finalizeCube_PM(processor); else STOP_HERE("unknown action");
    }

private:
	R_MultiDimArray<double> *paramMean; // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	R_MultiDimArray<double> *paramSize; // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]
	R_MultiDimArray<int> *status;   // nPNodes x nTNodes[0] x ... x nTNodes[nTDim-1]

	int *temp; // TODO: Any function that access this variable is not safe for multi-threading

	inline const int* getEntry(int pnode, int tnode[]){
		CHK_C_INDEX(pnode, nPNodes);
		temp[0] = pnode; for(int i=0; i<nTDim; i++) temp[i+1] = tnode[i];
		return temp;
	}

	int pNodeID;           // pNodeID and action decide which cube to be processed
	BetaBinomialStatus action;
	R_MultiDimArray<BetaBinMEStats*> *MEStats;
	R_MultiDimArray<BetaBinStats*> *PMStats;
	R_MultiDimArray<BetaBinSizeInfo> *sizeInfo;

	void initialize_ME(BottomUpCubeProcessor *processor);
    void processOneBaseSubset_ME(int *tnode, int ndim, BottomUpCubeProcessor *processor);
    void processOneNonbaseSubset_ME(int *tnode, int ndim, BottomUpCubeProcessor *processor);
    void finalizeCube_ME(BottomUpCubeProcessor *processor);

    void initialize_PM(BottomUpCubeProcessor *processor);
    void processOneBaseSubset_PM(int *tnode, int ndim, BottomUpCubeProcessor *processor);
    void processOneNonbaseSubset_PM(int *tnode, int ndim, BottomUpCubeProcessor *processor);
    void finalizeCube_PM(BottomUpCubeProcessor *processor);
};


#endif /* HIERARCHICAL_MODEL_H_ */

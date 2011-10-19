/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#ifndef CUBE_HPP_
#define CUBE_HPP_

#include "dag.hpp"
#include <limits.h>

class LatticeEnumeration;
class BottomUpCubeProcessor;

enum DimSchemaHierarchical_LevelType { ROOT_AT_LEVEL_0, LEAVES_AT_LEVEL_0 };

/**
 * Dimension schema for a tree-structured hierarchy.
 *
 * Each value/node in the hierarchy is indexed by a value ID in [0, N-1], where
 * N is the total number values/nodes in the hierarchy.
 * NodeID=0 must be the root.
 *
 * If levelType == ROOT_LEVEL_0,
 * 		The root is at level 0; the children of the root are at level 1, and so on.
 * 		Leaves may be at different levels.
 *
 * If levelType == LEAVES_LEVEL_0,
 *		All the leaf nodes are at level 0.
 * 		The last level corresponds to the root.
 * 		The tree does not need to be balanced.
 *
 * @author beechung
 *
 */
class DimSchemaHierarchical {
public:

	DimSchemaHierarchical(
		const int *nodeID, const int *parentID, const int nNodes, const int nEdges,
		const DimSchemaHierarchical_LevelType levelType=LEAVES_AT_LEVEL_0 // Do NOT change the default!!
	);
	~DimSchemaHierarchical(void);

	inline int numLevels(void) const { return nLevels; }
	inline int numNodes(void) const { return nNodes; }
	inline int numNodesAtLevel(int i) const { if(i<0 || i>=nLevels) STOP3("i=%d, nLevels=%d", i, nLevels); return nNodesAtLevel[i]; }
	inline const int* getNodeIDsAtLevel(int i) const { if(i<0 || i>=nLevels) STOP3("i=%d, nLevels=%d", i, nLevels); return nodeIDsAtLevel[i]; }
	inline int numChildren(int nodeID) const { return childrenIndex->numNextNodes(nodeID); }
	inline const int* getChildren(int nodeID) const { return childrenIndex->getNextNodes(nodeID); }
	inline int getParent(int nodeID) const { return (nodeID == 0 ? -1 : parentsIndex->getNextNodes(nodeID)[0]); }
	inline bool isBase(int nodeID) const { if(nodeID < 0 || nodeID >= nNodes) STOP_HERE("index out of bound"); return (level[nodeID] == 0);}
	inline int getLevel(int nodeID) const { if(nodeID < 0 || nodeID >= nNodes) STOP_HERE("index out of bound"); return level[nodeID]; }

private:
	SimpleDagIndex *parentsIndex;
	SimpleDagIndex *childrenIndex;
	int nNodes;
	int **nodeIDsAtLevel; // nodeIDsAtLevel[i] is an array of node IDs at level i
	int *nNodesAtLevel;   // nNodeAtLevel[i] is the number of nodes at level i
	int nLevels;
	int *level; // level[i] is the level of node i
};

/**
 * Base class of the object to be passed into BottomUpCubeProcessor.
 *
 *    initialize(): Initialize data structures (before cube computation).
 *    processOneBaseSubset(int *dimVID, ...): Process the base (finest-grained) subset that is
 *    			identified by dimVID. dimVID[i] is the value ID (the sequence number of
 *    			a node in the hierarchy) in the i-th hierarchical dimension.
 *    processOneNonbaseSubset(int *dimVID, ...): Process the non-base subset identified by dimVID.
 *    finalizeCube(): Perform the final processing of the data structures.
 */
class IBottomUpCubeComputable {
public:
    virtual void initialize(BottomUpCubeProcessor *processor)=0;
    virtual void processOneBaseSubset(int *dimVID, int ndim, BottomUpCubeProcessor *processor)=0;
    virtual void processOneNonbaseSubset(int *dimVID, int ndim, BottomUpCubeProcessor *processor)=0;
    virtual void finalizeCube(BottomUpCubeProcessor *processor)=0;
    virtual ~IBottomUpCubeComputable(void){};
};

/**
 * Process the cube in a bottom-up manner.
 *
 * Call compute(obj) to perform the bottom-up cube computation, where obj implements
 * IBottomUpCubeProcessor.
 * It will do the following in order:
 * 		obj.initialize(this),
 * 		obj.processOneBaseSubset(..,this),
 * 		obj.processOneNonbaseSubset(..,this),
 * 		finalizeCube(this).
 *
 * Useful methods:
 *
 * getDimValues(int[] dimVID): Get the String values from VIDs.
 * getSubsetIndex(int[] dimVID): Convert a dimVID tuple to a single index number running
 * 				from 0 to NumSubsets - 1.
 * getDimVIDs(int subsetIndex): Convert a subsetIndex to a dimVID tuple.
 *
 * @author Bee-Chung Chen
 *
 */
class BottomUpCubeProcessor {
public:
	/** the size of each dimension. */
    int *dimSize;

    /** dimNumLevels[i]: the number of levels of dim i. */
    int *dimNumLevels;

    /** dimSizeByLevel[i][k]: the k-th level of dim i. */
    int **dimSizeByLevel;

    int verbose;

    /** The dimensions in schema used to create the cube.
     *  Note: Only hierarhical dimensions in schema will be used.*/
    DimSchemaHierarchical **CubeDimension;

    int NumSubsets;
    int NumBaseSubsets;

    LatticeEnumeration *MultiDimLevels; // e.g., MultiDimLevels.level[0] = {0,0} is the base level.

    int nDim;

    IBottomUpCubeComputable *cube;

    BottomUpCubeProcessor(DimSchemaHierarchical **cubeDims, int ndim);
    ~BottomUpCubeProcessor(void);

    inline bool isBase(const int *dimVID) const {
    	for(int i=0; i<nDim; i++)
    		if(!CubeDimension[i]->isBase(dimVID[i])) return false;
    	return true;
    }

    inline bool isBase(int subsetIndex) const {
    	getDimVIDs(subsetIndex, temp_dimVID);
    	return isBase(temp_dimVID);
    }

    inline int getSubsetIndex(const int *dimVID) const {
    	return toLinearIndex(dimVID,dimSize,nDim);
    }

	inline void getDimVIDs(int subsetIndex, int *dimVID) const {
		if(subsetIndex >= NumSubsets) STOP_HERE("subsetIndex >= NumSubsets");
		fillInDimIndex(subsetIndex, dimVID, dimSize, nDim);
	}

    inline int getDimWithFewestChildren(const int *dimVID) const {
        // dimVID[i]: the value id on the i-th dimension
        int result = -1, min_num = INT_MAX;
        for(int i=0; i<nDim; i++){
            int nChildren = CubeDimension[i]->numChildren(dimVID[i]);
            if(nChildren < min_num){ result=i; min_num=nChildren;}
        }
        return result;
    }

    /**
     * Return N: the number of children (-1 if error)
     * childIndex[size] is allocated outside this function
     */
    int getBestChildSubsetIndices(const int *dimVID, int *childIndex, int size) const;

	void compute(IBottomUpCubeComputable *cube);

private:

	// TODO: Any function that accesses this variable is not safe for multi-threading
    int *temp_dimVID;

    void processBaseSubsets();
    void processNonbaseSubsets();

    /** levelSize[i]: the size of level[i] on dimension i.
     *  level[i]: the current level of dimension i.
     */
    inline void fillInLevelSize(int* levelSize, const int* level) const {
        for(int i=0; i<nDim; i++) levelSize[i] = dimSizeByLevel[i][level[i]];
    }

    /** range[i][0, ..., N-1]: the node IDs on dimension i at level[i].
     */
    inline void fillInRange(int** range, const int* level) const {
        for(int i=0; i<nDim; i++){
            range[i] = const_cast<int*>(CubeDimension[i]->getNodeIDsAtLevel(level[i]));
        }
    }

    inline int getStopIndex(const int* size) const {
        int stop = size[0];
        for(int i=1; i<nDim; i++) stop *= size[i];
        return stop;
    }

    inline void computeDimIndex(int linearIndex, int *dim, const int *levelSize, int **range){
        fillInDimIndex(linearIndex, dim, levelSize, nDim);
        for(int i=0; i<nDim; i++){
            dim[i] = range[i][dim[i]];
        }
    }

};

/**
 * Enumerate the levels in a bottom-up manner.
 * E.g., if the input is {3,2}, then the level will be:
 *  {{0,0},
 *   {0,1},{1,0},
 *   {1,1},{2,0},
 *   {2,1}
 *  }
 * numOfLevels[i] is the number of levels of the i-th dimension.
 */
class LatticeEnumeration {
public:
	int **level; // length x nDim
	int length;
	int nDim;

	LatticeEnumeration(int *numOfLevels, int ndim, bool display=false);
	~LatticeEnumeration(void);
};

class DimHierarchyArray {
public:
	DimSchemaHierarchical** hierarchy; // hierarchy[i]: pointer to the ith hierarchy
	int nDim;

	/**
	 *  INPUT: (dim[i], node_id[i], parent_id[i]) is a table
	 *         specifying that on the dim[i]-th dimension, node_id[i]
	 *         has parent parent_id[i].
	 *         This table MUST be sorted by dim.
	 */
	DimHierarchyArray(
		const int *dim,       // size: length
		const int *node_id,	  // size: length
		const int *parent_id, // size: length
		int nDim,
		int length
	);
	~DimHierarchyArray(void);
};


#endif /* CUBE_HPP_ */

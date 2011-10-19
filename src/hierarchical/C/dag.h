/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#ifndef DAG_H
#define DAG_H

typedef struct {
	int nNodes;		 // number of nodes
	int* startIndex; // length: nNodes
	int* parentNode; // parentNode[startIndex[i], ..., startIndex[i+1]-1] are
				     //   the parentNode IDs for node i
					 // startIndex[i] == startIndex[i+1] means node i has no parent
					 // Node IDs start from 0, NOT 1
	int length_parentNode;
} DAG_ParentIndex;

DAG_ParentIndex* create_DAG_ParentIndex(const int* nodeID, const int* parentID, const int nNodes, const int arrayLength);
void free_DAG_ParentIndex(DAG_ParentIndex* dag);

typedef struct {
	int nNodes;
	char* nodeData; // length: nNodes
} UpwardTraversalStruct;

UpwardTraversalStruct* create_UpwardTraversalStruct(const int nNodes);
void free_UpwardTraversalStruct(UpwardTraversalStruct* uts);


/**
 * IMPORTANT NOTE: This function should NOT be used in multi-threading!!
 * This function gets all the internal node IDs that are the ancestors of
 * the specified leaf.
 *
 * return N: Number of ancestors
 * ancestor[0, ..., N-1] are the output ancestor node IDs
 *
 * Let leaf=k. The traversal starts from the kth leaf; i.e.,
 * leafNodes->startIndex[leaf]
 *
 * If the number of ancestors > maxNumOfAncestors, only the first
 * maxNumOfAncestors will be in the output. "First" means close to
 * the leaf node.
 */
int get_allAncestorsFromLeaf(
	// OUTPUT:
	int* ancestor, // length: maxNumOfAncestors
	// INPUT
	const int leaf,
	const int maxNumOfAncestors,
	const DAG_ParentIndex* leafNodes, const DAG_ParentIndex* internalNodes,
	UpwardTraversalStruct* uts
);

/**
 * IMPORTANT NOTE: This function should NOT be used in multi-threading!!
 * This function gets all the node IDs that are the ancestors of
 * the input nodes.
 *
 * return N: Number of ancestors including the input nodes
 * ancestor[0, ..., N-1] are the output ancestor node IDs
 *
 * If the number of ancestors > maxNumOfAncestors, only the first
 * maxNumOfAncestors will be in the output. "First" means close to
 * the input nodes.
 */
int get_allAncestors(
	// OUTPUT:
	int* ancestor, // length: maxNumOfAncestors
	// INPUT
	const int *inputNodeID,
	const int nInputNodes,
	const int maxNumOfAncestors,
	const DAG_ParentIndex* internalNodes,
	UpwardTraversalStruct* uts,
	int debug
);

#endif

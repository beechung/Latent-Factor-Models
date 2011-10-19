/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/


#ifndef DAG_HPP
#define DAG_HPP

extern "C" {
#include "dag.h"
}

#include "util.h"

class SimpleDagIndex {

	DAG_ParentIndex* index;      // Naming convention: node1 -> node2 iff node2 is the parent of node1
	UpwardTraversalStruct *uts;  // Upward means from a node to all its ancestors

public:
	/**
	 * (nodeID[i] -> nextNodeID[i]), for i = 0, ..., nEdges-1, are the edges
	 * nodeIDs are between 0 and nNodes-1.
	 */
	SimpleDagIndex(const int* nodeID, const int* nextNodeID, const int nNodes, const int nEdges){
		index = create_DAG_ParentIndex(nodeID, nextNodeID, nNodes, nEdges);
		uts = NULL;
	}
	~SimpleDagIndex(void){
		free_DAG_ParentIndex(index);
		if(uts != NULL) free_UpwardTraversalStruct(uts);
	}
	int numNextNodes(const int nodeID) const {
		if(nodeID < 0 || nodeID >= index->nNodes) STOP_HERE("nodeID < 0 || nodeID >= index->nNodes");
		int start = index->startIndex[nodeID];
		int end = (nodeID+1 < index->nNodes) ? index->startIndex[nodeID+1] : index->length_parentNode;
		return end-start;
	}
	/**
	 * nextNodes(node)[0], ..., nextNodes(node)[numNextNodes(node)-1] are the next nodes from the
	 * input node.
	 */
	const int* getNextNodes(const int nodeID) const {
		if(nodeID < 0 || nodeID >= index->nNodes) STOP_HERE("nodeID < 0 || nodeID >= index->nNodes");
		return &(index->parentNode[index->startIndex[nodeID]]);
	}
	/**
	 * IMPORTANT NOTE: This function should NOT be used in multi-threading!!
	 * This function gets all the node IDs reachable from the input nodes,
	 * INCLUDING the input nodes.
	 * reachableNodes[0, ..., nInputNodes-1] = input node IDs
	 *
	 * return N: Number of reachable nodes stored in reachableNodes.
	 *           It stops at maxNumOutputNodes.
	 */
	int getReachableNodes(
		// OUTPUT:
		int* reachableNodes, // length: maxNumOutputNodes
		// INPUT
		const int *inputNodeIDs, const int nInputNodes, const int maxNumOutputNodes, bool debug=false
	){
		if(uts == NULL) uts = create_UpwardTraversalStruct(index->nNodes);
		return get_allAncestors(reachableNodes, inputNodeIDs, nInputNodes, maxNumOutputNodes, index, uts, debug ? 1 : 0);
	}
	int getReachableNodes(
		// OUTPUT:
		int* reachableNodes, // length: maxNumOutputNodes
		// INPUT
		const int nodeID, const int maxNumOutputNodes, bool debug=false
	){
		int node[1];
		node[0] = nodeID;
		return getReachableNodes(reachableNodes, node, 1, maxNumOutputNodes, debug);
	}

};

#endif

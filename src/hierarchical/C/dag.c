/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "dag.h"

#define STOP(msg) error(__FILE__, __LINE__, msg)
#define STOP2(msg,x) error(__FILE__, __LINE__, msg, x)
#define STOP3(msg,x1,x2) error(__FILE__, __LINE__, msg, x1, x2)
#define STOP4(msg,x1,x2,x3) error(__FILE__, __LINE__, msg, x1, x2, x3)

static void error(const char *filename, int lineno, const char *fmt, ...){
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	fprintf(stderr, "ERROR in file %s at line %d: ", filename, lineno);
	fprintf(stderr, buf);
	fprintf(stderr, "\n");
    exit(-1);
}

DAG_ParentIndex* create_DAG_ParentIndex(const int* nodeID, const int* parentID, const int nNodes, const int length){
	int i;

	int* start  = (int*)malloc(sizeof(int) * nNodes);
	int* parent = (int*)malloc(sizeof(int) * length);

	int* ind = (int*)malloc(sizeof(int) * nNodes);
	int* num = (int*)malloc(sizeof(int) * nNodes);
	for(i=0; i<nNodes; i++) num[i] = 0;

    for(i=0; i<length; i++){
        int id = nodeID[i];
        if(id >= nNodes){ STOP3("node id out of bound: %d >= %d", id, nNodes);}
        num[id]++;
    }

    start[0] = 0; ind[0] = 0;
    for(i=1; i<nNodes; i++){
        start[i] = start[i-1]+num[i-1];
        ind[i] = start[i];
    }

    for(i=0; i<length; i++){
        int id = nodeID[i];
        parent[ind[id]] = parentID[i];
        ind[id]++;
    }

    DAG_ParentIndex* out = (DAG_ParentIndex*)malloc(sizeof(DAG_ParentIndex));
    out->nNodes = nNodes;
    out->startIndex = start;
    out->parentNode = parent;
    out->length_parentNode = length;

#ifdef SANITY_CHECK //-------------------------------------------------------
	for(i=0; i<nNodes; i++){
		if(ind[i] != start[i]+num[i]){
			STOP("sanity check failed");
		}
	}
	if(start[nNodes-1]+num[nNodes-1] != length){
		STOP("sanity check failed");
	}
#endif //---------------------------------------------------------------------
    free(num);
    free(ind);

    return out;
}

void free_DAG_ParentIndex(DAG_ParentIndex* dag){
	if(dag != NULL){
		if(dag->startIndex != NULL) free(dag->startIndex);
		if(dag->parentNode != NULL) free(dag->parentNode);
	}
	free(dag);
}

UpwardTraversalStruct* create_UpwardTraversalStruct(const int nNodes){
	UpwardTraversalStruct* uts = (UpwardTraversalStruct*)malloc(sizeof(UpwardTraversalStruct));
	int i;
	uts->nNodes = nNodes;
	uts->nodeData = (char*)malloc(sizeof(char)*nNodes);
	for(i=0; i<nNodes; i++){
		uts->nodeData[i] = 0;
	}
	return uts;
}

void free_UpwardTraversalStruct(UpwardTraversalStruct* uts){
	if(uts != NULL){
		if(uts->nodeData != NULL){
			int i;
			for(i=0; i<uts->nNodes; i++){
				if(uts->nodeData[i] != 0){
					STOP("invalid state");
				}
			}
			free(uts->nodeData);
		}
		free(uts);
	}
}

static int add_an_ancestor(
	int* ancestor, // length: maxNumOfAncestors
	const int nodeID,
	const int nAncestors,
	const int maxNumOfAncestors,
	UpwardTraversalStruct* uts
){
	if(maxNumOfAncestors == nAncestors) return nAncestors;
	else if(nAncestors > maxNumOfAncestors){
		STOP("nAncestors > maxNumOfAncestors");
	}
	if(nodeID >= uts->nNodes){ STOP("nodeID >= uts->nNodes"); }
	if(uts->nodeData[nodeID] == 1) return nAncestors;
	else if(uts->nodeData[nodeID] != 0){ STOP("invalid state");}
	ancestor[nAncestors] = nodeID;
	uts->nodeData[nodeID] = 1;
	return nAncestors+1;
}

int get_allAncestorsFromLeaf(
	// OUTPUT:
	int* ancestor, // length: maxNumOfAncestors
	// INPUT
	const int leaf,
	const int maxNumOfAncestors,
	const DAG_ParentIndex* leafNodes, const DAG_ParentIndex* internalNodes,
	UpwardTraversalStruct* uts
){
	int i;
	if(leaf >= leafNodes->nNodes){ STOP("Input leaf number is out of bound");}
	if(internalNodes->nNodes != uts->nNodes){ STOP("internalNodes->nNodes != uts->nNodes");}
	int nAncestors = 0;
	int start = leafNodes->startIndex[leaf];
	int end = (leaf+1 < leafNodes->nNodes) ? leafNodes->startIndex[leaf+1] : leafNodes->length_parentNode;
	for(i=start; i<end; i++){
		// printf("%d:%d\n", nAncestors, leafNodes->parentNode[i]);
		nAncestors = add_an_ancestor(ancestor, leafNodes->parentNode[i], nAncestors, maxNumOfAncestors, uts);
		if(nAncestors == maxNumOfAncestors) break;
	}
	int current = 0;
	while(current < nAncestors && nAncestors < maxNumOfAncestors){
		int node = ancestor[current];
		if(node >= internalNodes->nNodes){ STOP("node ID is out of bound");}
		start = internalNodes->startIndex[node];
		end = (node+1 < internalNodes->nNodes) ? internalNodes->startIndex[node+1] : internalNodes->length_parentNode;
		for(i=start; i<end; i++){
			// printf("%d:%d\n", nAncestors, internalNodes->parentNode[i]);
			nAncestors = add_an_ancestor(ancestor, internalNodes->parentNode[i], nAncestors, maxNumOfAncestors, uts);
			if(nAncestors == maxNumOfAncestors) break;
		}
		current++;
	}
	for(i=0; i<nAncestors; i++){
		uts->nodeData[ancestor[i]] = 0;
	}
	return nAncestors;
}

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
){
	int i;
	if(uts == NULL){ STOP("uts == NULL");}
	if(debug){
		for(i=0; i<uts->nNodes; i++) if(uts->nodeData[i] != 0){ STOP("uts->nodeData[i] != 0");}
	}
	if(internalNodes->nNodes != uts->nNodes){ STOP("internalNodes->nNodes != uts->nNodes");}
	int nAncestors = 0;
	for(i=0; i<nInputNodes; i++){
		// printf("%d:%d\n", nAncestors, inputNodeID[i]);
		nAncestors = add_an_ancestor(ancestor, inputNodeID[i], nAncestors, maxNumOfAncestors, uts);
		if(nAncestors == maxNumOfAncestors) break;
	}
	int current = 0;
	while(current < nAncestors && nAncestors < maxNumOfAncestors){
		int node = ancestor[current];
		if(node >= internalNodes->nNodes){ STOP("node ID is out of bound");}
		int start = internalNodes->startIndex[node];
		int end = (node+1 < internalNodes->nNodes) ? internalNodes->startIndex[node+1] : internalNodes->length_parentNode;
		for(i=start; i<end; i++){
			// printf("%d:%d\n", nAncestors, internalNodes->parentNode[i]);
			nAncestors = add_an_ancestor(ancestor, internalNodes->parentNode[i], nAncestors, maxNumOfAncestors, uts);
			if(nAncestors == maxNumOfAncestors) break;
		}
		current++;
	}
	for(i=0; i<nAncestors; i++){
		uts->nodeData[ancestor[i]] = 0;
	}
	return nAncestors;
}

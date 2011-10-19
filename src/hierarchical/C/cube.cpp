/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
*/

#include "cube.hpp"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

LatticeEnumeration::LatticeEnumeration(int *numOfLevels, int ndim, bool display){

	int nLatLevel = 1;
    nDim = ndim;
    length = 1;
    for(int i=0; i<nDim; i++) nLatLevel += (numOfLevels[i]-1);
    for(int i=0; i<nDim; i++) length *= numOfLevels[i];

    level = new int*[length]; for(int i=0; i<length; i++) level[i] = new int[nDim];

    for(int i=0; i<nDim; i++) level[0][i] = 0;
    int beginIndex = 1;
    int availIndex = 1;
    int prevBeginIndex = 0;
    bool gotoEnd = false;

    for(int k=1; k<nLatLevel; k++){
        // Print previous lattice level
        if(display){
            for(int i=prevBeginIndex; i<beginIndex; i++){
                for(int d=0; d<nDim; d++) printf("%d", level[i][d]);
                printf(" ");
            }
            printf("\n");
        }

        // Generate this level from the previous level
        for(int p=prevBeginIndex; p<beginIndex; p++){
            for(int dim=0; dim<nDim; dim++){
                if(level[p][dim] >= numOfLevels[dim]-1) continue;
                // Generate candidate
                for(int d=0; d<nDim; d++) level[availIndex][d] = level[p][d];
                level[availIndex][dim]++;
                // Existence test
                bool exist = false;
                for(int i=beginIndex; i<availIndex; i++){
                    bool equal = true;
                    for(int d=0; d<nDim; d++){
                        if(level[i][d] != level[availIndex][d]){
                            equal = false; break;
                        }
                    }
                    if(equal){
                        exist = true; break;
                    }
                }
                // Update index
                if(!exist){
                    availIndex++;
                    if(availIndex == length){
                        for(int d=0; d<nDim; d++){
                            if(level[availIndex-1][d] != (numOfLevels[d]-1)) STOP_HERE("level[availIndex-1][d] != (numOfLevels[d]-1)");
                        }
                        gotoEnd = true;
                        break;
                    }
                }
            }
            if(gotoEnd) break;
        }
        if(gotoEnd) break;

        // Update indexes
        prevBeginIndex = beginIndex;
        beginIndex = availIndex;
    }

    if(availIndex != length) STOP_HERE("availIndex != length");
    if(display){
        for(int d=0; d<nDim; d++) printf("%d", level[availIndex-1][d]);
        printf("\n");
    }
}

LatticeEnumeration::~LatticeEnumeration(void){
	for(int i=0; i<length; i++)	 delete[] level[i];
	delete[] level;
}

BottomUpCubeProcessor::BottomUpCubeProcessor(DimSchemaHierarchical **cubeDims, int ndim){
	CubeDimension = cubeDims; nDim = ndim;

	verbose = 10; cube = NULL;

    dimSize = new int[nDim];
    dimNumLevels = new int[nDim];
    dimSizeByLevel = new int*[nDim];
    for(int i=0; i<nDim; i++){
        dimSize[i] = CubeDimension[i]->numNodes();
        dimNumLevels[i] = CubeDimension[i]->numLevels();
        dimSizeByLevel[i] = new int[dimNumLevels[i]];
        for(int j=0; j<dimNumLevels[i]; j++){
            dimSizeByLevel[i][j] = CubeDimension[i]->numNodesAtLevel(j);
        }
    }

    NumSubsets = 1;
    for(int i=0; i<nDim; i++) NumSubsets *= dimSize[i];
    NumBaseSubsets = 1;
    for(int i=0; i<nDim; i++) NumBaseSubsets *= dimSizeByLevel[i][0];

    MultiDimLevels = new LatticeEnumeration(dimNumLevels, nDim);

	temp_dimVID = new int[nDim];
}

BottomUpCubeProcessor::~BottomUpCubeProcessor(void){
	if(temp_dimVID != NULL) delete[] temp_dimVID;
	delete[] dimSize;
	delete[] dimNumLevels;
	for(int i=0; i<nDim; i++) delete[] dimSizeByLevel[i];
	delete[] dimSizeByLevel;
	delete MultiDimLevels;
}

void BottomUpCubeProcessor::compute(IBottomUpCubeComputable *cube){
	this->cube = cube;
    if(verbose>0){
        printf("Execute BottomUpCubeProcessor\n"
               "  prepare data structures ...\n");
    }
    if(verbose>0) printf("  initialize computation ...\n");
    cube->initialize(this);

    // Cube Computation:
    if(verbose>0) printf("  process base subsets ...\n");
    processBaseSubsets(); // process the base subsets
    if(verbose>0) printf("  process nonbase subsets ...\n");
    processNonbaseSubsets();// process the non-base subsets
    if(verbose>0) printf("  finalize score cube ...\n");
    cube->finalizeCube(this); // finalize cube
    this->cube = NULL;
}

int BottomUpCubeProcessor::getBestChildSubsetIndices(const int *dimVID, int *childIndex, int size) const {
	int bestDim = getDimWithFewestChildren(dimVID);
	int nChildren = CubeDimension[bestDim]->numChildren(dimVID[bestDim]);
	const int *children = CubeDimension[bestDim]->getChildren(dimVID[bestDim]);
	if(size < nChildren) STOP_HERE("size < nChildren");
	for(int i=0; i<nDim; i++) temp_dimVID[i] = dimVID[i];

	// loop through each child
	for(int c=0; c<nChildren; c++){
		temp_dimVID[bestDim] = children[c];
		childIndex[c] = getSubsetIndex(temp_dimVID);
	}
	return nChildren;
}

void BottomUpCubeProcessor::processBaseSubsets(){

    int* level = new int[nDim];
    // level[i]: the current level of the i-th dimension.
    for(int i=0; i<nDim; i++) level[i] = 0;

    int *levelSize = new int[nDim];
    int **range = new int*[nDim];
    int *dim = new int[nDim];
    int stop;

    fillInLevelSize(levelSize, level);
    fillInRange(range, level);
    stop = getStopIndex(levelSize);

    for(int counter=0; counter<stop; counter++){

        // Compute dimensional index (i.e., dim[])
        computeDimIndex(counter, dim, levelSize, range);

        // Call processOneBaseSubset(...)
        if(verbose>1){
            int index = toLinearIndex(dim,dimSize,nDim);
            printf("    processing cube[%d]: (", index);
            print_row_vector(dim, nDim, stdout, ", ");
            printf(")\n");
        }
        cube->processOneBaseSubset(dim, nDim, this);
    }

    delete[] level;
    delete[] levelSize;
    delete[] range;
    delete[] dim;
}

void BottomUpCubeProcessor::processNonbaseSubsets(){

    int* level;
         // level[i]: the current level of the i-th dimension.
    int *levelSize = new int[nDim];
         // leveSize[i]: the size of the current level at dimension i.
    int **range = new int*[nDim];
    int *dim = new int[nDim];
         // dim[i] the current value at dimension i.
    int inner_stop;

    for(int i=1; i<MultiDimLevels->length; i++){

        level = MultiDimLevels->level[i];

        fillInLevelSize(levelSize, level);
        fillInRange(range, level);
        inner_stop = getStopIndex(levelSize);

        for(int counter=0; counter<inner_stop; counter++){
            // Compute dimensional index (dim[])
            computeDimIndex(counter, dim, levelSize, range);
            if(verbose>1){
                int index = toLinearIndex(dim,dimSize,nDim);
                printf("    processing cube[%d]: (", index);
                print_row_vector(dim, nDim, stdout, ", ");
                printf(")\n");
            }
            // Call processOneNonbaseSubset(...)
            cube->processOneNonbaseSubset(dim, nDim, this);
        }
    }
    delete[] levelSize;
    delete[] range;
    delete[] dim;
}

DimSchemaHierarchical::DimSchemaHierarchical(
	const int *nodeID, const int *parentID, const int nNodes, const int nEdges,
	const DimSchemaHierarchical_LevelType levelType
){
	parentsIndex = new SimpleDagIndex(nodeID, parentID, nNodes, nEdges);
	childrenIndex= new SimpleDagIndex(parentID, nodeID, nNodes, nEdges);
	this->nNodes = nNodes;
	if(parentsIndex->numNextNodes(0) != 0) STOP_HERE("Node 0 is not the root");
	for(int i=1; i<nNodes; i++){
		if(parentsIndex->numNextNodes(i) != 1) STOP_HERE("The hierarchy is not a tree!!");
	}

	int* temp = new int[nNodes];  for(int i=0; i<nNodes; i++) temp[i]  = -1;
	level = new int[nNodes]; for(int i=0; i<nNodes; i++) level[i] = -1;
	int index=0;
	temp[0] = 0;
	nLevels = 0;
	// Root has level 0
	while(index >= 0){
		int node = temp[index]; index--;
		if(node < 0 || node >= nNodes) STOP_HERE("node < 0 || node >= nNodes");
		if(level[node] != -1) STOP_HERE("The hierarchy is not a tree!!");
		if(node == 0) level[node] = 0;
		else{
			const int *parents = parentsIndex->getNextNodes(node);
			int parent = parents[0];
			if(parent < 0 || parent >= nNodes) STOP_HERE("parent < 0 || parent >= nNodes");
			if(level[parent] < 0) STOP_HERE("level[parent] < 0");
			level[node] = level[parent]+1;
		}
		if(nLevels < level[node]+1) nLevels = level[node]+1;
		int nChildren = childrenIndex->numNextNodes(node);
		if(nChildren > 0){
			const int *children = childrenIndex->getNextNodes(node);
			for(int i=0; i<nChildren; i++){
				index++;
				temp[index] = children[i];
			}
		}
	}

	if(levelType == LEAVES_AT_LEVEL_0){
		// Convert to leaves having level 0
		for(int i=0; i<nNodes; i++){
			level[i] = nLevels - level[i] - 1;
			if(childrenIndex->numNextNodes(i) == 0) level[i] = 0;
		}
	}

	nNodesAtLevel = new int[nLevels]; for(int i=0; i<nLevels; i++) nNodesAtLevel[i] = 0;
	for(int i=0; i<nNodes; i++){
		if(level[i] == -1) STOP_HERE("The hierarchy is not a tree!!");
		if(level[i] < 0 || level[i] >= nLevels) STOP_HERE("level[i] < 0 || level[i] >= nLevels");
		nNodesAtLevel[level[i]]++;
	}

	nodeIDsAtLevel = new int*[nLevels];
	for(int i=0; i<nLevels; i++) nodeIDsAtLevel[i] = new int[nNodesAtLevel[i]];
	for(int i=0; i<nNodes; i++) temp[i] = 0;
	for(int i=0; i<nNodes; i++){
		int nodeLevel = level[i];
		if(temp[nodeLevel] < 0 || temp[nodeLevel] >= nNodesAtLevel[nodeLevel]) STOP_HERE("error");
		nodeIDsAtLevel[nodeLevel][temp[nodeLevel]] = i;
		temp[nodeLevel]++;
	}

	delete[] temp;
}

DimSchemaHierarchical::~DimSchemaHierarchical(void){
	delete parentsIndex;
	delete childrenIndex;
	delete[] nNodesAtLevel;
	for(int i=0; i<nLevels; i++) delete[] nodeIDsAtLevel[i];
	delete[] nodeIDsAtLevel;
	delete[] level;
}

DimHierarchyArray::DimHierarchyArray(
	const int *dim,       // size: length
	const int *node_id,	  // size: length
	const int *parent_id, // size: length
	int nDim,
	int length
){
	this->nDim = nDim;
	int *node_id_start   = NULL;
	int *parent_id_start = NULL;
	int nNodes = 0;
	int size = 0;
	int prevDim = -1;
	hierarchy = new DimSchemaHierarchical*[nDim];
	for(int i=0; i<length; i++){
		if(dim[i] != prevDim){
			if(dim[i] != prevDim+1) STOP_HERE("dim[i] != prevDim+1");
			if(i > 0){
				if(size == 0) STOP_HERE("size == 0");
				if(prevDim < 0 || prevDim >= nDim) STOP_HERE("prevDim < 0 || prevDim >= nDim");
				hierarchy[prevDim] = new DimSchemaHierarchical(node_id_start, parent_id_start, nNodes, size);
			}
			node_id_start   = &((const_cast<int*>(node_id))[i]);
			parent_id_start = &((const_cast<int*>(parent_id))[i]);
			size = 0;
			nNodes = 0;
		}
		prevDim = dim[i];
		if(node_id[i]+1 > nNodes) nNodes = node_id[i] + 1;
		if(parent_id[i]+1 > nNodes) nNodes = parent_id[i] + 1;
		size++;
	}

	if(size == 0) STOP_HERE("size == 0");
	if(prevDim != nDim-1) STOP_HERE("prevDim != nDim-1");
	hierarchy[prevDim] = new DimSchemaHierarchical(node_id_start, parent_id_start, nNodes, size);
}

DimHierarchyArray::~DimHierarchyArray(void){
	for(int i=0; i<nDim; i++) delete hierarchy[i];
	delete hierarchy;
}

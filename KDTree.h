/*
 * KDTree.h
 *
 *  Created on: Apr 27, 2017
 *      Author: gleit
 */

#ifndef KDTREE_H_
#define KDTREE_H_

#include "SPPoint.h"
#include "SPConfig.h"
#include "SPBPriorityQueue.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

typedef struct SPDKArray{
	SPPoint** arr;
	int size;
	int** mat;
}KD_ARRAY;

typedef struct sp_kd_tree_node KDTreeNode;

struct sp_kd_tree_node{
	int dim;
	double val;
	KDTreeNode* l;
	KDTreeNode* r;
	SPPoint* data;
};

KDTreeNode createNode(int dim, double val, KDTreeNode* l, KDTreeNode* r, SPPoint* data);

KDTreeNode createLeafNode(SPPoint* data);

bool isLeaf(KDTreeNode node);

int compareVal (const void * A, const void * B);

// A[i,j] has the order of the jth point comparing to the ith dimension of all points in P
// for all i rows: initialize row{sort by dim i}};
//sorting row O(nlog(n)) | all init O(d x nlog(n))
int** initMatrix(KD_ARRAY P);

char* concat(const char *s1, const char *s2);

char* convertIntToStr(int number);

char* buildAddress(int index);

char* buildFeatAddress(int index);

void writeFeatsToFile(FILE* f, SPPoint** feats,int index, int numOfFeats);

KD_ARRAY* InitArray(SPPoint** arr, int size);

void KNearestNeighborSearch(KDTreeNode curr, SPBPQueue* bpq, SPPoint* testPoint);

//this code implements the psuedo code suggested in the work instruction
SPBPQueue* InitKNearestNeighborSearch(KDTreeNode kdNode,SPPoint* testPoint);

int calcSpreadCoor(KD_ARRAY kdArr);

int determineCoor(KD_ARRAY kdArr);

double determineMedianValue(KD_ARRAY kdArr, int dim);

void SplitMatrix(int** A, int* X, int NumOfPoints, KD_ARRAY* SplittedArrays);

// X[i] has 0/1 value of the ith point: 0 if smaller than the median point of the splitting coor, 1 otherwise
//	signRowForSplitting{sign 0/1 = r/l for each cell in row coor};// we can do it since we dont corrupt to original arr)
int* initSortingArr(int** A, int coor, int NumOfPoints);

KD_ARRAY* BuildSplitedArrays (KD_ARRAY P, int* X, int NumOfPoints);

void destroyArr(KD_ARRAY* Arr);

//coor - by which coordinate to split
//returns a tupple of the right and left branches' roots
//O(d x n)
// d = spPCADimension
// n = total no of features of all img in the directory
KD_ARRAY* Split(KD_ARRAY kdArr, int coor);

KDTreeNode InitTree(KD_ARRAY arr);

//init of KDarray
//O(d x nlog(n))
void Init(SPPoint** arr, int size);


#endif /* KDTREE_H_ */

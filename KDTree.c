/*
 * KDTree.c
 *
 *  Created on: Apr 20, 2017
 *      Author: gleit
 */

#include "SPPoint.h"
#include "SPConfig.h"
#include "KDTree.h"
#include "SPBPriorityQueue.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

int SplitCoor = 0;

KDTreeNode createNode(int dim, double val, KDTreeNode* l, KDTreeNode* r, SPPoint* data){
	KDTreeNode* nodeP = (KDTreeNode*) malloc(sizeof(KDTreeNode));
	KDTreeNode node = *nodeP;
	node.dim = dim;
	node.val = val;
	node.l = l;
	node.r = r;
	node.data = data;
	return node;
}

KDTreeNode createLeafNode(SPPoint* data){
	KDTreeNode* nodeP = (KDTreeNode*) malloc(sizeof(KDTreeNode*));
	KDTreeNode node = *nodeP;
	node.dim = -1;
	node.val = -1;
	node.l = NULL;
	node.r = NULL;
	node.data = data;
	return node;
}

bool isLeaf(KDTreeNode node){
	//input kd tree node
	//checks if the node is a leaf
	return node.l == NULL && node.r == NULL;
}

int compareVal (const void * A, const void * B){
	int* arrA = (int*)A;
	int* arrB = (int*)B;
	int a = arrA[1];
	int b = arrB[1];
	if (a < b) return -1;
	if (a == b) return 0;
	 return 1; //if (a >  b)
}

// A[i,j] has the order of the jth point comparing to the ith dimension of all points in P
// for all i rows: initialize row{sort by dim i}};
//sorting row O(nlog(n)) | all init O(d x nlog(n))
int** initMatrix(KD_ARRAY P){
	SP_CONFIG_MSG* msg;
	int NumOfPoints = P.size;
	int dim = spConfigGetPCADim(publicConfig,msg);
	if(*msg != SP_CONFIG_SUCCESS)
		return NULL;///////////elaborate
	int** RowToSort = (int**)malloc(NumOfPoints*2); // will get the value of all the points' cur coor value, will be done for all coordinates

	int** A = (int**)malloc(dim * NumOfPoints);
	if (!A){
		free(RowToSort);
		return NULL;//////////////elaborate
	}

	for (int i=0; i < dim; i++){
		for (int j=0; j < NumOfPoints; j++){
			RowToSort[j][0] = j;       // keep the index of the point, will connect to the point after sorting
			RowToSort[j][1] = spPointGetAxisCoor(P.arr[j],i); //get the jth point's ith coor value
		}
		qsort(RowToSort, NumOfPoints, 2*sizeof(int), compareVal);//sort row with the original values
		for (int k=0; k < NumOfPoints; k++){
			A[i][RowToSort[k][0]] = k; // in the final matrix (A), in the current coordinate we initiate (line i),
			// we want to go point by point from the original point array (P),
			// and place in the matrix in that specific point's cell (A[i][j]),
			// the point's respective placement
			// (ex. "0" if that point has the lowest value in that coor, among all other points)
			// The placement of the jth point for the current coor in on RowToSort[k][1],
			// where curRowVal[k]=j.
		}
	}
	free(RowToSort);
	return A;
}

char* concat(const char *s1, const char *s2)
{
	//input: two strings
	//output: the strings concatinated
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

char* convertIntToStr(int number){
	//input: integer
	//output: a string built from the int
	int length = 0;
	int number1 = number;
	int number2 = number;
	char* strNum;
	while(number1 > 0){
		number1 = number1 / 10;
		length++;
	}
	strNum = (char*) malloc(sizeof(char)*length + 1);
	sprintf(strNum,"%d",number2);
	return strNum;
}

char* buildAddress(int index){
	char* tmp0 = concat(publicConfig.spImagesDirectory,publicConfig.spImagesPrefix);
	char* indexStr = convertIntToStr(index);
	char* tmp1 = concat(tmp0,indexStr);
	char* tmp2 = concat(tmp1,publicConfig.spImagesSuffix);
	free(tmp0);
	free(indexStr);
	free(tmp1);
	return tmp2;
}

char* buildFeatAddress(int index){
	char* tmp0 = concat(publicConfig.spImagesDirectory,publicConfig.spImagesPrefix);
	char* indexStr = convertIntToStr(index);
	char* tmp1 = concat(tmp0,indexStr);
	char* tmp2 = concat(tmp1,".feats");
	free(tmp0);
	free(indexStr);
	free(tmp1);
	return tmp2;
}

void writeFeatsToFile(FILE* f, SPPoint** feats,int index, int numOfFeats){
	// input: file, features, index of the img, number of features to write for that img index.
	// output: writes to file the features
	// the first line is the index, the second is number of features, the rest is the features data
	int i, j, dim;
	char* indexStr = convertIntToStr(index);
	char* numOfFeatsStr = convertIntToStr(numOfFeats);
	fputs(indexStr,f);
	fputs("\n",f);
	fputs(numOfFeatsStr,f);
	fputs("\n",f);
	free(indexStr);
	free(numOfFeatsStr);

	for(i = 0; i < numOfFeats; i++){
		dim = spPointGetDimension(feats[i]);
		for(j = 0; j < dim; j++){
			fprintf(f,"%lf",spPointGetAxisCoor(feats[i], j));
			fputs("\n",f);
		}
	}
}


//Vik
KD_ARRAY InitArray(SPPoint** pointArr, int size){
// input: gets an empty array pointer, mallocs it with the size 'size'
//and starts adding points
	KD_ARRAY kdArr;
	int numberOfImages = publicConfig.spNumOfImages;
	int dim = publicConfig.spPCADimension;
	int i = 0;
	int j = 0;
	int k = 0;
	int numOfFeatures;
	//pointArr = (SPPoint**)malloc(sizeof(SPPoint*)*size);// ??? did you allocate ok??? you are not sure about the only one * in the parenthasis
	FILE* f;
	char* featAddress;
	int index;
	int countPosition = 0; //helper index to write in the array
	for(i = 0; i < numberOfImages; i++){
		featAddress = buildFeatAddress(i);
		f = fopen(featAddress,"r");
		if(f == NULL){
			printf("could not open the file\n"); //need to change to the logger msgs
			printf("the feat Address is %s \n", featAddress);
			return kdArr;
		}
		fscanf(f,"%d",&index);
		fscanf(f,"%d",&numOfFeatures);
		for(j = 0; j < numOfFeatures; j++){
			double* data = (double*)malloc(sizeof(double)*numOfFeatures);
			pointArr[countPosition] = spPointCreate(data, dim,index);
			for(k = 0; k < dim; k++){
				fscanf(f,"%lf",&pointArr[countPosition]->data[k]);
			}
			countPosition++;
		}
		free(featAddress);
	}
	kdArr.size = size;
	kdArr.arr = pointArr;
	kdArr.mat = initMatrix(kdArr);
	return kdArr;
}


void KNearestNeighborSearch(KDTreeNode curr, SPBPQueue* bpq, SPPoint* testPoint){
	// input: kd tree with all the features, a test point, and priority queue to work with
	// output: update the priority queue, so that at the end it has the spknn nearest features
	if(curr == NULL){
		return;
	}
	if(isLeaf(curr)){
		enqueue(bpq,spPointGetIndex(curr->data),sqrt(spPointL2SquaredDistance(curr,testPoint)));
		return;
	}
	if(testPoint->data[curr->dim] <= curr->val){
		KNearestNeighborSearch(curr->l,bpq,testPoint);
	}
	else{
		KNearestNeighborSearch(curr->r,bpq,testPoint);
	}
	if(!spBPQueueIsFull(bpq) || abs(curr->val - testPoint->data[curr->dim]) < spBPQueueMaxValue(bpq)){
		return;// not sure about that. in the instruction they have sayed to search the other half of the tree.
		//this is the way to reach the other half right???
	}
}

//this code implements the psuedo code suggested in the work instruction
SPBPQueue* InitKNearestNeighborSearch(KDTreeNode kdNode,SPPoint* testPoint){
	//input: the kd tree, when it is ready. a point to look for close points to it
	//output: a priority queue with the closest points.
	// init the bpq -
	SPBPQueue* bpq = spBPQueueCreate(publicConfig->spKNN);
	//do the recursive search -
	KNearestNeighborSearch(kdNode,bpq,testPoint);
	return bpq;
	//------------------------------------------------
	//NOTE: YOU NEED TO DESTROY bpq AT SOME POINT IN TIME!!!
	//------------------------------------------------
}

int calcSpreadCoor(KD_ARRAY kdArr){
	//input: array of points
	//output: searches the max spreadand returns the corresponding coordinate
	int maxSpread = -1;
	int coorOfMaxSpread = 0;
	int min;
	int max;
	int tmpSpread;
	int tmpMaxMin;
	int i = 0, j = 0;
	for(i = 0; i < publicConfig->spPCADimension; i++){
		min = kdArr->arr->data[i];
		max = kdArr->arr->data[i];
		for(j = 0; j < kdArr->size; j++){
			tmpMaxMin = kdArr->arr->data[j];
			if( tmpMaxMin < min){
				min = tmpMaxMin;
			}else if(tmpMaxMin > max){
				max = tmpMaxMin;
			}
		}
		tmpSpread = abs(max- min);
		if(tmpSpread > maxSpread){
			maxSpread = tmpSpread;
			coorOfMaxSpread = i;
		}
	}
	return coorOfMaxSpread;
}

int determineCoor(KD_ARRAY kdArr){
	srand(time(NULL));   // should only be called once
	switch(publicConfig.spKDTreeSplitMethod){//do as get
	case 0:
		SplitCoor = calcSpreadCoor(kdArr);
		break;
	case 1:
		SplitCoor = rand() % publicConfig->spPCADimension;//random number between 0 and numer of dimension
		break;
	case 2:
		SplitCoor++;
		break;
		//cases max spread random incremental
	}
	return SplitCoor;
}

double determineMedianValue(KD_ARRAY kdArr, int dim){
	//input: kd array, requested coor of the points
	//output: the med
	int i;
	double median;
	SPBPQueue* bpq = spBPQueueCreate((kdArr->size/2) + kdArr->size%2);
	for(i = 0; i < kdArr->size; i++){
		spBPQueueEnqueue(bpq,i,kdArr->arr[i]->data[dim]);
	}
	median = spBPQueueMaxValue(bpq);
	spBPQueueDestroy(bpq);
	return median;
}

//KD_ARRAY copyArr(KD_ARRAY arr){
//	KD_ARRAY* array = (KD_ARRAY*)malloc(sizeof(KD_ARRAY));
//	KD_ARRAY res = *array;
//	for (int i=0; i < arr.size; i++)
//		res.arr[0]=spPointCopy(arr.arr[0]);
//	return res;
//}

void SplitMatrix(int** A, int* X, int NumOfPoints, KD_ARRAY* SplittedArrays){

	//Part 1: create the maps
	int* map1 = (int*)malloc(NumOfPoints);
	int* map2 = (int*)malloc(NumOfPoints);
	int l = 0; // iterates on P1, will represent it's ordered position
	int r = 0; // iterates on P2, will represent it's ordered position

	for (int i = 0; i < (NumOfPoints); i++) {
		if (X[i] == 0){         // the ith point is on the lth ordered position among P1
			map1[i] = l;
			map2[i] = -1;
			l++;
		}
		else {					// the ith point is on the rth ordered position among P2
			map1[i] = -1;
			map2[i] = r;
			r++;
		}
	}

	//Part 2: populate the matrices
	SP_CONFIG_MSG* msg;
	int dim = spConfigGetPCADim(publicConfig,msg);
	if(*msg != SP_CONFIG_SUCCESS)
		return;
	///////////elaborate

	int** A1 = (int**)malloc(dim * (NumOfPoints/2 + NumOfPoints%2));
	int** A2 = (int**)malloc(dim * NumOfPoints/2);

	l=0;
	r=0;
	for (int m=0; m < dim; m++){
		for (int n=0; n < NumOfPoints; n++){
			if(map1[n] != -1){            //suppose to be on A1
				A1[m][l] = A[m][n];
				l++;
			}
			else{                         //suppose to be on A2
				A2[m][r] = A[m][n];
				r++;
			}
		}

	}
	free(map1);
	free(map2);
	SplittedArrays[0].mat = A1;
	SplittedArrays[1].mat = A2;
}

// X[i] has 0/1 value of the ith point: 0 if smaller than the median point of the splitting coor, 1 otherwise
//	signRowForSplitting{sign 0/1 = r/l for each cell in row coor};// we can do it since we dont corrupt to original arr)
int* initSortingArr(int** A, int coor, int NumOfPoints){
	int* X = (int*)malloc(NumOfPoints);
	for (int j=0; j<NumOfPoints; j++) {
		if(A[coor][j] > (NumOfPoints/2+NumOfPoints%2)) // The smaller half of the images (rounded up) will be in the Left array, the bigger ones on the Right
			X[j]=1;
		else
			X[j]=0;
	}
	return X;
}

KD_ARRAY* BuildSplitedArrays (KD_ARRAY P, int* X, int NumOfPoints){

	KD_ARRAY* SplittedArr = (KD_ARRAY*)malloc(sizeof(KD_ARRAY)*2);

	KD_ARRAY* P1 = (KD_ARRAY*)malloc(sizeof(KD_ARRAY));
	SPPoint** p1 = (SPPoint**)malloc(sizeof(SPPoint*) * (NumOfPoints/2 + NumOfPoints%2));
	P1->arr = p1;

	KD_ARRAY* P2 = (KD_ARRAY*)malloc(sizeof(KD_ARRAY));
	SPPoint** p2 = (SPPoint**)malloc(sizeof(SPPoint*) * (NumOfPoints/2));
	P2->arr = p2;

	int l = 0; // iterate on the left array P1
	int r = 0; // iterate on the right array P2
	for (int j = 0; j < (NumOfPoints); j++) {
		if (X[j] == 0){         // the jth point's value at of coor <= median - will go on the left branch
			P1->arr[l] = spPointCopy(P.arr[j]);
			l++;
		}
		else {					// the jth point's value at of coor > median - will go on the left branch
			P2->arr[r] = spPointCopy(P.arr[j]);
			r++;
		}
	}
	SplittedArr[0]=*P1;
	SplittedArr[1]=*P2;
	return SplittedArr;
}

void destroyArr(KD_ARRAY* Arr){
	if(!Arr)
		return;
	KD_ARRAY arr=*Arr;
	if(*arr.arr){
		for (int i = 0; i < arr.size; i++)
			spPointDestroy(arr.arr[0]);
		free(arr.arr);
	}
	if(arr.mat)
		free(arr.mat);
	free(Arr);
}

//coor - by which coordinate to split
//returns a tupple of the right and left branches' roots
//O(d x n)
// d = spPCADimension
// n = total no of features of all img in the directory
KD_ARRAY* Split(KD_ARRAY kdArr, int coor){
	int* X = initSortingArr(kdArr.mat, coor, kdArr.size);     // X[i] has 0/1 value of the ith point: 0 if smaller than the median point of the splitting coor, 1 otherwise
	KD_ARRAY* splittedArr = BuildSplitedArrays(kdArr, X, kdArr.size);  // Splits the original Array into left and right
	SplitMatrix(kdArr.mat, X, kdArr.size, splittedArr);    //Assign the splitted arrays their respective matrix for the next run
	free(X);
	return splittedArr;
}

void destroyTree(KDTreeNode* node){
	spPointDestroy(node->data);
	if (node->l != NULL)
		destroyTree(node->l);
	if (node->r != NULL)
		destroyTree(node->r);
}

KDTreeNode InitTree(KD_ARRAY arr){

	//Recursion stopping condition
	if (arr.size == 1)
		return createLeafNode(arr.arr[0]);

	//Recursion step
	int dim = determineCoor(arr);
	double val = determineMedianValue(arr, dim);
	KD_ARRAY* splittedArr = Split(arr, dim);
	KDTreeNode l = InitTree(splittedArr[0]);
	KDTreeNode r = InitTree(splittedArr[1]);
	destroyArr(&splittedArr[0]);
	destroyArr(&splittedArr[1]);
	SPPoint* data = NULL;
	return createNode(dim, val, &l, &r, data);
}

//init of KDarray
//O(d x nlog(n))
void Init(SPPoint** arr, int size){
	KD_ARRAY* kdArr = InitArray(arr,size);
	InitTree(*kdArr);
}


#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../SPPoint.h"
#include "unit_test_util.h" //SUPPORTING MACROS ASSERT_TRUE/ASSERT_FALSE etc..
#include "../SPConfig.h"  //to delete

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

void printMatrix(KD_ARRAY arr){
	printf("/nMatrix:/n");
	SP_CONFIG_MSG* msg;
	int dim = spConfigGetPCADim(publicConfig,msg);
	for (int i = 0; i < dim ; ++i){
		for (int j = 0; j < arr.size; ++j)
			printf("%d, ",arr.mat[i][j]);
		printf("/n");
	}
}

void printPointsArr(SPPoint** arr, int size, int dim){
	printf("/nPoints array:/n");
	for (int i = 0; i < size; ++i){
		for (int j = 0; j < dim; ++j)
			printf("%d, ", spPointGetAxisCoor(arr[i], j));
		printf("/n");
	}
}

void PrintKDArrays(KD_ARRAY arr, char* name){
	printf("Array name: %s/n", name);
	SP_CONFIG_MSG* msg;
	int dim = spConfigGetPCADim(publicConfig,msg);
	printf("size = %d/n", arr.size);
	printPointsArr(arr.arr, arr.size, dim);
	printMatrix(arr);
	printf("-----/n/n");
	fflush(NULL);
}



bool IsIdenticalKDArrays(KD_ARRAY arr1, KD_ARRAY arr2){
	if(arr1.size != arr2.size)
		return false;
	for (int k = 0; k < arr1.size; ++k)
		if (spPointL2SquaredDistance(arr1.arr[k], arr2.arr[k]) != 0)
			return false;
	SP_CONFIG_MSG* msg;
	int dim = spConfigGetPCADim(publicConfig,msg);
	for (int i = 0; i < dim ; ++i)
		for (int j = 0; j < arr1.size; ++j)
			if(arr1.mat[i][j] != arr2.mat[i][j])
				return false;
	return true;
}

KD_ARRAY BuildTestArr(){

	SPPoint** pointArr = (SPPoint**)malloc(sizeof(SPPoint*)*3);
	KD_ARRAY kdArr = (KD_ARRAY)malloc(sizeof(KD_ARRAY));
	SPPoint** pointArrTest = (SPPoint**)malloc(sizeof(SPPoint*)*3);

	double data1[] = {3.000,0.100,178.000};
	double data2[] = {2.000,1.500,2.000};
	double data3[] = {1.000,1.500,1.000};

	pointArrTest[0] = spPointCreate(data1,3,1);
	pointArrTest[1] = spPointCreate(data2,3,1);
	pointArrTest[2] = spPointCreate(data3,3,1);

	kdArr.arr = pointArrTest;
	kdArr.size=3;
	kdArr.mat = initMatrix(kdArr);

	return kdArr;
}

KD_ARRAY BuildSmallTestArr(){
//	{1,1.5,1}
}

static bool SplitTest() {
	KD_ARRAY arrTest = BuildTestArr();
	KD_ARRAY* splittedarrays = Split(arrTest, 0);
	KD_ARRAY* splittedarraysTest = BuildTestSplittedArray();
	printKDArray(splittedarraysTest[0], "splittedarraysTest[0]");
	printKDArray(splittedarraysTest[1], "splittedarraysTest[1]");
	printKDArray(splittedarrays[0], "splittedarrays[0]");
	printKDArray(splittedarrays[1], "splittedarrays[1]");
	bool ID1 = IsIdenticalKDArrays(splittedarraysTest[0], splittedarrays[0]);
	bool ID2 = IsIdenticalKDArrays(splittedarraysTest[1], splittedarrays[1]);
	ASSERT_TRUE(ID1 && ID2);
	destroyArr(&arrTest);
	destroyArr(&splittedarrays[0]);
	destroyArr(&splittedarrays[1]);
	destroyArr(&splittedarraysTest[0]);
	destroyArr(&splittedarraysTest[1]);
	fflush(NULL);
	return true;
}

SPPoint* BuildTestPoint(){
	double data[3] = {1,1.5,1};
	return spPointCreate(data, 3, 1);
}

KDTreeNode PrintTree(KD_ARRAY arr){
	printKDArray(arr);

	//Recursion stopping condition
	if (arr.size == 1){
		printf("last node's point");
		printPointsArr(arr.arr[0]);
		return createLeafNode(arr.arr[0]);
	}

	//Recursion step
	int dim = determineCoor(arr);
	double val = determineMedianValue(arr, dim);
	KD_ARRAY* splittedArr = Split(arr, dim);
	printKDArray(splittedArr[0]);
	printKDArray(splittedArr[1]);
	KDTreeNode l = InitTree(splittedArr[0]);
	KDTreeNode r = InitTree(splittedArr[1]);
	destroyArr(&splittedArr[0]);
	destroyArr(&splittedArr[1]);
+	SPPoint* data = NULL;
+	return createNode(dim, val, &l, &r, data);
}

static bool InitTreeTest() {
	SPPoint* TestPoint = BuildTestPoint();         //{1,1.5,1}
	KD_ARRAY smallArrTest = BuildSmallTestArr();  //{1,1.5,1}
	KD_ARRAY arrTest = BuildTestArr();            //{3,0.1,178};{2,1.5,2};{1,1.5,1};
	KDTreeNode nodeS = InitTree(smallArrTest);    //{1,1.5,1}
	KDTreeNode node = InitTree(arrTest);          //{1,1.5,1} < {2,1.5,2} > {3,0.1,178}
	KDTreeNode T1 = PrintTree(smallArrTest);
	KDTreeNode T2 = PrintTree(arrTest);
	ASSERT_TRUE(nodeS.l == NULL && nodeS.r == NULL && spPointL2SquaredDistance(TestPoint, node.data) == 0);    //{1,1.5,1},{1,1.5,1}
	ASSERT_TRUE(nodeS.l != NULL && nodeS.r != NULL && spPointL2SquaredDistance(TestPoint, node.data) == 1);    //{1,1.5,1},{2,1.5,2}
	ASSERT_TRUE(node.l == NULL && node.r == NULL && spPointL2SquaredDistance(TestPoint, node.data) == 0);//{1,1.5,1},{1,1.5,1}
	fflush(NULL);
	spPointDestroy(TestPoint);
	destroyArr(&arrTest);
	destroyArr(&smallArrTest);
	destroyTree(&nodeS);
	destroyree(&node);
	destroyree(&T1);
	destroyree(&T2);
	return true;
}

static bool concatTest() {
	char* s1 = "Hello";
	char* s2 = "World";
	char* s3 = concat(s1,s2);
	ASSERT_TRUE(strcmp(s3,"HelloWorld") == 0);
	ASSERT_TRUE(strlen(s3)==10);
	return true;
}

static bool convertIntToStrTest(){
	int num1 = 123;
	int num2 = 3;
	int num3 = 1234567890;
	char* str1 = convertIntToStr(num1);
	char* str2 = convertIntToStr(num2);
	char* str3 = convertIntToStr(num3);
	ASSERT_TRUE(strcmp(str1,"123") == 0);
	ASSERT_TRUE(strcmp(str2,"3") == 0);
	ASSERT_TRUE(strcmp(str3,"1234567890") == 0);
	free(str1);
	free(str2);
	free(str3);
	return true;
}

static bool buildAddressTest(){
	publicConfig.spImagesDirectory = "images/";
	publicConfig.spImagesPrefix = "img";
	publicConfig.spImagesSuffix = ".jpg";
	char* address1 = buildAddress(1);
	char* address2 = buildAddress(2);
	char* address3 = buildAddress(3);
	ASSERT_TRUE(strcmp(address1,"images/img1.jpg") == 0);
	ASSERT_TRUE(strcmp(address2,"images/img2.jpg") == 0);
	ASSERT_TRUE(strcmp(address3,"images/img3.jpg") == 0);
	free(address1);
	free(address2);
	free(address3);
	return true;
}

static bool buildFeatAddressTest(){
	char* address1 = buildFeatAddress(1);
	char* address2 = buildFeatAddress(2);
	char* address3 = buildFeatAddress(3);
	ASSERT_TRUE(strcmp(address1,"images/img1.feats") == 0);
	ASSERT_TRUE(strcmp(address2,"images/img2.feats") == 0);
	ASSERT_TRUE(strcmp(address3,"images/img3.feats") == 0);
	free(address1);
	free(address2);
	free(address3);
	return true;
}

static bool writeFeatsToFileTest(){
	//allocating five points
	//point arr is the features
	int i;
	int j;
	FILE* f_feats = fopen("f_feats.txt","w+");
	SPPoint** pointArr = (SPPoint**)malloc(sizeof(SPPoint*)*5);
	for(i = 0; i < 5; i++){
		double* data = (double*)malloc(sizeof(double)*5);
		for(j = 0; j < 5; j++){
			data[j] =  i*j;
		}
		pointArr[i] = spPointCreate(data,5,2);
	}
	writeFeatsToFile(f_feats,pointArr,2,5);
//	identicalFiles("f_feats.txt","f_test.txt") ? printf("identical\n") : printf("not indentical\n");
//DELETE AFTER
	ASSERT_FALSE(identicalFiles("f_feats.txt","f_test.txt"));
	for(i = 0; i < 5; i++){
		spPointDestroy(pointArr[i]);
	}
	free(pointArr);
	return true;
}

static bool InitArrayTest(){
	//this test depend on a pre-written feats files
	publicConfig.spNumOfImages = 3;
	publicConfig.spPCADimension = 3;
	SPPoint** pointArr = (SPPoint**)malloc(sizeof(SPPoint*)*9);
	KD_ARRAY kdArr = InitArray(pointArr,9);
	SPPoint** pointArrTest = (SPPoint**)malloc(sizeof(SPPoint*)*9);
	double data1[] = {0.000,1.000,2.000};
	double data2[] = {3.000,4.000,5.000};
	double data3[] = {6.000,7.000,8.000};
	double data4[] = {0.000,1.000,2.000};
	double data5[] = {3.000,4.000,5.000};
	double data6[] = {6.000,7.000,8.000};
	double data7[] = {0.000,1.000,2.000};
	double data8[] = {3.000,4.000,5.000};
	double data9[] = {6.000,7.000,8.000};
	pointArrTest[0] = spPointCreate(data1,3,1);
	pointArrTest[1] = spPointCreate(data2,3,1);
	pointArrTest[2] = spPointCreate(data3,3,1);
	pointArrTest[3] = spPointCreate(data4,3,2);
	pointArrTest[4] = spPointCreate(data5,3,2);
	pointArrTest[5] = spPointCreate(data6,3,2);
	pointArrTest[6] = spPointCreate(data7,3,3);
	pointArrTest[7] = spPointCreate(data8,3,3);
	pointArrTest[8] = spPointCreate(data9,3,3);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[0],kdArr.arr[0])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[1],kdArr.arr[1])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[2],kdArr.arr[2])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[3],kdArr.arr[3])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[4],kdArr.arr[4])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[5],kdArr.arr[5])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[6],kdArr.arr[6])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[7],kdArr.arr[7])==0);
	ASSERT_TRUE(spPointL2SquaredDistance(pointArrTest[8],kdArr.arr[8])==0);
	//isIdentical()
	return true;
}

int mainTreeTest(){
	RUN_TEST(concatTest);
	RUN_TEST(convertIntToStrTest);
	RUN_TEST(buildAddressTest);
	RUN_TEST(buildFeatAddressTest);
	RUN_TEST(writeFeatsToFileTest);
	RUN_TEST(InitArrayTest);
	RUN_TEST(SplitTest);
	RUN_TEST(InitTreeTest);
	return 0;
}



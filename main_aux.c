#include <stdbool.h>
#include "main_aux.h"
#include "SPConfig.h"
#include "SPLogger.h"
#include "SPPoint.h"
#include <string.h>
#include "unit_tests\\unit_test_util.h"
#include <stdio.h>
#include <stdlib.h>
#include "SPBPriorityQueue.h"
#include "KDTree.h"

void runTests(){
	mainLoggerTest();
	mainConfigTest();
	mainTreeTest();
}

bool manageCMSG(SP_CONFIG_MSG* Cmsg){
	if (Cmsg)
		if(*Cmsg != SP_CONFIG_SUCCESS)
			return false;
	return true;
}

bool manageLMSG(SP_LOGGER_MSG* Smsg){
	if (Smsg) //always true if the code is correct
		if(*Smsg != SP_LOGGER_SUCCESS){
			printf("Logger creation failed");
			return false;
		}
	return true;
}

char* getConfigFileName(int argc, char* argv[]){
	if (argc > 1)
		return argv[1];
	return "spcbir.config";
}

void freeAll(SPLogger logger){
	spConfigDestroy(publicConfig);
	spLoggerDestroy(logger);
}

int spExtract(){
	FILE* f;
	int i = 0;
	int size = 0;
	int j;
	int l = getSPNumOfImages();
	int* numOfFeatures;
	SPPoint** tmpFeats;
	char* address;
	char* featAddress;
	ImageProc(publicConfig);//did I reach the function right???
		for(i = 0; i < l; i++){
			address = buildAddress(i);
			tmpFeats = getImageFeatures(address,i,numOfFeatures);//need to destroy that
			free(address);
			size = size + *numOfFeatures;
			//open a new file:
			featAddress = buildFeatAddress(i);
			f = fopen(featAddress,"w+");
			free(featAddress);
			if(f == NULL){
				//need to print here with guy's logger.
				return 0;
			}
			//write everything to a file:
			writeFeatsToFile(f,tmpFeats,i, *numOfFeatures); 
			//--------------------------------------------------------------
				for(i = 0; i < *numOfFeatures; i++){
					spPointDestroy(tmpFeats[i]);
			}
	//--------------------------------------------------------------
	free(tmpFeats);//not sure about this solution as it is ineffective but the best i came up with
	fclose(f);
	}
	return size;
}
/*
int receiveCommand(){
	char* command;
	printf("Enter 'search' to use the search. enter 'exit' to finish\n");
	scanf("%s", command);
	while(true){
		command = (char*)malloc(sizeof(char)*1024);
		if(strcmp(command,"exit") == 0){
			free(command);
			return 0;
		}
		else if(strcmp(command,"search")){
			free(command);
			return 1;
		}
		else{
			printf("The command you entered was not valid\nPlease try again -\n");
			free(command);
		}
	}
	return 0;
}*/

int receiveQuery(char* buffer){
	printf("please enter an image path:\n");
	scanf("%s", buffer);
	if(strcmp(buffer,"<>") == 0){
		return 1;
	}
	return 0;
}

int* findSimilarImages(KDTreeNode curr,char* query){
	int* numOfFeatures;
	SPPoint** features = getImageFeatures(query,0,numOfFeatures);
	SPBPQueue* bpq = spBPQueueCreate(getSpKNN());
	int* imgInstance = (int*)malloc(sizeof(int)*getSPNumOfImages());
	int i = 0;
	BPQueueElement* res;
	for(i = 0; i < *numOfFeatures; i++){
		KNearestNeighborSearch(curr, bpq, features[i]);
		while(!spBPQueueIsEmpty(bpq)){
			res = (BPQueueElement*)malloc(sizeof(BPQueueElement));
			spBPQueuePeek(bpq,res);
			imgInstance[res->index]--;//when we dequeue we get the minimum number. but we want the max. so I reverse the array.
			free(res);
			spBPQueueDequeue(bpq);
		}
	}
	i = 0;
	for(i = 0; i < *numOfFeatures; i++){
		free(features[i]);
	}
	free(features);
	spBPQueueDestroy(bpq);
	return imgInstance;
}

void showImages(bool showOrNot,int* imgSimilarityRankArr, char* queryPath){
	SPBPQueue* bpq = spBPQueueCreate(getSpKNN());
	BPQueueElement* res;
	char* addr;

	int j;
	for(j = 0; j < getSPNumOfImages(); j++){
		spBPQueueEnqueue(bpq,j, imgSimilarityRankArr[j]);
	}

	if(showOrNot){
		int i = 0;
		int l = getSPNumOfSimilarImages();
		for(i = 0; i < l; i++){
			res = (BPQueueElement*)malloc(sizeof(res));
			spBPQueuePeek(bpq,res);
			spBPQueueDequeue(bpq);
			addr = buildAddress(res->index*(-1));//remember that all of the numbers are with minus sign, since the priority queue, is minimum priority queue
			showImage(addr);
			free(addr);
			free(res);
		}
	}else{
		printf("Best candidates for - %s - are:\n",queryPath);
		int i = 0;
		int l = getSPNumOfSimilarImages();
		for(i = 0; i < l; i++){
			res = (BPQueueElement*)malloc(sizeof(res));
			spBPQueuePeek(bpq,res);
			spBPQueueDequeue(bpq);

			addr = buildAddress(res->index*(-1));//remember that all of the numbers are with minus sign, since the priority queue, is minimum priority queue
			printf("%s\n",addr);
			free(addr);
			free(res);
		}
	}
}

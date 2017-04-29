#include <stdbool.h>
#include "main_aux.h"

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

void freeAll(SPConfig config, SPLogger logger){
	spConfigDestroy(config);
	spLoggerDestroy(logger);
}

int spExtract(){
	FILE* f;
	int i = 0;
	int size = 0;
	int j;
	int l = publicConfig->spNumOfImages;
	int* numOfFeatures;
	SPPoint** tmpFeats;
	char* address;
	char* featAddress;
	ImageProc(publicConfig);
		for(i = 0; i < l; i++){
			address = buildAddress(i);
			tmpFeats = getImageFeatures(address,i,numOfFeatures);//need to destroy that
			free(address);
			size = size + *numOfFeatures;
			//open a new file:
			featAddress = buildFeatAddress(i);
			f = fopen(featAddress,"w+");
			free(featAddress);
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

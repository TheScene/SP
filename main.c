//#include "sp_image_proc_util.h"
#include "SPConfig.h"
#include "main_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main_aux.h"
#include <stdbool.h>
#include "SPPoint.h"

SP_CONFIG_MSG* Cmsg;
SP_LOGGER_MSG* Lmsg;

#define IsConfigErrExit 	if (!manageCMSG(&Cmsg)){printf("Exiting..."); return 1;}
#define IsLoggerErrExit 	if (!manageCMSG(&Lmsg)){printf("Exiting..."); return 1;}


int main(int argc, char** argv){

	// PART A: Tests
	//runTests();
	//return 0;

	// PART B - initiating
	// B1: config
	char* ConfigFileName =  argv[1];
	SPConfig* config = spConfigCreate(ConfigFileName, Cmsg);


	// B2: logger
	//int loggerLvl = (int)getSPLoggerLevel();
	*Lmsg = spLoggerCreate(getSPLoggerFilename(),getSPLoggerLevel());

	//IsLoggerErrExit

	// B3: features //Vik to verify
	int size;
	if (getSPExtractionMode()){
		size = spExtract(config, Cmsg); //including verifying DONE and saving to directory
	}
	SPPoint** arr = (SPPoint**)malloc(sizeof(SPPoint*)*size);
	
	//Init(arr,size);//inits the kdTree
	//IsConfigErrExit
	//initDataStructures(config, Cmsg); should be inside the init
	//IsConfigErrExit
	int* resImages;
	int result;
	// PART C - Query
	while(true){//why is it error???
		char* query = (char*)malloc(sizeof(char)*1024);
		result = receiveQuery(query);
		if(result == 1){
			free(query);
			break;
		}
		resImages = findSimilarImages(query);
		showImages(getSPMinimalGUI(), resImages,query);
	}

	// PART D - free
	freeAll(logger);
	return 0;
}


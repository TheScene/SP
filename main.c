//#include "sp_image_proc_util.h"
#include "SPConfig.h"
#include "main_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main_aux.h"
#include <stdbool.h>

SP_CONFIG_MSG* Cmsg;
SP_LOGGER_LEVEL* Lmsg;

#define IsConfigErrExit 	if (!manageCMSG(&Cmsg)){printf("Exiting..."); return 1;}
#define IsLoggerErrExit 	if (!manageCMSG(&Lmsg)){printf("Exiting..."); return 1;}


int main(int argc, char** argv){

	/* PART A: Tests */
	runTests();
	return 0;

	/* PART B - initiating */
	// B1: config
	char* ConfigFileName =  getConfigFileName(argc, argv);
	SPConfig* config = spConfigCreate(ConfigFileName, Cmsg);
	IsConfigErrExit

	// B2: logger
	Lmsg = spLoggerCreate(config->spLoggerFilename, config->spLoggerLevel);
	IsLoggerErrExit

	// B3: features //Vik to verify
	int size;
	if (config->spExtractionMode)
		size = spExtract(config, Cmsg); //including verifying DONE and saving to directory
	SPPoint** arr = (SPPoint**)malloc((SPPoint*)*size);
	
	Init(arr,size);//inits the kdTree
	IsConfigErrExit
	initDataStructures(config, Cmsg);
	IsConfigErrExit

	/* PART C - Query */

	While (true){
//		comand = receiveComand(msg);
//		manageMSG;
//		comand? continue : return 0;
//		resImages = findSimilarImages(query, msg);
//		manageMSG;
//		showImages(config.spMinimalGUI, resImages, msg);
//		manageMSG;
	}

	/* PART D - free */
	freeAll(config, logger);
	return 0;
}

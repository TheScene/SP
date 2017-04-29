/*
 * spKDArray_test_util.h
 *
 *  Created on: Apr 29, 2017
 *      Author: viktor
 */

#include "SPPoint.h"
#include <stdbool.h>

typedef struct SPDKArray{
	SPPoint** arr;
	int size;
	int** mat;
}KD_ARRAY;

typedef enum sp_config_msg_t {
	SP_CONFIG_MISSING_DIR,
	SP_CONFIG_MISSING_PREFIX,
	SP_CONFIG_MISSING_SUFFIX,
	SP_CONFIG_MISSING_NUM_IMAGES,
	SP_CONFIG_CANNOT_OPEN_FILE,
	SP_CONFIG_ALLOC_FAIL,
	SP_CONFIG_INVALID_INTEGER,
	SP_CONFIG_INVALID_STRING,
	SP_CONFIG_INVALID_ARGUMENT,
	SP_CONFIG_INDEX_OUT_OF_RANGE,
	SP_CONFIG_NOT_A_VARIABLE,
	SP_CONFIG_SUCCESS
} SP_CONFIG_MSG;

typedef struct sp_config_t{
	char* spImagesDirectory;
	char* spImagesPrefix;
	char* spImagesSuffix;
	int spNumOfImages;
	int spPCADimension;
	char* spPCAFilename;
	int spNumOfFeatures;
	bool spExtractionMode;
	int spNumOfSimilarImages;
	//SP_KDTREE_SPLIT_METHOD spKDTreeSplitMethod;
	int spKNN;
	bool spMinimalGUI;
	//SP_LOGGER_LEVEL spLoggerLevel;
	char* spLoggerFilename;
}SPConfig;

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "unit_test_util.h" //SUPPORTING MACROS ASSERT_TRUE/ASSERT_FALSE etc..
#include "../SPConfig.h"

char* imagePath;
char* pcaPath;
const char* ValidFilename = "configTest5.config"; //a valid config file, included in the submitted file
const char* NonValidFilename = "configTest6.config"; //a config file with non-valid content, included in the submitted file
const char* NotAFilename = "blaaa.config"; //not such file exists in the submitted file


//Checks basic errors and their msg
static bool basicConfigCreateTest() {
	SPConfig config;
	SP_CONFIG_MSG Omsg = SP_CONFIG_SUCCESS;
	SP_CONFIG_MSG* msg = &Omsg;

	config = spConfigCreate(NULL, msg);
	ASSERT_TRUE(*msg == SP_CONFIG_INVALID_ARGUMENT);
	spConfigDestroy(config);

	config = spConfigCreate(NotAFilename, msg);
	ASSERT_TRUE(*msg == SP_CONFIG_CANNOT_OPEN_FILE);
	spConfigDestroy(config);
	return true;
}


//Checks the getters and assignment of default values
static bool ConfigValidVarTest() {
	SPConfig config;
	SP_CONFIG_MSG Omsg = SP_CONFIG_SUCCESS;
	SP_CONFIG_MSG* msg = &Omsg;
	config = spConfigCreate(ValidFilename, msg);

	ASSERT_TRUE(spConfigIsExtractionMode(config, msg) == true);
	ASSERT_TRUE(*msg == SP_CONFIG_SUCCESS);

	ASSERT_TRUE(spConfigMinimalGui(config, msg) == false);
	ASSERT_TRUE(*msg == SP_CONFIG_SUCCESS);

	ASSERT_TRUE(spConfigGetNumOfImages(config, msg) == 867);
	ASSERT_TRUE(*msg == SP_CONFIG_SUCCESS);

	ASSERT_TRUE(spConfigGetNumOfFeatures(config, msg) == 100);
	ASSERT_TRUE(*msg == SP_CONFIG_SUCCESS);

	ASSERT_TRUE(spConfigGetPCADim(config, msg) == 15);
	ASSERT_TRUE(*msg == SP_CONFIG_SUCCESS);

	ASSERT_TRUE(spConfigGetImagePath(imagePath, config, 1) == SP_CONFIG_SUCCESS);
	ASSERT_TRUE(strcmp(imagePath,"abcd/BFF.mhbjn.jpg") == 0);

	ASSERT_TRUE(spConfigGetPCAPath(pcaPath, config) == SP_CONFIG_SUCCESS);
	ASSERT_TRUE(strcmp(pcaPath,"abcd/pca.yml") == 0);

	spConfigDestroy(config);
	return true;
}


//Checks the behavior when given a config file with non valid variables
static bool ConfigNonValidVarTest() {
	SPConfig config;
	SP_CONFIG_MSG Omsg = SP_CONFIG_SUCCESS;
	SP_CONFIG_MSG* msg = &Omsg;
	config = spConfigCreate(NonValidFilename, msg);

	ASSERT_TRUE(*msg == SP_CONFIG_INVALID_INTEGER);

	spConfigDestroy(config);
	return true;
}


int main() {
	RUN_TEST(basicConfigCreateTest);
	RUN_TEST(ConfigValidVarTest);
	RUN_TEST(ConfigNonValidVarTest);
	return 0;
}

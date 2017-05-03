// He who uses the functions below should deal with the resulting value stored in msg

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include "SPConfig.h"

#define MAX_LEN 1025 //1024 string including the terminating null character
#define numberOfArguements 14 

//-------------------------------------------------------------
//TO DO - need to delete it later
int k = 0;
#define pr 	{printf("%d", k); k++; fflush(NULL);}
#define cu 	printf("cur"); fflush(NULL);
//-------------------------------------------------------------

FILE* configFile;
char* sysVarsControl[numberOfArguements] = {"spImagesDirectory","spImagesPrefix","spImagesSuffix","spNumOfImages","spPCADimension","spPCAFilename","spNumOfFeatures","spExtractionMode","spNumOfSimilarImages","spKDTreeSplitMethod","spKNN","spMinimalGUI","spLoggerLevel","spLoggerFilename" };
int sysVarsFlag[numberOfArguements] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int currentVarIndex;




/** A type used for defining the configuration file**/
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
	SP_KDTREE_SPLIT_METHOD spKDTreeSplitMethod;
	int spKNN;
	bool spMinimalGUI;
	SP_LOGGER_LEVEL spLoggerLevel;
	char* spLoggerFilename;
}*SPConfig;


SPConfig* publicConfig;

//needed only for debuggin - since we work on windows and here we need to write a new "getline"
size_t getline(char **lineptr, size_t *n, FILE *stream) {
    char *bufptr = NULL;
    char *p = bufptr;
    size_t size;
    int c;

    if (lineptr == NULL) {
        return -1;
    }
    if (stream == NULL) {
        return -1;
    }
    if (n == NULL) {
        return -1;
    }
    bufptr = *lineptr;
    size = *n;

    c = fgetc(stream);
    if (c == EOF) {
        return -1;
    }
    if (bufptr == NULL) {
        bufptr = malloc(128);
        if (bufptr == NULL) {
            return -1;
        }
        size = 128;
    }
    p = bufptr;
    while(c != EOF) {
        if ((p - bufptr) > (size - 1)) {
            size = size + 128;
            bufptr = realloc(bufptr, size);
            if (bufptr == NULL) {
                return -1;
            }
        }
        *p++ = c;
        if (c == '\n') {
            break;
        }
        c = fgetc(stream);
    }

    *p++ = '\0';
    *lineptr = bufptr;
    *n = size;

    return p - bufptr - 1;
}

bool isLineComment(char* line){
	//input: line in the config
	//output: checks if the line is a comment. return true if yes. false otherwise.
	int i = 0;
	while(line[i] == ' ' && i<strlen(line)){
		i++;
	}
	if(line[i] == '#'){
		return true;
	}
	return false;
}

bool isLineEmpty(char* line){
	//input not a comment line
	//output true if line is empty
	int i = 0;
	if(line[i] == '\n'){
		return true;
	}
	while( i<strlen(line) && line[i] == ' '){
		i++;
	}
	if(i == strlen(line)){ //the second condition is weird but it meant to work when line is uninitialized
		return true;
	}
	return false;
}

bool spacesIsValid(char* line){
	//input: not an empty line in the config, and not a comment
	//output: true if the line has the right format
	int i = 0;
	int length = strlen(line);
	while(i<length && line[i]==' '){
		i++;
	}
	while(i<length && line[i]!=' ' && line[i]!='='){
		i++;
	}
	if(i == length){
		return false;
	}
	while(i<length && line[i] == ' '){
		i++;
	}
	if(i == length){
		return false;
	}
	if(line[i]!='='){
		return false;
	}
	i++;
	while(i<length && line[i] == ' '){
		i++;
	}
	if(i == length){
		return false;
	}
	while(i<length && line[i] != ' '){
		i++;
	}
	if(i == length){
		return true;
	}
	while(i<length && line[i] == ' '){
		i++;
	}
	if(i == length){
		return true;
	}
	if(line[i] !=' '){
		return false;
	}
	return true;
}

char* getTheVariable(char* line){
	//input: a valid spaces line from the config
	//output: gets the variable name and returns it
	int i = 0, length, start, end;
	char* var;
	while(line[i]==' '){
		i++;
	} 
	start = i;
	while(line[i]!=' ' && line[i] != '='){
		i++;
	}
	end = i-1;
	length = end + 1 - start;
	var = (char*)malloc(length * sizeof(char) + 1);
	if(var == NULL){
		// -------------------------------------------
		//printf("var alloc fail");
		// -------------------------------------------
		return NULL;
	}
	memcpy(var,&line[start],length);
	var[length] = '\0';
	/*
	for(i = start; i <= end; i++){
		var[i - start] = line[i];
	}
	var[i] = '\0';//not sure for that line
	// ---------------------------
	//	printf("the var is %s \n",var);
	// ---------------------------
	*/
	return var;
}

bool checkTheVariableIsOk(char* var){
	//input: gets a string representing an arguement
	//output: check if it is one of the valid arguements. if yes return true. else return false.
	int i;
	for(i = 0; i < numberOfArguements; i++){
		if(strcmp(var,sysVarsControl[i])==0){
			// -------------------------
			//printf("found the var in the control \n");
			// -------------------------
			sysVarsFlag[i] = 1;
			currentVarIndex = i;
			return true;
		}
	}
	return false;
}

char* getTheValue(char* line){
	//input: a line in the config file, that was checked for being valid
	//output: get the value assigned in that line and return it as a string
	int i = 0, length, start;
	char* val;
	while(line[i] != '='){
		i++;
	}
	i++;
	while(line[i] == ' '){
		i++;
	}
	start = i;
	while(i <strlen(line) && line[i] != ' ' && line[i] != '\n'){
		i++;
	}
	length = i - start;
	val = (char*)malloc((sizeof(char)*length) + 1);
	memcpy(val, &line[start],length);
	val[length] = '\0';
	return val; 
	/*
	length = end + 1 - start;
	val = (char*)malloc(length * sizeof(char)+1);
	// ------------------------------------------------------
	// NOT SURE IF TO DELETE LATER
	if(val == NULL){
		printf("malloc fail\n");
		return NULL;
	}
	for(i = start; i <= end; i++){
		// -------------------------------------------------
		printf("%c \n", line[i]);
//		printf("flag\n");
//		printf("%d \n",(i-start));
		// --------------------------------------------------
		val[i-start] = line[i];
	}
	val[i] = '\0';
	
	// ----------------------------------
	//printf("the length is %d",(int)strlen(val));
	//printf("the value is %s", val);
	//printf("1%s1",val);
	// ----------------------------------
	return val;
	*/
}

SP_CONFIG_MSG checkTheValueIsOk(char* val){
	//input: value
	//output: SP_CONFIG_SUCCESS if the val holds the constraints
	int num;
	// ---------------------------------------------------------
	//printf("current var index %d \n", currentVarIndex); //DELETE LATER
	// ---------------------------------------------------------
	switch(currentVarIndex){
		case 0:
			return SP_CONFIG_SUCCESS; // I add those lines because the default return is false.
			break;
		case 1:
			return SP_CONFIG_SUCCESS;
			break;
		case 2:
			// ------------------------------------------------------------
			//printf("1%s1\n",val);
			//printf("the val is %s \n", val);
			//printf("the strcmp(.jpg,val) is %d \n", strcmp(val,".jpg"));
			// ------------------------------------------------------------
			if (strcmp(val,".jpg")!=0 && strcmp(val,".png")!=0 && strcmp(val,".bmp")!=0 && strcmp(val,".gif")!=0 ){
				// -----------------------------
				//printf("flagggggg"); // DELETE
				// -----------------------------
				return SP_CONFIG_INVALID_STRING;
			}
			else{
				return SP_CONFIG_SUCCESS;
			}
			break;
		case 3:
			return atoi(val) > 0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_INTEGER;
			break;
		case 4:
			num = atoi(val);	
			return num >= 10 && num <= 28 ? SP_CONFIG_SUCCESS : SP_CONFIG_INDEX_OUT_OF_RANGE; 
			break;
		case 5:
			return SP_CONFIG_SUCCESS;
			break;
		case 6:
			return atoi(val) > 0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_INTEGER;
			break;
		case 7:
			return strcmp(val,"true")==0 || strcmp(val,"false")==0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_ARGUMENT;
			break;
		case 8:
			return atoi(val) > 0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_INTEGER;
			break;
		case 9:
			return strcmp(val,"RANDOM")==0 || strcmp(val,"MAX_SPREAD")==0 || strcmp(val,"INCREMENTAL")==0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_ARGUMENT; 
			break;
		case 10:
			return atoi(val) > 0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_INTEGER;
			break;
		case 11:
			return strcmp(val,"true")==0 || strcmp(val,"false")==0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INVALID_ARGUMENT;
			break;
		case 12:
			return strcmp(val,"1")==0 || strcmp(val,"2")==0 || strcmp(val,"3")==0 || strcmp(val,"4")==0 ? SP_CONFIG_SUCCESS : SP_CONFIG_INDEX_OUT_OF_RANGE;
			break;
		case 13:
			return SP_CONFIG_SUCCESS;
			break;
	}
	return SP_CONFIG_INVALID_ARGUMENT;	
}

SP_CONFIG_MSG checkConstraints(char* line){
	//input: a line that is not a comment, not an empty line, and has valid spaces (according to formmat)
	//output: true if the constraints hold (the variable is valid, and the value holds the constraint)
	
	char* var;
	char* val;
	SP_CONFIG_MSG msg;
	
	var = getTheVariable(line);
	if(var == NULL){
		return SP_CONFIG_ALLOC_FAIL;
	}
	// ------------------------------
	//printf("the var is %s \n", var);
	// ------------------------------
	if(!checkTheVariableIsOk(var)){
		free(var);
		return SP_CONFIG_NOT_A_VARIABLE;
	}
	val = getTheValue(line);
	if(val == NULL){
		return SP_CONFIG_ALLOC_FAIL;
	}
	msg = checkTheValueIsOk(val);
	free(val);
	free(var);
	return msg;
}

SP_CONFIG_MSG checkLineIsValid(const char* filename, char* line, int lineNumber){
	//input: line in the config file
	//output: appropriate enum msg. also prints an error if there is one
	SP_CONFIG_MSG msg;
	if(isLineEmpty(line)){
		return SP_CONFIG_SUCCESS;
	}
	if(isLineComment(line)){
		return SP_CONFIG_SUCCESS;
	}	
	if(!spacesIsValid(line)){
		//guysFormatPrint(&filename,&(char)lineNumber,"Invalid configuration line",false);//check with guy what is the name of this function and change it!!!
		printf("File: %s\nLine: %d\nMessage: Invalid configuration line\n", filename,lineNumber);
		return SP_CONFIG_INVALID_ARGUMENT;
	}
	// -------------------------------
	//printf("flag1\n"); //DELETE
	//printf("the current var %d \n",currentVarIndex);
	// -------------------------------
	msg = checkConstraints(line);
	if(msg != SP_CONFIG_SUCCESS){
		// ------------------------------
		//printf("flag2");
		// ------------------------------ 
		printf("File: %s\nLine: %d\nMessage: Invalid value - constraints not met\n", filename, lineNumber);
		return msg;
	}
	return msg;
}

SP_CONFIG_MSG checkTheNecessaryVarsControl(const char* filename, int lineNumber){
	//input: config name and line number
	//output: appropriate enum msg. also prints an error if there is one
	// ------------------------------------------
	// delete this!
	//int i = 0;
	//printf("\n\n");
	//for(i = 0; i < 14; i++){
	//		printf("%d ", sysVarsFlag[i]);
	//}
	//printf("\n\n");
	// ------------------------------------------
	if(sysVarsFlag[0] == 0){
		printf("File: %s\nLine: %d\nMessage: Parameter spImagesDirectory is not set\n", filename,lineNumber);
		return SP_CONFIG_MISSING_DIR;
	}
	if(sysVarsFlag[1] == 0){
		printf("File: %s\nLine: %d\nMessage: Parameter spImagesPrefix is not set\n", filename,lineNumber);
		return SP_CONFIG_MISSING_PREFIX;
	}
	if(sysVarsFlag[2] == 0){
		printf("File: %s\nLine: %d\nMessage: Parameter spImagesSuffix is not set\n", filename,lineNumber);
		return SP_CONFIG_MISSING_SUFFIX;
	}
	if(sysVarsFlag[3] == 0){
		printf("File: %s\nLine: %d\nMessage: Parameter spNumOfImages is not set\n",filename,lineNumber);
		return SP_CONFIG_MISSING_NUM_IMAGES;
	}
	return SP_CONFIG_SUCCESS;
}

SP_CONFIG_MSG checkConfigFileValid(const char* filename){
	//input: gets a location of the file 
	//output: checks if the config file is valid and return appropriate enum msg
	char* line = (char*)malloc(sizeof(char)*1024);
	char* ptr = line;
	size_t len = 1024;
	ssize_t read;
	SP_CONFIG_MSG msg = SP_CONFIG_SUCCESS;
	int rowCount = 0;//need it to print the line where there is an error.
	
	if(filename == NULL){
		return SP_CONFIG_CANNOT_OPEN_FILE;
	}

	configFile = fopen(filename,"r");

	if(configFile == NULL){
		return SP_CONFIG_CANNOT_OPEN_FILE;
	}

	while((read = getline(&line,&len,configFile)) != -1){
		// ------------------------------------------------------------------------------
		//int m = 0;
		//for(m = 0; m<14; m++){
		//	printf(" %d ", sysVarsFlag[m]);
		//}
		//printf("\n");
		//printf("the line is %s", line);//note to self!!!: DELETE THIS!!!
		// -----------------------------------------------------------------------------
		rowCount++;//the first line in the config file is 1
		msg = checkLineIsValid(filename,line,rowCount);//also prints the regular msgs for the error. msgs number 1 or 2. (see the work instructions handed in class).
		if(msg != SP_CONFIG_SUCCESS){
			fclose(configFile);
			free(ptr);
			return msg;
		}
	}
	free(ptr);
	msg = checkTheNecessaryVarsControl(filename,rowCount);//also prints the regular msgs for the error. msgs 3 only. (see the work instructions handed in class).
	fclose(configFile);
	return msg;
}

void assignvariablesfromConfig(){
	publicConfig->spPCADimension = 20;
	publicConfig->spPCAFilename = "pca.yml";
	publicConfig->spNumOfFeatures = 100;
	publicConfig->spExtractionMode = true;
	publicConfig->spNumOfSimilarImages = 1;
	publicConfig->spKDTreeSplitMethod = MAX_SPREAD;
	publicConfig->spKNN = 1;
	publicConfig->spMinimalGUI = false;
	publicConfig->spLoggerLevel = 3;
	publicConfig->spLoggerFilename = "stdout";
	char* line = (char*)malloc(sizeof(char)*1024);
	char* ptr = line;
	char* val;
	char* var;
	int i = 0;
	int length;
	size_t len = 1024;
	ssize_t read;
	while((read = getline(&line,&len,configFile)) != -1){
		if(!isLineEmpty(line)){
			if(!isLineComment(line)){
				val = getTheValue(line);
				var = getTheVariable(line);
				for(i = 0; i < numberOfArguements; i++){
					if(strcmp(sysVarsControl[i],var) == 0){
						length = strlen(val);
						switch(i){
							case 0://need to check that the pointer is ok (the struct of SPConfig is pointer struct)
								publicConfig->spImagesDirectory = (char*)malloc((sizeof(char)*length) + 1);
								memcpy(publicConfig->spImagesDirectory, &val[0],length);
								publicConfig->spImagesDirectory[length] = '\0';
								break;
							case 1:
								publicConfig->spImagesPrefix = (char*)malloc((sizeof(char)*length) + 1);
								memcpy(publicConfig->spImagesPrefix, &val[0],length);
								publicConfig->spImagesPrefix[length] = '\0';
								break;
							case 2:
								publicConfig->spImagesSuffix = (char*)malloc((sizeof(char)*length) + 1);
								memcpy(publicConfig->spImagesSuffix, &val[0],length);
								publicConfig->spImagesSuffix[length] = '\0';
								break;
							case 3:
								publicConfig->spNumOfImages = atoi(val);
								break;
							case 4:
								publicConfig->spPCADimension = atoi(val);
								break;
							case 5:
								publicConfig->spPCAFilename = (char*)malloc((sizeof(char)*length) + 1);
								memcpy(publicConfig->spPCAFilename, &val[0],length);
								publicConfig->spPCAFilename[length] = '\0';
								break;
							case 6:
								publicConfig->spNumOfFeatures = atoi(val);
								break;
							case 7:
								publicConfig->spExtractionMode = strcmp(val,"true") == 0 ? true : false; //notice that it is case sensitive
								break;
							case 8:
								publicConfig->spNumOfSimilarImages = atoi(val);
								break;
							case 9:
								if(strcmp(val,"RANDOM") == 0){
									publicConfig->spKDTreeSplitMethod = RANDOM;
								}
								if(strcmp(val," MAX_SPREAD") == 0){
									publicConfig->spKDTreeSplitMethod =  MAX_SPREAD;
								}
								if(strcmp(val,"INCREMENTAL") == 0){
									publicConfig->spKDTreeSplitMethod = INCREMENTAL;
								}
								break;
							case 10:
								publicConfig->spKNN = atoi(val);
								break;
							case 11:
								publicConfig->spMinimalGUI = strcmp(val,"true") == 0 ? true : false;
								break;
							case 12:
								if(strcmp(val,"1") == 0){
									publicConfig->spLoggerLevel = SP_LOGGER_ERROR_LEVEL;
								}
								if(strcmp(val,"2") == 0){
									publicConfig->spLoggerLevel = SP_LOGGER_WARNING_ERROR_LEVEL;
								}
								if(strcmp(val,"3") == 0){
									publicConfig->spLoggerLevel = SP_LOGGER_INFO_WARNING_ERROR_LEVEL;
								}
								if(strcmp(val,"4") == 0){
									publicConfig->spLoggerLevel = SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL;
								}
								break;
							case 13:
								publicConfig->spLoggerFilename = (char*)malloc((sizeof(char)*length) + 1);
								memcpy(publicConfig->spLoggerFilename, &val[0],length);
								publicConfig->spLoggerFilename[length] = '\0';
								break;
						}
					}
				}
			}
		}
	}
	fclose(configFile);
	free(ptr);
	//printf("%d 2 %d 3", config.spNumOfImages, config.spNumOfFeatures);
}

/**
 * Creates a new system configuration struct. The configuration struct
 * is initialized based on the configuration file given by 'filename'.
 *
 * @param filename - the name of the configuration file
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return NULL in case an error occurs. Otherwise, a pointer to a struct which
 * 		   contains all system configuration.
 *
 * The resulting value stored in msg is as follow:
 * - SP_CONFIG_INVALID_ARGUMENT - if filename == NULL
 * - SP_CONFIG_CANNOT_OPEN_FILE - if the configuration file given by filename cannot be open
 * - SP_CONFIG_ALLOC_FAIL - if an allocation failure occurred
 * - SP_CONFIG_INVALID_INTEGER - if a line in the config file contains invalid integer
 * - SP_CONFIG_INVALID_STRING - if a line in the config file contains invalid string
 * - SP_CONFIG_MISSING_DIR - if spImagesDirectory is missing
 * - SP_CONFIG_MISSING_PREFIX - if spImagesPrefix is missing
 * - SP_CONFIG_MISSING_SUFFIX - if spImagesSuffix is missing
 * - SP_CONFIG_MISSING_NUM_IMAGES - if spNumOfImages is missing
 * - SP_CONFIG_SUCCESS - in case of success
 *
 *
 */
SPConfig* spConfigCreate(const char* filename, SP_CONFIG_MSG* msg){
	assert (msg); 
	//malloc
	publicConfig = (SPConfig*)malloc(sizeof(SPConfig));
	if (!publicConfig){
		*msg = SP_CONFIG_ALLOC_FAIL;
		return NULL;
	}

	//verify there's a file name
	if (!filename){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		spConfigDestroy(*publicConfig);
		return NULL;
	}

	//verify file opens
	configFile = fopen(filename, "r");
	if (!configFile){
		*msg = SP_CONFIG_CANNOT_OPEN_FILE;
		spConfigDestroy(*publicConfig);
		return NULL;
	}

	//verify file is valid
	*msg = checkConfigFileValid(filename);
	if (*msg != SP_CONFIG_SUCCESS){
		spConfigDestroy(*publicConfig);
		return NULL;
	}

	//if all works - get the variables
	assignvariablesfromConfig(); //Viktor: including msg (ex. info msg "used default value "20" for spPCADimension")

	if(configFile)
		fclose(configFile);

	return publicConfig;
}


/*
 * Returns true if spExtractionMode = true, false otherwise.
 * Assuming there is a valid value of ExtractionMode
 *
 * @param config - the configuration structure
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return true if spExtractionMode = true, false otherwise.
 *
 * - SP_CONFIG_INVALID_ARGUMENT - if config == NULL
 * - SP_CONFIG_SUCCESS - in case of success
 */
bool spConfigIsExtractionMode(const SPConfig config, SP_CONFIG_MSG* msg){
	assert(*msg); //err
	*msg = SP_CONFIG_SUCCESS;
	if (!publicConfig){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		return false;
	}
	return publicConfig->spExtractionMode;
}


/*
 * Returns true if spMinimalGUI = true, false otherwise.
 * Assuming there is a valid value of MinimalGui
 *
 * @param config - the configuration structure
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return true if spExtractionMode = true, false otherwise.
 *
 * - SP_CONFIG_INVALID_ARGUMENT - if config == NULL
 * - SP_CONFIG_SUCCESS - in case of success
 */
bool spConfigMinimalGui(const SPConfig config, SP_CONFIG_MSG* msg){
	assert(*msg); //err
	*msg = SP_CONFIG_SUCCESS;
	if (!publicConfig){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		return false;
	}
	return publicConfig->spMinimalGUI;
}


/*
 * Returns the number of images set in the configuration file, i.e the value
 * of spNumOfImages.
 * Assuming there is a valid value of NumOfImages
 *
 * @param config - the configuration structure
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return positive integer in success, negative integer otherwise.
 *
 * - SP_CONFIG_INVALID_ARGUMENT - if config == NULL
 * - SP_CONFIG_SUCCESS - in case of success
 */
int spConfigGetNumOfImages(const SPConfig config, SP_CONFIG_MSG* msg){
	assert(*msg); //err
	*msg = SP_CONFIG_SUCCESS;
	if (!publicConfig){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		return -1;
	}
	//printf("%d 1", config->spNumOfImages);
	return publicConfig->spNumOfImages;
}


/*
 * Returns the number of features to be extracted. i.e the value
 * of spNumOfFeatures.
 * Assuming there is a valid value of NumOfFeatures
 *
 * @param config - the configuration structure
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return positive integer in success, negative integer otherwise.
 *
 * - SP_CONFIG_INVALID_ARGUMENT - if config == NULL
 * - SP_CONFIG_SUCCESS - in case of success
 */
int spConfigGetNumOfFeatures(const SPConfig config, SP_CONFIG_MSG* msg){
	assert(*msg); //err
	*msg = SP_CONFIG_SUCCESS;
	if (!publicConfig){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		return -1;
	}
	return publicConfig->spNumOfImages;
}


/**
 * Returns the dimension of the PCA. i.e the value of spPCADimension.
 * Assuming there is a valid value of PCADim
 *
 * @param config - the configuration structure
 * @assert msg != NULL
 * @param msg - pointer in which the msg returned by the function is stored
 * @return positive integer in success, negative integer otherwise.
 *
 * - SP_CONFIG_INVALID_ARGUMENT - if config == NULL
 * - SP_CONFIG_SUCCESS - in case of success
 */
int spConfigGetPCADim(const SPConfig config, SP_CONFIG_MSG* msg){
	assert(*msg); //err
	*msg = SP_CONFIG_SUCCESS;
	if (!publicConfig){
		*msg = SP_CONFIG_INVALID_ARGUMENT;
		return -1;
	}
	return publicConfig->spPCADimension;
}


/**
 * Given an index 'index' the function stores in imagePath the full path of the
 * ith image.
 * assuming there is a valid value of Directory, Prefix, Suffix, N. Images, index
 * assuming that the total length in shorter than MAX_LEN (1025)
 * assuming imagePath's can contain MAX_LEN (1025)
 * assuming the resulting path will be checked as valid before use
 * not assuming imagePath is an empty string
 *
 * For example:
 * Given that the value of:
 *  spImagesDirectory = "./images/"
 *  spImagesPrefix = "img"
 *  spImagesSuffix = ".png"
 *  spNumOfImages = 17
 *  index = 10
 *
 * The functions stores "./images/img10.png" to the address given by imagePath.
 * Thus the address given by imagePath must contain enough space to
 * store the resulting string.
 *
 * @param imagePath - an address to store the result in, it must contain enough space.
 * @param config - the configuration structure
 * @param index - the index of the image.
 *
 * @return
 * - SP_CONFIG_INVALID_ARGUMENT - if imagePath == NULL or config == NULL
 * - SP_CONFIG_INDEX_OUT_OF_RANGE - if index >= spNumOfImages
 * - SP_CONFIG_SUCCESS - in case of success
 */
SP_CONFIG_MSG spConfigGetImagePath(char* imagePath, const SPConfig config,
		int index){
	if (!imagePath || publicConfig)
		return SP_CONFIG_INVALID_ARGUMENT;
	if (index >= publicConfig->spNumOfImages)
		return SP_CONFIG_INDEX_OUT_OF_RANGE;
	char Sindex[MAX_LEN];
	sprintf(Sindex, "%d", index);
	strcpy(imagePath, "");// in case it wasn't empty
	strcat(imagePath, publicConfig->spImagesDirectory);
	strcat(imagePath, publicConfig->spImagesPrefix);
	strcat(imagePath, Sindex);
	strcat(imagePath, publicConfig->spImagesSuffix);
	return SP_CONFIG_SUCCESS;
}


/**
 * The function stores in pcaPath the full path of the pca file.
 * assuming there is a valid value of Directory and File Name
 * assuming that the total length in shorter than MAX_LEN (1025)
 * assuming pcaPath's can contain MAX_LEN (1025)
 * assuming the resulting path will be checked as valid before use
 * not assuming imagePath is an empty string
 *
 * For example given the values of:
 *  spImagesDirectory = "./images/"
 *  spPcaFilename = "pca.yml"
 *
 * The functions stores "./images/pca.yml" to the address given by pcaPath.
 * Thus the address given by pcaPath must contain enough space to
 * store the resulting string.
 *
 * @param imagePath - an address to store the result in, it must contain enough space.
 * @param config - the configuration structure
 * @return
 *  - SP_CONFIG_INVALID_ARGUMENT - if imagePath == NULL or config == NULL
 *  - SP_CONFIG_SUCCESS - in case of success
 */
SP_CONFIG_MSG spConfigGetPCAPath(char* pcaPath, const SPConfig config){
	if (!pcaPath) //|| !config*/
		return SP_CONFIG_INVALID_ARGUMENT;
	strcpy(pcaPath, "");// in case it wasn't empty
	strcat(pcaPath, publicConfig->spImagesDirectory);
	strcat(pcaPath, publicConfig->spPCAFilename);
	return SP_CONFIG_SUCCESS;
}


/**
 * Frees all memory resources associate with config.
 * If config == NULL nothig is done.
 */
void spConfigDestroy(SPConfig* config){
	if (publicConfig){
		if (publicConfig->spImagesDirectory)
			free(publicConfig->spImagesDirectory);
		if (publicConfig->spImagesPrefix)
			free(publicConfig->spImagesPrefix);
		if (publicConfig->spImagesSuffix)
			free(publicConfig->spImagesSuffix);
		if (publicConfig->spPCAFilename)
			free(publicConfig->spPCAFilename);
		if (publicConfig->spLoggerFilename)
			free(publicConfig->spLoggerFilename);
		free(publicConfig);
	}
	if (configFile)
		fclose(configFile);
}

char* getSPLoggerFilename(){
		return publicConfig->spLoggerFilename;
}

SP_LOGGER_LEVEL getSPLoggerLevel(){
	return publicConfig->spLoggerLevel;
}

bool getSPExtractionMode(){
	return publicConfig->spExtractionMode;
}

int getSPNumOfImages(){
	return publicConfig->spNumOfImages;
}

int getSpKNN(){
	return publicConfig->spKNN;
}

bool getSPMinimalGUI(){
	return publicConfig->spMinimalGUI;
}

bool getSPNumOfSimilarImages(){
	return publicConfig->spNumOfSimilarImages;
}

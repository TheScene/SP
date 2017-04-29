#include "SPLogger.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#define SP_LOGGER_OPEN_MODE "w" //File open mode
#define MAX_LEN 1025 //1024 string including the terminating null character


/** A type used for defining the logger**/
typedef struct sp_logger_t{
	FILE* outputChannel; //The logger file
	bool isStdOut; //Indicates if the logger is stdout
	SP_LOGGER_LEVEL level; //Indicates the level
}*SPLogger;


/**
 * Creates a logger. This function should be called once, prior
 * to the usage of any SP Logger print functions. It is the responsibility
 * of the user to create the logger prior to usage, and the logger
 * must be destroyed at the end of usage.
 *
 * @param filename - The name of the log file, if not specified stdout is used
 * 					 as default.
 * @param level - The level of the logger prints
 * @return
 * SP_LOGGER_DEFINED 			- The logger has been defined
 * SP_LOGGER_OUT_OF_MEMORY 		- In case of memory allocation failure
 * SP_LOGGER_CANNOT_OPEN_FILE 	- If the file given by filename cannot be opened
 * SP_LOGGER_SUCCESS 			- In case the logger has been successfully opened
 */
SP_LOGGER_MSG spLoggerCreate(const char* filename, SP_LOGGER_LEVEL level) {
	if (logger != NULL) { //Already defined
		return SP_LOGGER_DEFINED;
	}
	logger = (SPLogger) malloc(sizeof(*logger));
	if (logger == NULL) { //Allocation failure
		return SP_LOGGER_OUT_OF_MEMORY;
	}
	logger->level = level; //Set the level of the logger
	if (filename == NULL) { //In case the filename is not set use stdout
		logger->outputChannel = stdout;
		logger->isStdOut = true;
	} else { //Otherwise open the file in write mode
		logger->outputChannel = fopen(filename, SP_LOGGER_OPEN_MODE);
		if (logger->outputChannel == NULL) { //Open failed
			spLoggerDestroy(logger);
			logger = NULL;
			return SP_LOGGER_CANNOT_OPEN_FILE;
		}
		logger->isStdOut = false;
	}
	return SP_LOGGER_SUCCESS;
}


/**
 * Frees all memory allocated for the logger. If the logger is not defined
 * then nothing happens.
 */
void spLoggerDestroy() {
	if (!logger) {
		return;
	}
	if (!logger->isStdOut) {//Close file only if not stdout
		fclose(logger->outputChannel);
	}
	free(logger);//free allocation
	logger = NULL;
}


/*
 * prints the msg by this format:
 * ---<msg type>---
 *  - file: <file>
 *  - function: <function>
 *  - line: <line>
 *  - message: <msg>
 *
 *   <msg Type> - ERROR, WARNING, INFO, DEBUG
 * 	<file> 	   - is the string given by file, it represents the file in which
 * 		   		 the error print call occurred.
 * 	<function> - is the string given by function, it represents the function in which
 * 			   	 the error print call occurred.
 * 	<line> 	   - is the string given by line, it represents the line in which
 * 		   		 the error print call occurred
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	returns the final fomated msg including /n between the lines and at the end (always needed)
 */
char* formatMSG(const char* msg, const char* file,	const char* function, const int line, char* msgType){
	char* formatedMSG = calloc(1, MAX_LEN*5 + strlen(msg));
	char* temp = calloc(1, MAX_LEN*5 + strlen(msg));
	if (!formatedMSG || !temp)
		spLoggerPrintMsg("Allocation error");
	sprintf(temp, "---%s---\n", msgType);
	strcat(formatedMSG, temp);
	sprintf(temp, "- file: %s\n", file);
	strcat(formatedMSG, temp);
	sprintf(temp, "- function: %s\n", function);
	strcat(formatedMSG, temp);
	sprintf(temp, "- line: %d\n", line);
	strcat(formatedMSG, temp);
	sprintf(temp, "- message: %s\n", msg);
	strcat(formatedMSG, temp);
	free(temp);
	return formatedMSG;
}


/*
 * prints the msg by this format:
 * ---INFO---
 *  - message: <msg>
 *
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	returns the final fomated msg including \n between the lines and at the end (always needed)
 */
char* formatInfoMSG(const char* msg){
	char* formatedMSG = calloc(1, MAX_LEN*5 + strlen(msg));
	char* temp = calloc(1, MAX_LEN*5 + strlen(msg));
	if (!formatedMSG || !temp)
		spLoggerPrintMsg("Allocation error");
	sprintf(formatedMSG, "---INFO---\n");
	sprintf(temp, "- message: %s\n", msg);
	strcat(formatedMSG, temp);
	free(temp);
	return formatedMSG;
}


/**
 * 	Prints error message. The error message format is given below:
 * 	---ERROR---
 * 	- file: <file>
 *  - function: <function>
 *  - line: <line>
 *  - message: <msg>
 *
 * 	<file> 	   - is the string given by file, it represents the file in which
 * 		   		 the error print call occurred.
 * 	<function> - is the string given by function, it represents the function in which
 * 			   	 the error print call occurred.
 * 	<line> 	   - is the string given by line, it represents the line in which
 * 		   		 the error print call occurred
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	Error messages will be printed at levels:
 *
 * 	SP_LOGGER_ERROR_LEVEL,
 *	SP_LOGGER_WARNING_ERROR_LEVEL,
 *	SP_LOGGER_INFO_WARNING_ERROR_LEVEL,
 *	SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL
 *
 * 	A new line will be printed after the print call.
 *
 * @param msg     	- The message to printed
 * @param file    	- A string representing the filename in which spLoggerPrintError call occurred
 * @param function 	- A string representing the function name in which spLoggerPrintError call ocurred
 * @param line		- A string representing the line in which the function call occurred
 * @return
 * SP_LOGGER_UNDIFINED 			- If the logger is undefined
 * SP_LOGGER_INVAlID_ARGUMENT	- If any of msg or file or function are null or line is negative
 * SP_LOGGER_WRITE_FAIL			- If Write failure occurred
 * SP_LOGGER_SUCCESS			- otherwise
 */
SP_LOGGER_MSG spLoggerPrintError(const char* msg, const char* file,
		const char* function, const int line){
	int success;
	char* formatedMSG;
	if(!logger) // if logger is NULL
		return SP_LOGGER_UNDIFINED;
	// no need to check the logger.SP_LOGGER_LEVEL in this case, but in the next functions it's needed
	if (!msg || !file || !function || line < 0) // if msg / file / func / line invalid
		return SP_LOGGER_INVAlID_ARGUMENT;
	formatedMSG = formatMSG(msg, file, function, line, "ERROR"); // get the msg to be printed
	if (logger->isStdOut == true)
		success = printf("%s", formatedMSG);
	else
		success = fprintf(logger->outputChannel, "%s", formatedMSG); // print to stdout/file
	free(formatedMSG);
	if(success < 0)
		return SP_LOGGER_WRITE_FAIL;
	else
		return SP_LOGGER_SUCCESS;
}


/**
 * 	Prints warning message. The warning message format is given below:
 * 	---WARNING---
 * 	- file: <file>
 *  - function: <function>
 *  - line: <line>
 *  - message: <msg>
 *
 * 	<file> 	   - is the string given by file, it represents the file in which
 * 		   		 the warning print call occurred.
 * 	<function> - is the string given by function, it represents the function in which
 * 			   	 the warning print call occurred.
 * 	<line> 	   - is the string given by line, it represents the line in which
 * 		   		 the warning print call occurred
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	Warning messages will be printed at levels:
 *
 *	SP_LOGGER_WARNING_ERROR_LEVEL,
 *	SP_LOGGER_INFO_WARNING_ERROR_LEVEL,
 *	SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL
 *
 *	A new line will be printed after the print call.
 *
 * @param msg     	- The message to printed
 * @param file    	- A string representing the filename in which spLoggerPrintWarning call occurred
 * @param function 	- A string representing the function name in which spLoggerPrintWarning call ocurred
 * @param line		- A string representing the line in which the spLoggerPrintWarning call occurred
 * @return
 * SP_LOGGER_UNDIFINED 			- If the logger is undefined
 * SP_LOGGER_INVAlID_ARGUMENT	- If any of msg or file or function are null or line is negative
 * SP_LOGGER_WRITE_FAIL			- If write failure occurred
 * SP_LOGGER_SUCCESS			- otherwise
 */
SP_LOGGER_MSG spLoggerPrintWarning(const char* msg, const char* file,
		const char* function, const int line){
	char* formatedMSG;
	if(!logger) // if logger is NULL
		return SP_LOGGER_UNDIFINED;
	if (logger->level == SP_LOGGER_ERROR_LEVEL)
		return SP_LOGGER_SUCCESS;
	if (!msg || !file || !function || line < 0) // if msg / file / func / line invalid
		return SP_LOGGER_INVAlID_ARGUMENT;
	formatedMSG = formatMSG(msg, file, function, line, "WARNING"); // get the msg to be printed
	if (logger->isStdOut == true)
		printf("%s", formatedMSG);
	else
		fprintf(logger->outputChannel, "%s", formatedMSG); // print to stdout/file
	free(formatedMSG);
	return SP_LOGGER_SUCCESS;
}


/**
 * 	Prints Info message. The info message format is given below:
 * 	---INFO---
 *  - message: <msg>
 *
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	Info messages will be printed at levels:
 *
 *	SP_LOGGER_INFO_WARNING_ERROR_LEVEL,
 *	SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL
 *
 * 	A new line will be printed after the print call.
 *
 * @param msg     	- The message to printed
 * @return
 * SP_LOGGER_UNDIFINED 			- If the logger is undefined
 * SP_LOGGER_INVAlID_ARGUMENT	- If msg is null
 * SP_LOGGER_WRITE_FAIL			- If Write failure occurred
 * SP_LOGGER_SUCCESS			- otherwise
 */
SP_LOGGER_MSG spLoggerPrintInfo(const char* msg){
	int success;
	char* formatedMSG;
	if(!logger) // if logger is NULL
		return SP_LOGGER_UNDIFINED;
	if (logger->level == SP_LOGGER_ERROR_LEVEL || logger->level == SP_LOGGER_WARNING_ERROR_LEVEL)
		return SP_LOGGER_SUCCESS;
	if (!msg) // if msg / file / func / line invalid
		return SP_LOGGER_INVAlID_ARGUMENT;
	formatedMSG = formatInfoMSG(msg); // get the msg to be printed
	if (logger->isStdOut == true)
		success = printf("%s", formatedMSG);
	else
		success = fprintf(logger->outputChannel, "%s", formatedMSG); // print to stdout/file
	free(formatedMSG);
	if(success < 0)
		return SP_LOGGER_WRITE_FAIL;
	else
		return SP_LOGGER_SUCCESS;
}


/**
 * 	Prints the debug message. The debug message format is given below:
 * 	---DEBUG---
 * 	- file: <file>
 *  - function: <function>
 *  - line: <line>
 *  - message: <msg>
 *
 * 	<file> 	   - is the string given by file, it represents the file in which
 * 		   		 the debug print call occurred.
 * 	<function> - is the string given by function, it represents the function in which
 * 			   	 the debug print call occurred.
 * 	<line> 	   - is the string given by line, it represents the line in which
 * 		   		 the debug print call occurred
 * 	<msg> 	   - is the string given by msg, it contains the msg given by the user
 *
 * 	Debug messages will be printed at level:
 *
 *	SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL
 *
 * 	A new line will be printed after the print call.
 *
 * @param msg     	- The message to printed
 * @param file    	- A string representing the filename in which spLoggerPrintWarning call occurred
 * @param function 	- A string representing the function name in which spLoggerPrintWarning call ocurred
 * @param line		- A string representing the line in which the function call occurred
 * @return
 * SP_LOGGER_UNDIFINED 			- If the logger is undefined
 * SP_LOGGER_INVAlID_ARGUMENT	- If any of msg or file or function are null or line is negative
 * SP_LOGGER_WRITE_FAIL			- If Write failure occurred
 * SP_LOGGER_SUCCESS			- otherwise
 */
SP_LOGGER_MSG spLoggerPrintDebug(const char* msg, const char* file,
		const char* function, const int line){
	int success;
	char* formatedMSG;
	if(!logger) // if logger is NULL
		return SP_LOGGER_UNDIFINED;
	if (logger->level != SP_LOGGER_DEBUG_INFO_WARNING_ERROR_LEVEL)
		return SP_LOGGER_SUCCESS;
	if (!msg || !file || !function || line < 0) // if msg / file / func / line invalid
		return SP_LOGGER_INVAlID_ARGUMENT;
	formatedMSG = formatMSG(msg, file, function, line, "DEBUG"); // get the msg to be printed
	if (logger->isStdOut == true)
		success = printf("%s", formatedMSG);
	else
		success = fprintf(logger->outputChannel, "%s", formatedMSG); // print to stdout/file
	free(formatedMSG);
	if(success < 0)
		return SP_LOGGER_WRITE_FAIL;
	else
		return SP_LOGGER_SUCCESS;
}


/**
 * The given message is printed. A new line is printed at the end of msg
 * The message will be printed in all levels.
 *
 * @param msg - The message to be printed
 * @return
 * SP_LOGGER_UNDIFINED 			- If the logger is undefined
 * SP_LOGGER_INVAlID_ARGUMENT	- If msg is null
 * SP_LOGGER_WRITE_FAIL			- If Write failure occurred
 * SP_LOGGER_SUCCESS			- otherwise
 */
SP_LOGGER_MSG spLoggerPrintMsg(const char* msg){
	int success;
	if(!logger) // if logger is NULL
		return SP_LOGGER_UNDIFINED;
	if (!msg) // if msg invalid
		return SP_LOGGER_INVAlID_ARGUMENT;
	success = printf("%s\n", msg);
	if(success < 0)
		return SP_LOGGER_WRITE_FAIL;
	else
		return SP_LOGGER_SUCCESS;
}

/*
 * main_aux.h
 *
 *  Created on: Mar 30, 2017
 *      Author: gleit
 */

#ifndef MAIN_AUX_H_
#define MAIN_AUX_H_
#include <stdbool.h>
#include "unit_tests\\unit_test_util.h"
#include "SPConfig.h"
#include "SPLogger.h"

void runTests();

char* getConfigFileName(int argc, char** agrv);

bool manageCMSG(SP_CONFIG_MSG* Cmsg);

bool manageLMSG(SP_LOGGER_MSG* Smsg);

void freeAll(SPConfig config, SPLogger logger);


/*
 * this function is activated when we are in extraction mode
 * it takes all the pictures, and extracts the features
 * into the appropriate files. img1.jpg feats to feats1
 * img2.jpg to feats2 and so on.
 *
 * @param -
 *
 * @return - void. creates feats files, and extract to them
 * the featurs.
 *
 * */
int spExtract();

#endif /* MAIN_AUX_H_ */

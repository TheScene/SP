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
#include "KDTree.h"

void runTests();

char* getConfigFileName(int argc, char** agrv);

bool manageCMSG(SP_CONFIG_MSG* Cmsg);

bool manageLMSG(SP_LOGGER_MSG* Smsg);

void freeAll(SPLogger logger);


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

/*
 * recives the command from the user.
 * this is the user interface function
 *
 * @param - none
 * @return - an integer. 0 if exit, 1 if searching
 * for an image
 * */
int receiveCommand();

/*
 * recieves a query image address from the user
 *
 * @parm - a place to hold the string
 * @return- the place to hold the string is updated and 1 if exit. 0 otherwise
 *
 * */
int receiveQuery(char* buffer);

/*
 * recives the kdTree and a picture, and finds the
 * SPKNN images the are the most similar to it.
 *
 * @param - a string representing thpath of the query image
 * and the root of the kdtree. which is curr.
 *
 * @return - integer array with the size of the number of images
 * each cell in the array i has the value of how close the image i
 * to the query image.
 * */
int* findSimilarImages(KDTreeNode* curr,char* query);


/*
 * this function shows the images. with pictures or addresses
 *
 * @param - showOrNot bool that is true when we want to show pictures
 * and false when we want to show addresses. imgSimilarityRankArr is an
 * integer array that has the rank of how each photo is similar to the
 * query image (see findSimilarImages())
 * */
void showImages(bool showOrNot,int* imgSimilarityRankArr, char* queryPath);

#endif  MAIN_AUX_H_

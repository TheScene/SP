/*
 * SPPoint.c
 *
 *  Created on: Dec 16, 2016
 *      Author: Guy Leitersdorf
 */
#include "SPPoint.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

SPPoint* spPointCreate(double* data, int dim, int index){
	if (dim<=0 || index < 0)
		return NULL;
	if(data == NULL){
		return NULL;
	}
	SPPoint* newPoint = (SPPoint*)malloc(sizeof(SPPoint));
	newPoint->data = (double*)malloc(dim*sizeof(double));
	if (newPoint == NULL)
		return NULL;
	if (newPoint->data == NULL)
		return NULL;
	for (int i=0; i<dim; i++){
			newPoint->data[i] = data[i];
		}
	newPoint->dim = dim;
	newPoint->index = index;
	return newPoint;
}



SPPoint* spPointCopy(SPPoint* source){
	assert(source != NULL);
	return spPointCreate(source->data, source->dim, source->index);
	}


void spPointDestroy(SPPoint* point){
	if (point != NULL){
		free(point->data);
		free(point);
	}
	return;
}

int spPointGetDimension(SPPoint* point){
	assert(point!=NULL);
	return point->dim;
}

int spPointGetIndex(SPPoint* point){
	assert(point!=NULL);
	return point->index;
}

double spPointGetAxisCoor(SPPoint* point, int axis){
	assert(point!=NULL && axis < point->dim);
	return point->data[axis];
}

double spPointL2SquaredDistance(SPPoint* p, SPPoint* q){
	assert(p!=NULL && q!=NULL && p->dim == q->dim);
	double sum = 0;
	double dif;
	for (int i = 0; i < p->dim; i++) {
		dif = spPointGetAxisCoor(p,i) - spPointGetAxisCoor(q,i);
		sum += dif*dif;
	}
	return sum;
}


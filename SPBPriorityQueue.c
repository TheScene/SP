/*
 * SPBPriorityQueue.c
 *
 *  Created on: Dec 20, 2016
 *      Author: Guy Leitersdorf
 */
#include "SPBPriorityQueue.h"
#include "SPPoint.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
//Element functions

/*
 * element create helper function
 * */
BPQueueElement* spBPElementcreate(int index, double value){
	BPQueueElement* element = (BPQueueElement*)malloc(sizeof(BPQueueElement)); //outofmem?
	element->index = index;
	element->value = value;
	return element;
}

/*
 * element copy helper function
 * */
BPQueueElement* spBPElementCopy(BPQueueElement* source){
	return spBPElementcreate(source->index, source->value);
}

/*
 * element destroy helper function
 * */
void spBPElementDestroy(BPQueueElement* source) {
	free(source);
}

// End of Element functions









// Queue functions

struct sp_bp_queue_t{
	int first;
	int size;
	int maxSize;
	BPQueueElement** array;
};


/*
 * helper function that tells us if
 * source queue is null
 * */
bool spBPQueueIsNULL (SPBPQueue* source){
	return (source == NULL);
}


SPBPQueue* spBPQueueCreate(int maxSize){
	assert(maxSize>0);
	SPBPQueue* newQueue = (SPBPQueue*)malloc(sizeof(SPBPQueue));
	newQueue->size=0;
	newQueue->maxSize=maxSize; //assert maxsize >0
	newQueue->array = (BPQueueElement**)malloc(sizeof(BPQueueElement**) * maxSize);
	return newQueue;
}


void spBPQueueDestroy(SPBPQueue* source){
	spBPQueueClear(source);
	free (source);
}


void spBPQueueClear(SPBPQueue* source){
	if (source == NULL)
		return;
	for (int i = 0; i < source->maxSize; i++) {
		spBPElementDestroy(source->array[i]);
	}


	source->size = 0;
}


int spBPQueueSize(SPBPQueue* source){
	assert(source!=NULL);
	return source->size;
}


int spBPQueueGetMaxSize(SPBPQueue* source){
	assert(source!=NULL);
	return source->maxSize;
}

/*
 * get the index in the array that holds the last
 * element in the queue
 *
 * @param source
 * @assert source!=null
 * @return the place of last element in array
 * */
int spBPQueueGetLast(SPBPQueue* source){
	assert(source!=NULL);
	return ((source->first + source->size - 1)%source->maxSize);
}


/*
 * finds the right place in the queue and
 * inserts the element there
 * @param source queue and an element
 * @assert source != null && element != null
 * @return inserts the element de facto
 */
void spBPQueueEnqueueLegal(SPBPQueue* source, BPQueueElement* element){
	assert(source!=NULL && element!=NULL);
	int last = spBPQueueGetLast(source);
	int i = last;
	source->array[i] = element;
	for (int j=0 ; j<source->size-1; j++) {
			if (source->array[i]->value  >  source->array[(i-1)%(source->maxSize)]->value){
				return;
			}
			if((source->array[i]->value == source->array[(i-1)%(source->maxSize)]->value) && (source->array[i]->index  >  source->array[(i-1)%(source->maxSize)]->index)){
				return;
			}
			source->array[i] = source->array[(i-1)%(source->maxSize)];
			source->array[(i-1)%(source->maxSize)] = element;
			i = (i-1)%source->maxSize;
		}
}


SP_BPQUEUE_MSG spBPQueueEnqueue(SPBPQueue* source, int index, double value){
	assert(source!=NULL);
	BPQueueElement* element = spBPElementcreate(index, value);
	if (element == NULL) //case malloc failed
		return SP_BPQUEUE_OUT_OF_MEMORY;
	if(spBPQueueIsEmpty(source)){ //case list empty
		source->first = 0;
		source->size++;
		source->array[0] = element;
		return SP_BPQUEUE_SUCCESS;
	}
	BPQueueElement* last = source->array[spBPQueueGetLast(source)];
	if (spBPQueueIsFull(source)){//case element rejected
		if(last->value < value){
			return SP_BPQUEUE_FULL;
		}
		if (last->value == value && last->index < index){
			return SP_BPQUEUE_FULL;
		}
	}
	if (spBPQueueIsFull(source)) //case list too big
		spBPElementDestroy(last);
	else{
		source->size++;
	}
	spBPQueueEnqueueLegal(source, element);// finally - enqueue
	return SP_BPQUEUE_SUCCESS;
}



SP_BPQUEUE_MSG spBPQueueDequeue(SPBPQueue* source){
	if (spBPQueueIsNULL(source))
		return SP_BPQUEUE_INVALID_ARGUMENT;
	if (spBPQueueIsEmpty(source))
		return SP_BPQUEUE_EMPTY;
	spBPElementDestroy(source->array[source->first]);
	source->first = ((source->first +1)%source->maxSize);
	source->size--;
	return SP_BPQUEUE_SUCCESS;
}


/*
 * helper function to return the appropriate
 * msg for the peek function
 *
 * @param source queue and element res
 * @return the appropriate msg
 * */
SP_BPQUEUE_MSG spBPQueuePeekGetMSG(SPBPQueue* source, BPQueueElement* res){
	if (res == NULL || spBPQueueIsNULL(source))
		return SP_BPQUEUE_INVALID_ARGUMENT;
	if (spBPQueueIsEmpty(source))
		return SP_BPQUEUE_EMPTY;
	return SP_BPQUEUE_SUCCESS;
}


SP_BPQUEUE_MSG spBPQueuePeek(SPBPQueue* source, BPQueueElement* res){
	SP_BPQUEUE_MSG msg = spBPQueuePeekGetMSG(source, res);
	if (msg == SP_BPQUEUE_SUCCESS){
		res->value = source->array[source->first]->value;
		res->index = source->array[source->first]->index;
	}
	return msg;
}


SP_BPQUEUE_MSG spBPQueuePeekLast(SPBPQueue* source, BPQueueElement* res){
	SP_BPQUEUE_MSG msg = spBPQueuePeekGetMSG(source, res);
	if (msg == SP_BPQUEUE_SUCCESS){
		res->value = source->array[spBPQueueGetLast(source)]->value;
		res->index = source->array[spBPQueueGetLast(source)]->index;
	}
	return msg;
}


SPBPQueue* spBPQueueCopy(SPBPQueue* source){
	if (spBPQueueIsNULL(source))
		return NULL;
	SPBPQueue* newQueue = spBPQueueCreate(source->maxSize);
	for (int i = source->first; i < source->size; i=((i+1)%source->maxSize))
		newQueue->array[i]=source->array[i];
	newQueue->first = source->first;
	newQueue->size = source->size;
	return newQueue;
}


double spBPQueueMinValue(SPBPQueue* source){
	assert(source!=NULL && !spBPQueueIsEmpty(source));
	return source->array[source->first]->value;
}


double spBPQueueMaxValue(SPBPQueue* source){
	assert(source!=NULL && !spBPQueueIsEmpty(source));
	return source->array[spBPQueueGetLast(source)]->value;
}


bool spBPQueueIsEmpty(SPBPQueue* source){
	assert(source!=NULL);
	return (source->size == 0);
}


bool spBPQueueIsFull(SPBPQueue* source){
	assert(source!=NULL);
	if (source->size == source->maxSize)
		return true;
	return false;
}


// End of Queue functions

// EOF

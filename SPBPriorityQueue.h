#ifndef SPBPRIORITYQUEUE_H_
#define SPBPRIORITYQUEUE_H_
#include <stdbool.h>

/**
 * SP Bounded Priority Queue summary
 *
 * the Bounded Priority Queue that we implemented
 * uses an array, to keep all the elements.
 * we are using the array in a circular manner,
 * that is, when we Dequeue an element from the
 * array, we just take out and delete the most
 * left element in the array, and when the array
 * has no longer place in the right, we simply
 * put elements in the most left side of the array
 * (given that the most left side is empty).
 */


/** type used to define Bounded priority queue **/
typedef struct sp_bp_queue_t SPBPQueue;

typedef struct sp_bpq_element_t {
	int index;
	double value;
} BPQueueElement;

/** type for error reporting **/
typedef enum sp_bp_queue_msg_t {
	SP_BPQUEUE_OUT_OF_MEMORY,
	SP_BPQUEUE_FULL,
	SP_BPQUEUE_EMPTY,
	SP_BPQUEUE_INVALID_ARGUMENT,
	SP_BPQUEUE_SUCCESS
} SP_BPQUEUE_MSG;

/**
 * creates a new priority queue. we are using
 * circular arrays to implement it.
 *
 * @param maxSize - the maximum size of the queue
 * @assert maxSize>0
 * @return returns a new, empty priority queue
 * with the capacity ability of maxSize
 */
SPBPQueue* spBPQueueCreate(int maxSize);

/**
 * copies a Queue to a new memory block
 *
 * @param source queue
 * @return the new address of the copy. null if
 * source is null
 */
SPBPQueue* spBPQueueCopy(SPBPQueue* source);

/**
 * empties and frees the memory from the queue
 *
 * @param source queue
 * @return void function
 */
void spBPQueueDestroy(SPBPQueue* source);

/*
 * removes the elements from the queue
 *
 * @param source queue
 * @return clears the queue. the queue still exists
 * but without its elements.
 */
void spBPQueueClear(SPBPQueue* source);

/**
 * @param source queue
 * @assert source!=null
 * @return the size of the queue
 */
int spBPQueueSize(SPBPQueue* source);

/**
 * @param source queue
 * @assert source != null
 * @return the maximum  size of the queue
 */
int spBPQueueGetMaxSize(SPBPQueue* source);

/**
 * inserts an element to the queue
 *
 * @param source queue, index and a value to insert as an element
 * @assert source!=NULL
 * @return a msg of success or faliure in accordance to what happened
 */
SP_BPQUEUE_MSG spBPQueueEnqueue(SPBPQueue* source, int index, double value);

/**
 * gets the element with lowest value out of the queue
 *
 * @param source queue
 * @return a msg of success or failure in accordance to what happened
 */
SP_BPQUEUE_MSG spBPQueueDequeue(SPBPQueue* source);

/**
 * put a copy of the lowest element inside res
 *
 * @param source queue and res to put inside it the copy.
 * @return a msg of success or failure in accordance to what happened
 */
SP_BPQUEUE_MSG spBPQueuePeek(SPBPQueue* source, BPQueueElement* res);

/**
 *  put a copy with of the highest element inside res
 *
 *  @param source queue and res to put inside it the copy.
 *  @return a msg of success or failure in accordance to what happened
 */
SP_BPQUEUE_MSG spBPQueuePeekLast(SPBPQueue* source, BPQueueElement* res);

/**
 * returns the minimum value in the queue
 *
 * @param source queue
 * @assert source != null && source is not empty
 * @return the minimum value
 */
double spBPQueueMinValue(SPBPQueue* source);

/**
 * returns the maximum value in the queue
 *
 * @param source queue
 * @assert source != null && source is not empty
 * @return the maximum value
 */
double spBPQueueMaxValue(SPBPQueue* source);

/**
 * returns true if source is empty
 *
 * @param source queue
 * @assert source != NULL
 * @return true if source empty
 */
bool spBPQueueIsEmpty(SPBPQueue* source);

/**
 * returns true if source is full
 *
 * @param source queue
 * @assert source != NULL
 * @return true if source full
 */
bool spBPQueueIsFull(SPBPQueue* source);

#endif

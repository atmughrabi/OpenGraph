#ifndef DYNAMICQUEUE_H
#define DYNAMICQUEUE_H

#include <stdint.h>
#include "graphConfig.h"

// A linked list (LL) node to store a queue entry
struct  QNode
{

    uint32_t key;
    struct QNode *next;
};

// The queue, front stores the front node of LL and rear stores ths
// last node of LL
struct  DynamicQueue
{
    uint32_t size;
    struct QNode *front, *rear;
};

// A utility function to create a new linked list node.
struct QNode *newQNode(uint32_t k);


// A utility function to create an empty queue
struct DynamicQueue *newDynamicQueue();


// The function to add a key k to q
void enQueue(struct DynamicQueue *q, uint32_t k);


// Function to remove a key from given queue q
struct QNode *deQueue(struct DynamicQueue *q);

struct QNode *frontDynamicQueue(struct DynamicQueue *q);

uint8_t isEmptyDynamicQueue (struct DynamicQueue *q);

#endif
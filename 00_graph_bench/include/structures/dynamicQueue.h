#ifndef DYNAMICQUEUE_H
#define DYNAMICQUEUE_H

#include <linux/types.h>
#include "graphConfig.h"

// A linked list (LL) node to store a queue entry
struct  QNode
{

    __u32 key;
    struct QNode *next;
};

// The queue, front stores the front node of LL and rear stores ths
// last node of LL
struct  DynamicQueue
{
    __u32 size;
    struct QNode *front, *rear;
};

// A utility function to create a new linked list node.
struct QNode *newQNode(__u32 k);


// A utility function to create an empty queue
struct DynamicQueue *newDynamicQueue();


// The function to add a key k to q
void enQueue(struct DynamicQueue *q, __u32 k);


// Function to remove a key from given queue q
struct QNode *deQueue(struct DynamicQueue *q);

struct QNode *frontDynamicQueue(struct DynamicQueue *q);

__u8 isEmptyDynamicQueue (struct DynamicQueue *q);

#endif
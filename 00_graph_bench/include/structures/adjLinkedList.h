#ifndef ADJLINKEDLIST_H
#define ADJLINKEDLIST_H

#include <linux/types.h>
#include "graphConfig.h"
#include "edgeList.h"


// A structure to represent an adjacency list node
struct  AdjLinkedListNode
{

    __u32 dest;
    // __u32 src;

#if WEIGHTED
    __u32 weight;
#endif

    struct AdjLinkedListNode *next;

};

// A structure to represent an adjacency list
struct  AdjLinkedList
{

    __u8 visited;
    __u32 out_degree;
    struct AdjLinkedListNode *outNodes;

#if DIRECTED
    __u32 in_degree;
    struct AdjLinkedListNode *inNodes;
#endif

};


// A utility function to create a new adjacency list node
struct AdjLinkedListNode *newAdjLinkedListOutNode(__u32 dest);
struct AdjLinkedListNode *newAdjLinkedListInNode( __u32 src);

#endif
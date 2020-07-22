#ifndef ADJLINKEDLIST_H
#define ADJLINKEDLIST_H

#include <stdint.h>
#include "graphConfig.h"
#include "edgeList.h"


// A structure to represent an adjacency list node
struct  AdjLinkedListNode
{

    uint32_t dest;
    // uint32_t src;

#if WEIGHTED
    float weight;
#endif

    struct AdjLinkedListNode *next;

};

// A structure to represent an adjacency list
struct  AdjLinkedList
{

    uint8_t visited;
    uint32_t out_degree;
    struct AdjLinkedListNode *outNodes;

#if DIRECTED
    uint32_t in_degree;
    struct AdjLinkedListNode *inNodes;
#endif

};


// A utility function to create a new adjacency list node
struct AdjLinkedListNode *newAdjLinkedListOutNode(uint32_t dest);
struct AdjLinkedListNode *newAdjLinkedListInNode( uint32_t src);

#endif
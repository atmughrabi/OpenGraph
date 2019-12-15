#ifndef GRAPHADJLINKEDLIST_H
#define GRAPHADJLINKEDLIST_H

#include <linux/types.h>
#include <omp.h>

#include "graphConfig.h"
#include "adjLinkedList.h"
#include "edgeList.h"

// A structure to represent a GraphAdjLinkedList. A GraphAdjLinkedList
// is an array of adjacency lists.
// Size of array will be V (number of vertices
// in GraphAdjLinkedList)
struct  GraphAdjLinkedList
{
    __u32 num_vertices;
    __u32 num_edges;
  
#if WEIGHTED
    __u32 max_weight;
#endif

    struct AdjLinkedList *vertices;

};


// A utility function that creates a GraphAdjLinkedList of V vertices
struct GraphAdjLinkedList *graphAdjLinkedListGraphNew(__u32 V);
struct GraphAdjLinkedList *graphAdjLinkedListEdgeListNew(struct EdgeList *edgeList);
void graphAdjLinkedListPrint(struct GraphAdjLinkedList *graphAdjLinkedList);
void graphAdjLinkedListFree(struct GraphAdjLinkedList *graphAdjLinkedList);
void adjLinkedListAddEdge(struct GraphAdjLinkedList *graphAdjLinkedList, struct EdgeList *edge, __u32 i, omp_lock_t *vertex_lock);
void   graphAdjLinkedListPrintMessageWithtime(const char *msg, double time);
struct GraphAdjLinkedList *graphAdjLinkedListPreProcessingStep (struct Arguments *arguments);

#endif



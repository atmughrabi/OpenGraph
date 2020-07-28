#ifndef GRAPHADJLINKEDLIST_H
#define GRAPHADJLINKEDLIST_H

#include <stdint.h>
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
    uint32_t num_vertices;
    uint32_t num_edges;
  	uint32_t avg_degree;
#if WEIGHTED
    float max_weight;
#endif

    struct AdjLinkedList *vertices;

};


// A utility function that creates a GraphAdjLinkedList of V vertices
struct GraphAdjLinkedList *graphAdjLinkedListGraphNew(uint32_t V);
struct GraphAdjLinkedList *graphAdjLinkedListEdgeListNew(struct EdgeList *edgeList);
void graphAdjLinkedListPrint(struct GraphAdjLinkedList *graphAdjLinkedList);
void graphAdjLinkedListFree(struct GraphAdjLinkedList *graphAdjLinkedList);
void adjLinkedListAddEdge(struct GraphAdjLinkedList *graphAdjLinkedList, struct EdgeList *edge, uint32_t i, omp_lock_t *vertex_lock);
void   graphAdjLinkedListPrintMessageWithtime(const char *msg, double time);
struct GraphAdjLinkedList *graphAdjLinkedListPreProcessingStep (struct Arguments *arguments);

#endif



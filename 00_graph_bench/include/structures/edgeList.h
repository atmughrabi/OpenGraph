#ifndef EDGELIST_H
#define EDGELIST_H

#include <stdint.h>
#include "graphConfig.h"
#include "edgeList.h"

struct  EdgeList
{

    uint32_t num_edges;
    uint32_t num_vertices;
#if WEIGHTED
    float max_weight;
    float *edges_array_weight;
#endif
    uint32_t *edges_array_src;
    uint32_t *edges_array_dest;
};


uint32_t maxTwoIntegers(uint32_t num1, uint32_t num2);
void edgeListPrint(struct EdgeList *edgeList);
void freeEdgeList( struct EdgeList *edgeList);
char *readEdgeListstxt(const char *fname, uint32_t weighted);
struct EdgeList *readEdgeListsbin(const char *fname, uint8_t inverse, uint32_t symmetric, uint32_t weighted);
struct EdgeList *readEdgeListsMem( struct EdgeList *edgeListmem,  uint8_t inverse, uint32_t symmetric, uint32_t weighted);
struct EdgeList *newEdgeList(uint32_t num_edges);
void writeEdgeListToTXTFile(struct EdgeList *edgeList, const char *fname);
struct EdgeList *removeDulpicatesSelfLoopEdges( struct EdgeList *edgeList);

#endif
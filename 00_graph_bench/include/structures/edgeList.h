#ifndef EDGELIST_H
#define EDGELIST_H

#include <linux/types.h>
#include "graphConfig.h"
#include "edgeList.h"

struct  EdgeList
{

    __u32 num_edges;
    __u32 num_vertices;
#if WEIGHTED
    __u32 max_weight;
    __u32 *edges_array_weight;
#endif
    __u32 *edges_array_src;
    __u32 *edges_array_dest;
};


__u32 maxTwoIntegers(__u32 num1, __u32 num2);
void edgeListPrint(struct EdgeList *edgeList);
void freeEdgeList( struct EdgeList *edgeList);
char *readEdgeListstxt(const char *fname, __u32 weighted);
struct EdgeList *readEdgeListsbin(const char *fname, __u8 inverse, __u32 symmetric, __u32 weighted);
struct EdgeList *readEdgeListsMem( struct EdgeList *edgeListmem,  __u8 inverse, __u32 symmetric, __u32 weighted);
struct EdgeList *newEdgeList(__u32 num_edges);
void writeEdgeListToTXTFile(struct EdgeList *edgeList, const char *fname);
struct EdgeList *removeDulpicatesSelfLoopEdges( struct EdgeList *edgeList);

#endif
#ifndef ADJARRAYLIST_H
#define ADJARRAYLIST_H

#include <stdint.h>
#include "graphConfig.h"
#include "edgeList.h"




// A structure to represent an adjacency list
struct  AdjArrayList
{

    uint32_t out_degree;
    struct EdgeList *outNodes;

#if DIRECTED
    uint32_t in_degree;
    struct EdgeList *inNodes;
#endif
};


// A utility function to create a new adjacency list node
void adjArrayListPrint(struct AdjArrayList *adjArrayList);
struct AdjArrayList *adjArrayListNew();
struct AdjArrayList *adjArrayListCreateNeighbourList(struct AdjArrayList *adjArrayList);
struct AdjArrayList *adjArrayListCreateNeighbourListOutNodes(struct AdjArrayList *adjArrayList);
struct AdjArrayList *adjArrayListCreateNeighbourListInNodes(struct AdjArrayList *adjArrayList);
void adjArrayListFree(struct AdjArrayList *adjArrayList);


#endif
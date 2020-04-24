#ifndef DFS_H
#define DFS_H

#include <stdint.h>

#include "graphConfig.h"
#include "arrayStack.h"
#include "bitmap.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"


// ********************************************************************************************
// ***************					Stats DataStructure							 **************
// ********************************************************************************************

struct DFSStats
{
    uint32_t *distances;
    int *parents;
    uint32_t  processed_nodes;
    uint32_t  num_vertices;
    double time_total;
};

struct DFSStats *newDFSStatsGraphCSR(struct GraphCSR *graph);
struct DFSStats *newDFSStatsGraphGrid(struct GraphGrid *graph);
struct DFSStats *newDFSStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct DFSStats *newDFSStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeDFSStats(struct DFSStats *stats);


// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct DFSStats  *depthFirstSearchGraphCSR(uint32_t source, struct GraphCSR *graph);
struct DFSStats  *depthFirstSearchGraphCSRBase(uint32_t source, struct GraphCSR *graph);
struct DFSStats  *pDepthFirstSearchGraphCSR(uint32_t source, struct GraphCSR *graph);
void parallelDepthFirstSearchGraphCSRTask(uint32_t source, struct GraphCSR *graph, struct DFSStats *stats);

#endif
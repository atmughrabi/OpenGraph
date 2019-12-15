#ifndef DFS_H
#define DFS_H

#include <linux/types.h>

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
    __u32 *distances;
    int *parents;
    __u32  processed_nodes;
    __u32  num_vertices;
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

struct DFSStats  *depthFirstSearchGraphCSR(__u32 source, struct GraphCSR *graph);
struct DFSStats  *depthFirstSearchGraphCSRBase(__u32 source, struct GraphCSR *graph);
struct DFSStats  *pDepthFirstSearchGraphCSR(__u32 source, struct GraphCSR *graph);
void parallelDepthFirstSearchGraphCSRTask(__u32 source, struct GraphCSR *graph, struct DFSStats *stats);

#endif
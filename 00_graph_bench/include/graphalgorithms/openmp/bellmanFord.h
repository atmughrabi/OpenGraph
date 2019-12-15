#ifndef BELLMANFORD_H
#define BELLMANFORD_H

#include <linux/types.h>

#include "graphConfig.h"
#include "edgeList.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"


// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct BellmanFordStats
{
    __u32 *distances;
    __u32 *parents;
    __u32  processed_nodes;
    __u32 num_vertices;
    double time_total;
};

struct BellmanFordStats *newBellmanFordStatsGraphCSR(struct GraphCSR *graph);
struct BellmanFordStats *newBellmanFordStatsGraphGrid(struct GraphGrid *graph);
struct BellmanFordStats *newBellmanFordStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct BellmanFordStats *newBellmanFordStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeBellmanFordStats(struct BellmanFordStats *stats);


// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************
__u32 bellmanFordAtomicMin(__u32 *dist, __u32 new);
__u32 bellmanFordCompareDistanceArrays(struct BellmanFordStats *stats1, struct BellmanFordStats *stats2);
int bellmanFordAtomicRelax(__u32 src, __u32 dest, __u32 weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext);
int bellmanFordRelax(__u32 src, __u32 dest, __u32 weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext);
void durstenfeldShuffle(__u32 *vertices, __u32 size);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphGrid(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphGrid *graph);

struct BellmanFordStats *bellmanFordPullRowGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph);
struct BellmanFordStats *bellmanFordPushColumnGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph);



// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphCSR(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphCSR *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordRandomizedDataDrivenPushGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph);
void bellmanFordSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjArrayList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjLinkedList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph);

struct BellmanFordStats *bellmanFordPullGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph);
struct BellmanFordStats *bellmanFordPushGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph);

#endif
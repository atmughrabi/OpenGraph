#ifndef TRIANGLECOUNT_H
#define TRIANGLECOUNT_H

#include <linux/types.h>
#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct TCStats
{
    __u32 num_vertices;
    __u64 *counts;
    __u64 total_counts;
    double time_total;
};

struct TCStats *newTCStatsGraphCSR(struct GraphCSR *graph);
struct TCStats *newTCStatsGraphGrid(struct GraphGrid *graph);
struct TCStats *newTCStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct TCStats *newTCStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
void freeTCStats(struct TCStats *stats);

// ********************************************************************************************
// ***************                  Helper Functions                             **************
// ********************************************************************************************

__u32 minTwoNodes(__u32 node_v, __u32 node_u, __u32 degree_v, __u32 degree_u);
__u32 maxTwoNodes(__u32 node_v, __u32 node_u, __u32 degree_v, __u32 degree_u);
__u32 countIntersectionsBinarySearch(__u32 u, __u32 v, struct GraphCSR *graph);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct TCStats *triangleCountGraphCSR(__u32 pushpull, struct GraphCSR *graph);
struct TCStats *triangleCountBasicGraphCSR(struct GraphCSR *graph);
struct TCStats *triangleCountPullGraphCSR(struct GraphCSR *graph);
struct TCStats *triangleCountPushGraphCSR(struct GraphCSR *graph);
struct TCStats *triangleCountBinaryIntersectionGraphCSR(struct GraphCSR *graph);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct TCStats *triangleCountGraphGrid(__u32 pushpull, struct GraphGrid *graph);
struct TCStats *triangleCountRowGraphGrid(struct GraphGrid *graph);
struct TCStats *triangleCountColumnGraphGrid(struct GraphGrid *graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjArrayList(__u32 pushpull, struct GraphAdjArrayList *graph);
struct TCStats *triangleCountPullGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct TCStats *triangleCountPushGraphAdjArrayList(struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjLinkedList(__u32 pushpull, struct GraphAdjLinkedList *graph);
struct TCStats *triangleCountPullGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
struct TCStats *triangleCountPushGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

#endif
#ifndef SPMV_H
#define SPMV_H

#include <linux/types.h>

#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"


// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************


struct SPMVStats
{

    __u32 iterations;
    __u32 num_vertices;
    float *vector_output;
    float *vector_input;
    double time_total;
};

struct SPMVStats *newSPMVStatsGraphCSR(struct GraphCSR *graph);
struct SPMVStats *newSPMVStatsGraphGrid(struct GraphGrid *graph);
struct SPMVStats *newSPMVStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct SPMVStats *newSPMVStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeSPMVStats(struct SPMVStats *stats);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphGrid(__u32 iterations, __u32 pushpull, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowGraphGrid(__u32 iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnGraphGrid(__u32 iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowFixedPointGraphGrid(__u32 iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnFixedPointGraphGrid(__u32 iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphCSR(__u32 iterations, __u32 pushpull, struct GraphCSR *graph);
struct SPMVStats *SPMVPullGraphCSR(__u32 iterations, struct GraphCSR *graph);
struct SPMVStats *SPMVPushGraphCSR(__u32 iterations, struct GraphCSR *graph);

struct SPMVStats *SPMVPullFixedPointGraphCSR(__u32 iterations, struct GraphCSR *graph);
struct SPMVStats *SPMVPushFixedPointGraphCSR(__u32 iterations, struct GraphCSR *graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjArrayList(__u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPullGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjLinkedList(__u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPullGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);

#endif
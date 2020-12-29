#ifndef SPMV_H
#define SPMV_H

#include <stdint.h>

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

    uint32_t iterations;
    struct Arguments *arguments;
    uint32_t num_vertices;
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

struct SPMVStats *SPMVGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowFixedPointGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnFixedPointGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct SPMVStats *SPMVPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct SPMVStats *SPMVPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

struct SPMVStats *SPMVPullFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct SPMVStats *SPMVPushFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);


// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

#endif
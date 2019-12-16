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

struct SPMVStats *SPMVGraphGrid(uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowGraphGrid(uint32_t iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnGraphGrid(uint32_t iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPullRowFixedPointGraphGrid(uint32_t iterations, struct GraphGrid *graph);
struct SPMVStats *SPMVPushColumnFixedPointGraphGrid(uint32_t iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphCSR(uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph);
struct SPMVStats *SPMVPullGraphCSR(uint32_t iterations, struct GraphCSR *graph);
struct SPMVStats *SPMVPushGraphCSR(uint32_t iterations, struct GraphCSR *graph);

struct SPMVStats *SPMVPullFixedPointGraphCSR(uint32_t iterations, struct GraphCSR *graph);
struct SPMVStats *SPMVPushFixedPointGraphCSR(uint32_t iterations, struct GraphCSR *graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjArrayList(uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPullGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjLinkedList(uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPullGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);

struct SPMVStats *SPMVPullFixedPointGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);
struct SPMVStats *SPMVPushFixedPointGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);

#endif
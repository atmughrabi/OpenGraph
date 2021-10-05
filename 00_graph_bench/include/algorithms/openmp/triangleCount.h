#ifndef TRIANGLECOUNT_H
#define TRIANGLECOUNT_H

#include <stdint.h>
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
    uint32_t num_vertices;
    uint64_t *counts;
    uint64_t total_counts;
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

uint32_t minTwoNodes(uint32_t node_v, uint32_t node_u, uint32_t degree_v, uint32_t degree_u);
uint32_t maxTwoNodes(uint32_t node_v, uint32_t node_u, uint32_t degree_v, uint32_t degree_u);
uint32_t countIntersectionsBinarySearch(uint32_t u, uint32_t v, struct GraphCSR *graph);

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct TCStats *triangleCountGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct TCStats *triangleCountBasicGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct TCStats *triangleCountPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct TCStats *triangleCountPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct TCStats *triangleCountBinaryIntersectionGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct TCStats *triangleCountGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct TCStats *triangleCountRowGraphGrid(struct GraphGrid *graph);
struct TCStats *triangleCountColumnGraphGrid(struct GraphGrid *graph);

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct TCStats *triangleCountPullGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct TCStats *triangleCountPushGraphAdjArrayList(struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct TCStats *triangleCountPullGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
struct TCStats *triangleCountPushGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

#endif
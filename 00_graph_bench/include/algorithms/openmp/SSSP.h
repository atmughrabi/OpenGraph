#ifndef SSSP_H
#define SSSP_H

#include <stdint.h>

#include "graphConfig.h"
#include "edgeList.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct SSSPStats
{
    uint32_t *distances;
    uint32_t *parents;
    uint32_t *buckets_map;
    uint32_t  bucket_counter;
    uint32_t  bucket_current;
    uint32_t  buckets_total;
    uint32_t  processed_nodes;
    uint32_t  delta;
    uint32_t num_vertices;
    double time_total;
};

struct SSSPStats *newSSSPStatsGeneral(uint32_t num_vertices, uint32_t delta);
struct SSSPStats *newSSSPStatsGraphCSR(struct GraphCSR *graph, uint32_t delta);
struct SSSPStats *newSSSPStatsGraphGrid(struct GraphGrid *graph, uint32_t delta);
struct SSSPStats *newSSSPStatsGraphAdjArrayList(struct GraphAdjArrayList *graph, uint32_t delta);
struct SSSPStats *newSSSPStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph, uint32_t delta);

void freeSSSPStats(struct SSSPStats *stats);

// ********************************************************************************************
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************
uint32_t SSSPAtomicMin(uint32_t *dist, uint32_t new);
uint32_t SSSPCompareDistanceArrays(struct SSSPStats *stats1, struct SSSPStats *stats2);
int SSSPAtomicRelax(uint32_t src, uint32_t dest, float weight, struct SSSPStats *stats);
int SSSPRelax(uint32_t src, uint32_t dest, float weight, struct SSSPStats *stats);

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);

struct SSSPStats *SSSPPullRowGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct SSSPStats *SSSPPushColumnGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);



// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

struct SSSPStats *SSSPDataDrivenPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct SSSPStats *SSSPDataDrivenPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct SSSPStats *SSSPDataDrivenSplitPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
void SSSPSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus, uint32_t delta);

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

struct SSSPStats *SSSPDataDrivenPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct SSSPStats *SSSPDataDrivenPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

struct SSSPStats *SSSPPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct SSSPStats *SSSPPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);


#endif
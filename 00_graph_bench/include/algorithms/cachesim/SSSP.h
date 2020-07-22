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
void durstenfeldShuffle(uint32_t *vertices, uint32_t size);

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphGrid(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph, uint32_t delta);

struct SSSPStats *SSSPPullRowGraphGrid(uint32_t source,  uint32_t iterations, struct GraphGrid *graph, uint32_t delta);
struct SSSPStats *SSSPPushColumnGraphGrid(uint32_t source,  uint32_t iterations, struct GraphGrid *graph, uint32_t delta);



// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphCSR(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph, uint32_t delta);

struct SSSPStats *SSSPDataDrivenPullGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta);
struct SSSPStats *SSSPDataDrivenPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta);
struct SSSPStats *SSSPDataDrivenSplitPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta);
void SSSPSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus, uint32_t delta);

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphAdjArrayList(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph, uint32_t delta);

struct SSSPStats *SSSPDataDrivenPullGraphAdjArrayList(uint32_t source,  uint32_t iterations, struct GraphAdjArrayList *graph, uint32_t delta);
struct SSSPStats *SSSPDataDrivenPushGraphAdjArrayList(uint32_t source,  uint32_t iterations, struct GraphAdjArrayList *graph, uint32_t delta);

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphAdjLinkedList(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph, uint32_t delta);

struct SSSPStats *SSSPPullGraphAdjLinkedList(uint32_t source,  uint32_t iterations, struct GraphAdjLinkedList *graph, uint32_t delta);
struct SSSPStats *SSSPPushGraphAdjLinkedList(uint32_t source,  uint32_t iterations, struct GraphAdjLinkedList *graph, uint32_t delta);


#endif
#ifndef SSSP_H
#define SSSP_H

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

struct SSSPStats
{
    __u32 *distances;
    __u32 *parents;
    __u32 *buckets_map;
    __u32  bucket_counter;
    __u32  bucket_current;
    __u32  buckets_total;
    __u32  processed_nodes;
    __u32  delta;
    __u32 num_vertices;
    double time_total;
};

struct SSSPStats *newSSSPStatsGeneral(__u32 num_vertices, __u32 delta);
struct SSSPStats *newSSSPStatsGraphCSR(struct GraphCSR *graph, __u32 delta);
struct SSSPStats *newSSSPStatsGraphGrid(struct GraphGrid *graph, __u32 delta);
struct SSSPStats *newSSSPStatsGraphAdjArrayList(struct GraphAdjArrayList *graph, __u32 delta);
struct SSSPStats *newSSSPStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph, __u32 delta);

void freeSSSPStats(struct SSSPStats *stats);

// ********************************************************************************************
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************
__u32 SSSPAtomicMin(__u32 *dist, __u32 new);
__u32 SSSPCompareDistanceArrays(struct SSSPStats *stats1, struct SSSPStats *stats2);
int SSSPAtomicRelax(__u32 src, __u32 dest, __u32 weight, struct SSSPStats *stats);
int SSSPRelax(__u32 src, __u32 dest, __u32 weight, struct SSSPStats *stats);
void durstenfeldShuffle(__u32 *vertices, __u32 size);

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct SSSPStats * SSSPGraphGrid(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphGrid *graph, __u32 delta);

struct SSSPStats *SSSPPullRowGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph, __u32 delta);
struct SSSPStats *SSSPPushColumnGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph, __u32 delta);



// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct SSSPStats *SSSPGraphCSR(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphCSR *graph, __u32 delta);

struct SSSPStats *SSSPDataDrivenPullGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph, __u32 delta);
struct SSSPStats *SSSPDataDrivenPushGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph, __u32 delta);
void SSSPSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus, __u32 delta);

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct SSSPStats * SSSPGraphAdjArrayList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph, __u32 delta);

struct SSSPStats *SSSPDataDrivenPullGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph, __u32 delta);
struct SSSPStats *SSSPDataDrivenPushGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph, __u32 delta);

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct SSSPStats * SSSPGraphAdjLinkedList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph, __u32 delta);

struct SSSPStats *SSSPPullGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph, __u32 delta);
struct SSSPStats *SSSPPushGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph, __u32 delta);


#endif
#ifndef EPOCHREORDER_H
#define EPOCHREORDER_H

#include <linux/types.h>
#include "graphCSR.h"
#include "BFS.h"
#include "bitmap.h"
#include "arrayQueue.h"

struct EpochReorder
{
    __u32 softcounter;
    __u32 hardcounter;
    __u32 softThreshold;
    __u32 hardThreshold;
    __u32 rrIndex;
    __u32 numCounters;  //frequncy[numcounters][numverticies]
    __u32 numVertices;
    __u32 reusecounter;
 	 struct Bitmap *recencyBits;
    __u32 **frequency;
    __u32 **reuse;
    __u32 *base_reuse;
   

};

struct EpochReorder *newEpochReoder( __u32 softThreshold, __u32 hardThreshold, __u32 numCounters, __u32 numVertices);
void freeEpochReorder(struct EpochReorder *epochReorder);

__u32 *epochReorderPageRank(struct GraphCSR *graph);
__u32 *epochReorderRecordBFS(struct GraphCSR *graph);
float *epochReorderPageRankPullGraphCSR(struct EpochReorder *epochReorder, double epsilon,  __u32 iterations, struct GraphCSR *graph);

void epochReorderBreadthFirstSearchGraphCSR(struct EpochReorder *epochReorder, __u32 source, struct GraphCSR *graph);
__u32 epochReorderBottomUpStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
__u32 epochReorderTopDownStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);

__u32 *epochReorderCreateLabels(struct EpochReorder *epochReorder);
void epochReorderIncrementCounters(struct EpochReorder *epochReorder, __u32 v);
void atomicEpochReorderIncrementCounters(struct EpochReorder *epochReorder, __u32 v);
__u32 epochAtomicMin(__u32 *dist, __u32 newValue);

void printEpochs(struct EpochReorder *epochReorder);
void radixSortCountSortEdgesByEpochs (__u32 **histValues, __u32 **histValuesTemp, __u32 **histMaps, __u32 **histMapsTemp, __u32 **labels, __u32 **labelsTemp, __u32 radix, __u32 buckets, __u32 *buckets_count, __u32 num_vertices);
__u32 *radixSortEdgesByEpochs (__u32 *histValues, __u32 *histMaps, __u32 *labels, __u32 num_vertices);

#endif
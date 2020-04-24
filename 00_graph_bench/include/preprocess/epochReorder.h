#ifndef EPOCHREORDER_H
#define EPOCHREORDER_H

#include <stdint.h>
#include "graphCSR.h"
#include "BFS.h"
#include "bitmap.h"
#include "arrayQueue.h"

struct EpochReorder
{
    uint32_t softcounter;
    uint32_t hardcounter;
    uint32_t softThreshold;
    uint32_t hardThreshold;
    uint32_t rrIndex;
    uint32_t numCounters;  //frequncy[numcounters][numverticies]
    uint32_t numVertices;
    uint32_t reusecounter;
 	 struct Bitmap *recencyBits;
    uint32_t **frequency;
    uint32_t **reuse;
    uint32_t *base_reuse;
   

};

struct EpochReorder *newEpochReoder( uint32_t softThreshold, uint32_t hardThreshold, uint32_t numCounters, uint32_t numVertices);
void freeEpochReorder(struct EpochReorder *epochReorder);

uint32_t *epochReorderPageRank(struct GraphCSR *graph);
uint32_t *epochReorderRecordBFS(struct GraphCSR *graph);
float *epochReorderPageRankPullGraphCSR(struct EpochReorder *epochReorder, double epsilon,  uint32_t iterations, struct GraphCSR *graph);

void epochReorderBreadthFirstSearchGraphCSR(struct EpochReorder *epochReorder, uint32_t source, struct GraphCSR *graph);
uint32_t epochReorderBottomUpStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
uint32_t epochReorderTopDownStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);

uint32_t *epochReorderCreateLabels(struct EpochReorder *epochReorder);
void epochReorderIncrementCounters(struct EpochReorder *epochReorder, uint32_t v);
void atomicEpochReorderIncrementCounters(struct EpochReorder *epochReorder, uint32_t v);
uint32_t epochAtomicMin(uint32_t *dist, uint32_t newValue);

void printEpochs(struct EpochReorder *epochReorder);
void radixSortCountSortEdgesByEpochs (uint32_t **histValues, uint32_t **histValuesTemp, uint32_t **histMaps, uint32_t **histMapsTemp, uint32_t **labels, uint32_t **labelsTemp, uint32_t radix, uint32_t buckets, uint32_t *buckets_count, uint32_t num_vertices);
uint32_t *radixSortEdgesByEpochs (uint32_t *histValues, uint32_t *histMaps, uint32_t *labels, uint32_t num_vertices);

#endif
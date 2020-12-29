#ifndef INCREMENTALAGGREGATION_H
#define INCREMENTALAGGREGATION_H

#include <stdint.h>

#include "graphConfig.h"
#include "arrayQueue.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

struct __attribute__((__packed__)) MyPair
{
    uint32_t degree;
    uint32_t child;
};


union Atom
{
    uint64_t atomicPair;
    struct MyPair pair;
};
// typedef struct pair { void *a[2]; } pair;

// inline
// void pair_swap(_Atomic(pair) *myPair) {
//   pair actual = { 0 };
//   pair future = { 0 };

//   while (!atomic_compare_exchange_weak(myPair, &actual, future)) {
//       future.a[0] = actual.a[1];
//       future.a[1] = actual.a[0];
//   }
// }

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct IncrementalAggregationStats
{
    uint32_t *vertices;
    uint32_t *degrees;
    uint32_t num_clusters;
    //dendogram
    uint32_t *atomDegree;
    uint32_t *atomChild;
    union Atom *atom;

    uint32_t *sibling;
    uint32_t *dest;
    uint32_t *weightSum;
    uint32_t *labels;
    uint32_t num_vertices;
    double totalQ;
    double time_total;

};

struct IncrementalAggregationStats *newIncrementalAggregationStatsGraphCSR(struct GraphCSR *graph);
struct IncrementalAggregationStats *newIncrementalAggregationStatsGraphGrid(struct GraphGrid *graph);
struct IncrementalAggregationStats *newIncrementalAggregationStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct IncrementalAggregationStats *newIncrementalAggregationStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeIncrementalAggregationStats(struct IncrementalAggregationStats *stats);

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct IncrementalAggregationStats *incrementalAggregationGraphCSR(struct GraphCSR *graph);
void findBestDestination(struct ArrayQueue *Neighbors, struct ArrayQueue *reachableSet, float *deltaQ, uint32_t *u, uint32_t degreeVout, uint32_t v, struct IncrementalAggregationStats *stats, struct GraphCSR *graph);
void traversDendrogramReachableSetDFS(uint32_t v, union Atom *atom, uint32_t *sibling, struct ArrayQueue *reachableSet);
void printSet(struct ArrayQueue *Set);
void returnReachableSetOfNodesFromDendrogram(uint32_t v, union Atom *atom, uint32_t *sibling, struct ArrayQueue *reachableSet);

uint32_t *returnLabelsOfNodesFromDendrogram(struct ArrayQueue *reachableSet, union Atom *atom, uint32_t *sibling, uint32_t num_vertices);
void traversDendrogramLabelsDFS(uint32_t *newLablesCounter, uint32_t *newLables, uint32_t v, union Atom *atom, uint32_t *sibling);

#endif
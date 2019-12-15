#ifndef INCREMENTALAGGREGATION_H
#define INCREMENTALAGGREGATION_H

#include <linux/types.h>

#include "graphConfig.h"
#include "arrayQueue.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

struct __attribute__((__packed__)) MyPair
{
    __u32 degree;
    __u32 child;
};


union Atom {
    __u64 atomicPair;
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
    __u32 *vertices;
    __u32 *degrees;
    __u32 num_clusters;
    //dendogram
    __u32 *atomDegree;
    __u32 *atomChild;
    union Atom *atom;

    __u32 *sibling;
    __u32 *dest;
    __u32 *weightSum;
    __u32 *labels;
    __u32 num_vertices;
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
void findBestDestination(struct ArrayQueue *Neighbors, struct ArrayQueue *reachableSet, float *deltaQ, __u32 *u, __u32 degreeVout, __u32 v, struct IncrementalAggregationStats *stats, struct GraphCSR *graph);
void traversDendrogramReachableSetDFS(__u32 v, union Atom *atom, __u32 *sibling, struct ArrayQueue *reachableSet);
void printSet(struct ArrayQueue *Set);
void returnReachableSetOfNodesFromDendrogram(__u32 v, union Atom *atom, __u32 *sibling, struct ArrayQueue *reachableSet);

__u32 *returnLabelsOfNodesFromDendrogram(struct ArrayQueue *reachableSet, union Atom *atom, __u32 *sibling, __u32 num_vertices);
void traversDendrogramLabelsDFS(__u32 *newLablesCounter, __u32 *newLables, __u32 v, union Atom *atom, __u32 *sibling);

#endif
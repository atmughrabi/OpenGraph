#ifndef PAGERANK_H
#define PAGERANK_H

#include <linux/types.h>

#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#define Damp 0.85f

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct PageRankStats
{
	float damp;
	float base_pr;
    __u32 iterations;
    __u32 num_vertices;
    __u32 *realRanks;
    float *pageRanks;
    double time_total;
    double error_total;
};

struct PageRankStats *newPageRankStatsGraphCSR(struct GraphCSR *graph);
struct PageRankStats *newPageRankStatsGraphGrid(struct GraphGrid *graph);
struct PageRankStats *newPageRankStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct PageRankStats *newPageRankStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freePageRankStats(struct PageRankStats *stats);


// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************

void addAtomicFixedPoint(__u64 *num, __u64 value);
void addAtomicFloat(float *num, float value);
void addAtomicDouble(double *num, double value);
void setAtomic(__u64 *num, __u64 value);

void pageRankPrint(float *pageRankArray, __u32 num_vertices);
void pageRankCompare(float *pageRankArrayOp1, float *pageRankArrayOp2);
void swapWorkLists (__u8 **workList1, __u8 **workList2);
void resetWorkList(__u8 *workList, __u32 size);
void setWorkList(__u8 *workList,  __u32 size);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphGrid(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphGrid *graph);
struct PageRankStats *pageRankPullRowGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPushColumnGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPullRowFixedPointGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPushColumnFixedPointGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphCSR(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphCSR *graph);
struct PageRankStats *pageRankPullGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);

struct PageRankStats *pageRankPullFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);

struct PageRankStats *pageRankDataDrivenPullGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankDataDrivenPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph);

// void pageRankDataDrivenPullFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph);
// void pageRankDataDrivenPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph);
// void pageRankDataDrivenPullPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphAdjArrayList(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPullGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);

struct PageRankStats *pageRankPullFixedPointGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPushFixedPointGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);

struct PageRankStats *pageRankDataDrivenPullGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankDataDrivenPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphAdjLinkedList(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPullGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);

struct PageRankStats *pageRankPullFixedPointGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPushFixedPointGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);

struct PageRankStats *pageRankDataDrivenPullGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankDataDrivenPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph);

#endif
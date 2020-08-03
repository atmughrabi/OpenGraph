#ifndef PAGERANK_H
#define PAGERANK_H

#include <stdint.h>

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
    uint32_t iterations;
    uint32_t num_vertices;
    uint32_t activeVertices;
    uint32_t *realRanks;
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
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************

void addAtomicFixedPoint(uint64_t *num, uint64_t value);
void addAtomicFloat(float *num, float value);
void addAtomicDouble(double *num, double value);
void setAtomic(uint64_t *num, uint64_t value);

void pageRankPrint(float *pageRankArray, uint32_t num_vertices);
void pageRankCompare(float *pageRankArrayOp1, float *pageRankArrayOp2);
void swapWorkLists (uint8_t **workList1, uint8_t **workList2);
void resetWorkList(uint8_t *workList, uint32_t size);
void setWorkList(uint8_t *workList,  uint32_t size);

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphGrid(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph);
struct PageRankStats *pageRankPullRowGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPushColumnGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPullRowFixedPointGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph);
struct PageRankStats *pageRankPushColumnFixedPointGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphCSR(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph);
struct PageRankStats *pageRankPullGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);

struct PageRankStats *pageRankPullFixedPoint64BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPullFixedPoint32BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPullFixedPoint16BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPullFixedPoint8BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);

struct PageRankStats *pageRankPullQuant32BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPullQuant16BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPullQuant8BitGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankPushQuantGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);

struct PageRankStats *pageRankDataDrivenPullGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankDataDrivenPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph);

// void pageRankDataDrivenPullFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph);
// void pageRankDataDrivenPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph);
// void pageRankDataDrivenPullPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph);

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphAdjArrayList(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPullGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);

struct PageRankStats *pageRankPullFixedPointGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankPushFixedPointGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);

struct PageRankStats *pageRankDataDrivenPullGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankDataDrivenPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct PageRankStats *pageRankGraphAdjLinkedList(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPullGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);

struct PageRankStats *pageRankPullFixedPointGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankPushFixedPointGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);

struct PageRankStats *pageRankDataDrivenPullGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankDataDrivenPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);
struct PageRankStats *pageRankDataDrivenPullPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph);

#endif
#ifndef BELLMANFORD_H
#define BELLMANFORD_H

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

struct BellmanFordStats
{
    uint32_t *distances;
    uint32_t *parents;
    uint32_t  processed_nodes;
    uint32_t num_vertices;
    double time_total;
};

struct BellmanFordStats *newBellmanFordStatsGraphCSR(struct GraphCSR *graph);
struct BellmanFordStats *newBellmanFordStatsGraphGrid(struct GraphGrid *graph);
struct BellmanFordStats *newBellmanFordStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct BellmanFordStats *newBellmanFordStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeBellmanFordStats(struct BellmanFordStats *stats);


// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************
uint32_t bellmanFordAtomicMin(uint32_t *dist, uint32_t new);
uint32_t bellmanFordCompareDistanceArrays(struct BellmanFordStats *stats1, struct BellmanFordStats *stats2);
int bellmanFordAtomicRelax(uint32_t src, uint32_t dest, float weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext);
int bellmanFordRelax(uint32_t src, uint32_t dest, float weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext);
void durstenfeldShuffle(mt19937state *mt19937var, uint32_t *vertices, uint32_t size);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);

struct BellmanFordStats *bellmanFordPullRowGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct BellmanFordStats *bellmanFordPushColumnGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);



// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordRandomizedDataDrivenPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
void bellmanFordSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

struct BellmanFordStats *bellmanFordPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct BellmanFordStats *bellmanFordPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

#endif
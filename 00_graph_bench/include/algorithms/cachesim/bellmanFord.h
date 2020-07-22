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
void durstenfeldShuffle(uint32_t *vertices, uint32_t size);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphGrid(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph);

struct BellmanFordStats *bellmanFordPullRowGraphGrid(uint32_t source,  uint32_t iterations, struct GraphGrid *graph);
struct BellmanFordStats *bellmanFordPushColumnGraphGrid(uint32_t source,  uint32_t iterations, struct GraphGrid *graph);



// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphCSR(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph);
struct BellmanFordStats *bellmanFordRandomizedDataDrivenPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph);
void bellmanFordSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjArrayList(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph);

struct BellmanFordStats *bellmanFordDataDrivenPullGraphAdjArrayList(uint32_t source,  uint32_t iterations, struct GraphAdjArrayList *graph);
struct BellmanFordStats *bellmanFordDataDrivenPushGraphAdjArrayList(uint32_t source,  uint32_t iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjLinkedList(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph);

struct BellmanFordStats *bellmanFordPullGraphAdjLinkedList(uint32_t source,  uint32_t iterations, struct GraphAdjLinkedList *graph);
struct BellmanFordStats *bellmanFordPushGraphAdjLinkedList(uint32_t source,  uint32_t iterations, struct GraphAdjLinkedList *graph);

#endif
#ifndef BFS_H
#define BFS_H

#include <stdint.h>

#include "graphConfig.h"
#include "arrayQueue.h"
#include "bitmap.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

// ********************************************************************************************
// ***************					Stats DataStructure							 **************
// ********************************************************************************************

struct BFSStats
{
    uint32_t *distances;
    int *parents;
    uint32_t  processed_nodes;
    uint32_t  num_vertices;
    uint32_t iteration;
    double time_total;
};

struct BFSStats *newBFSStatsGraphCSR(struct GraphCSR *graph);
struct BFSStats *newBFSStatsGraphGrid(struct GraphGrid *graph);
struct BFSStats *newBFSStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct BFSStats *newBFSStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeBFSStats(struct BFSStats *stats);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphCSR(uint32_t source,uint32_t pushpull, struct GraphCSR *graph);

struct BFSStats *breadthFirstSearchPullGraphCSR(uint32_t source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushGraphCSR(uint32_t source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphCSR(uint32_t source, struct GraphCSR *graph);

uint32_t topDownStepGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
uint32_t bottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);

// ********************************************************************************************
// ***************		CSR DataStructure/Bitmap Frontiers						 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchPushBitmapGraphCSR(uint32_t source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushDirectionOptimizedBitmapGraphCSR(uint32_t source, struct GraphCSR *graph);
uint32_t topDownStepUsingBitmapsGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct BFSStats *stats);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphGrid(uint32_t source, uint32_t pushpull, struct GraphGrid *graph);

struct BFSStats *breadthFirstSearchRowGraphGrid(uint32_t source, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGrid(uint32_t source, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGrid(struct GraphGrid *graph, struct Partition *partition, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue *localFrontierQueue, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitions(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue);


// ********************************************************************************************
// ***************					GRID DataStructure/Bitmap Frontiers			 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchRowGraphGridBitmap(uint32_t source, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGridBitmap(uint32_t source, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGridBitmap(struct GraphGrid *graph, struct Partition *partition, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitionsBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmap);


// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************


struct BFSStats *breadthFirstSearchGraphAdjArrayList(uint32_t source, uint32_t pushpull, struct GraphAdjArrayList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph);

uint32_t bottomUpStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
uint32_t topDownStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);


// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphAdjLinkedList(uint32_t source, uint32_t pushpull, struct GraphAdjLinkedList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph);

uint32_t bottomUpStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
uint32_t topDownStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);

#endif
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
    uint32_t *distances_DualOrder;
    int *parents_DualOrder;
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

struct BFSStats *breadthFirstSearchGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

struct BFSStats *breadthFirstSearchPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);

uint32_t topDownStepGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
uint32_t bottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);

// ********************************************************************************************
// ***************		CSR DataStructure/Bitmap Frontiers						 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchPushBitmapGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushDirectionOptimizedBitmapGraphCSR(struct Arguments *arguments, struct GraphCSR *graph);
uint32_t topDownStepUsingBitmapsGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct BFSStats *stats);

// ********************************************************************************************
// ***************					CSR DataStructure DualOrder				     **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphCSRDualOrder(struct Arguments *arguments, struct GraphCSR *graph);

struct BFSStats *breadthFirstSearchPullGraphCSRDualOrder(struct Arguments *arguments, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushGraphCSRDualOrder(struct Arguments *arguments, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphCSRDualOrder(struct Arguments *arguments, struct GraphCSR *graph);

uint32_t topDownStepGraphCSRDualOrder(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
uint32_t bottomUpStepGraphCSRDualOrder(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);

struct BFSStats *breadthFirstSearchRowGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGrid(struct Arguments *arguments, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGrid(struct GraphGrid *graph, struct Partition *partition, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue *localFrontierQueue, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitions(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue);


// ********************************************************************************************
// ***************					GRID DataStructure/Bitmap Frontiers			 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchRowGraphGridBitmap(struct Arguments *arguments, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGridBitmap(struct Arguments *arguments, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGridBitmap(struct GraphGrid *graph, struct Partition *partition, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitionsBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmap);


// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************


struct BFSStats *breadthFirstSearchGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph);

uint32_t bottomUpStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
uint32_t topDownStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);


// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph);

uint32_t bottomUpStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
uint32_t topDownStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void syncDualOrderParentArrays(int **parents, int **parents_DualOrder, uint32_t *labels, uint32_t num_vertices);
void syncDualOrderDistancesArrays(uint32_t *distances, uint32_t *distances_DualOrder, uint32_t *labels, uint32_t num_vertices);

#endif
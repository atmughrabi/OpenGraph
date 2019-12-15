#ifndef BFS_H
#define BFS_H

#include <linux/types.h>

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
    __u32 *distances;
    int *parents;
    __u32  processed_nodes;
    __u32  num_vertices;
    __u32 iteration;
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

struct BFSStats *breadthFirstSearchGraphCSR(__u32 source,__u32 pushpull, struct GraphCSR *graph);

struct BFSStats *breadthFirstSearchPullGraphCSR(__u32 source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushGraphCSR(__u32 source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphCSR(__u32 source, struct GraphCSR *graph);

__u32 topDownStepGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
__u32 bottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);

// ********************************************************************************************
// ***************		CSR DataStructure/Bitmap Frontiers						 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchPushBitmapGraphCSR(__u32 source, struct GraphCSR *graph);
struct BFSStats *breadthFirstSearchPushDirectionOptimizedBitmapGraphCSR(__u32 source, struct GraphCSR *graph);
__u32 topDownStepUsingBitmapsGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct BFSStats *stats);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphGrid(__u32 source, __u32 pushpull, struct GraphGrid *graph);

struct BFSStats *breadthFirstSearchRowGraphGrid(__u32 source, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGrid(__u32 source, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGrid(struct GraphGrid *graph, struct Partition *partition, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue *localFrontierQueue, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitions(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue);


// ********************************************************************************************
// ***************					GRID DataStructure/Bitmap Frontiers			 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchRowGraphGridBitmap(__u32 source, struct GraphGrid *graph);
struct BFSStats *breadthFirstSearchColumnGraphGridBitmap(__u32 source, struct GraphGrid *graph);

void breadthFirstSearchStreamEdgesRowGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchStreamEdgesColumnGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchPartitionGraphGridBitmap(struct GraphGrid *graph, struct Partition *partition, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats);
void breadthFirstSearchSetActivePartitionsBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmap);


// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************


struct BFSStats *breadthFirstSearchGraphAdjArrayList(__u32 source, __u32 pushpull, struct GraphAdjArrayList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjArrayList(__u32 source, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjArrayList(__u32 source, struct GraphAdjArrayList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjArrayList(__u32 source, struct GraphAdjArrayList *graph);

__u32 bottomUpStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
__u32 topDownStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);


// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphAdjLinkedList(__u32 source, __u32 pushpull, struct GraphAdjLinkedList *graph);

struct BFSStats *breadthFirstSearchPullGraphAdjLinkedList(__u32 source, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchPushGraphAdjLinkedList(__u32 source, struct GraphAdjLinkedList *graph);
struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(__u32 source, struct GraphAdjLinkedList *graph);

__u32 bottomUpStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats);
__u32 topDownStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);

#endif
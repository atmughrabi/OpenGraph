#ifndef CONNECTEDCOMPONENTS_H
#define CONNECTEDCOMPONENTS_H

#include <linux/types.h>
#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

#ifndef HASHSIZE
#define HASHSIZE (1 << 16) // hash table size 256
#endif 
#define JUDYERROR_SAMPLE 1 

struct CCStats
{	
    __u32 neighbor_rounds;
    __u32 iterations;
    __u32 num_vertices;
    __u32 *components;
    __u32 *counts;
    __u32 *labels;
    double time_total;
};

struct CCStats *newCCStatsGraphCSR(struct GraphCSR *graph);
struct CCStats *newCCStatsGraphGrid(struct GraphGrid *graph);
struct CCStats *newCCStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct CCStats *newCCStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
void printCCStats(struct CCStats *stats);
void freeCCStats(struct CCStats *stats);


void printComponents(struct CCStats *stats);
// ********************************************************************************************
// ***************						 Helper Functions						 **************
// ********************************************************************************************
__u32 atomicMin(__u32 *oldValue, __u32 newValue);
void addSample(__u32 id);
void linkNodes(__u32 u, __u32 v, __u32 *components);
void compressNodes(__u32 num_vertices, __u32 *components);
__u32 sampleFrequentNode(__u32 num_vertices, __u32 num_samples, __u32 *components);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphCSR(__u32 iterations, __u32 pushpull, struct GraphCSR *graph);
struct CCStats *connectedComponentsAfforestGraphCSR(__u32 iterations, struct GraphCSR *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphCSR(__u32 iterations, struct GraphCSR *graph);
struct CCStats *connectedComponentsWeaklyGraphCSR( __u32 iterations, struct GraphCSR *graph);
__u32 connectedComponentsVerifyGraphCSR(struct CCStats *stats, struct GraphCSR *graph);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphGrid(__u32 iterations, __u32 pushpull, struct GraphGrid *graph);
struct CCStats *connectedComponentsAfforestGraphGrid(__u32 iterations, struct GraphGrid *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphGrid(__u32 iterations, struct GraphGrid *graph);
struct CCStats *connectedComponentsWeaklyGraphGrid(__u32 iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjArrayList(__u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsAfforestGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphAdjArrayList(__u32 iterations, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsWeaklyGraphAdjArrayList( __u32 iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjLinkedList(__u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsAfforestGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphAdjLinkedList(__u32 iterations, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsWeaklyGraphAdjLinkedList( __u32 iterations, struct GraphAdjLinkedList *graph);



#endif
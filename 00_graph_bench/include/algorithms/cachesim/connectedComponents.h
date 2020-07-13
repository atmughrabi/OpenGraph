#ifndef CONNECTEDCOMPONENTS_H
#define CONNECTEDCOMPONENTS_H

#include <stdint.h>
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
    uint32_t neighbor_rounds;
    uint32_t iterations;
    uint32_t num_vertices;
    uint32_t *components;
    uint32_t *counts;
    uint32_t *labels;
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
uint32_t atomicMin(uint32_t *oldValue, uint32_t newValue);
void addSample(uint32_t id);
void linkNodes(uint32_t u, uint32_t v, uint32_t *components);
void compressNodes(uint32_t num_vertices, uint32_t *components);
uint32_t sampleFrequentNode(uint32_t num_vertices, uint32_t num_samples, uint32_t *components);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphCSR(uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph);
struct CCStats *connectedComponentsAfforestGraphCSR(uint32_t iterations, struct GraphCSR *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphCSR(uint32_t iterations, struct GraphCSR *graph);
struct CCStats *connectedComponentsWeaklyGraphCSR( uint32_t iterations, struct GraphCSR *graph);
uint32_t connectedComponentsVerifyGraphCSR(struct CCStats *stats, struct GraphCSR *graph);

// ********************************************************************************************
// ***************					GRID DataStructure							 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphGrid(uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph);
struct CCStats *connectedComponentsAfforestGraphGrid(uint32_t iterations, struct GraphGrid *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphGrid(uint32_t iterations, struct GraphGrid *graph);
struct CCStats *connectedComponentsWeaklyGraphGrid(uint32_t iterations, struct GraphGrid *graph);

// ********************************************************************************************
// ***************					ArrayList DataStructure					     **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjArrayList(uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsAfforestGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph);
struct CCStats *connectedComponentsWeaklyGraphAdjArrayList( uint32_t iterations, struct GraphAdjArrayList *graph);

// ********************************************************************************************
// ***************					LinkedList DataStructure					 **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjLinkedList(uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsAfforestGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsShiloachVishkinGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph);
struct CCStats *connectedComponentsWeaklyGraphAdjLinkedList( uint32_t iterations, struct GraphAdjLinkedList *graph);



#endif
#ifndef GRAPHRUN_H
#define GRAPHRUN_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "graphConfig.h"
#include "timer.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjLinkedList.h"
#include "graphAdjArrayList.h"

#include "BFS.h"
#include "DFS.h"
#include "pageRank.h"
#include "incrementalAggregation.h"
#include "bellmanFord.h"
#include "SSSP.h"
#include "SPMV.h"
#include "connectedComponents.h"


// Random root helper functions
__u32 generateRandomRootGraphCSR(struct GraphCSR *graph);
__u32 generateRandomRootGraphGrid(struct GraphGrid *graph);
__u32 generateRandomRootGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
__u32 generateRandomRootGraphAdjArrayList(struct GraphAdjArrayList *graph);
__u32 generateRandomRootGeneral(void *graph, struct Arguments *arguments);

void freeGraphDataStructure(void *graph, __u32 datastructure);
void freeGraphStatsGeneral(void *stats, __u32 algorithm);

void writeSerializedGraphDataStructure(struct Arguments *arguments);
void readSerializeGraphDataStructure(struct Arguments *arguments);

void generateGraphPrintMessageWithtime(const char *msg, double time);
void *generateGraphDataStructure(struct Arguments *arguments);

void runGraphAlgorithms(void *graph, struct Arguments *arguments);

struct BFSStats *runBreadthFirstSearchAlgorithm(void *graph, __u32 datastructure, int root, __u32 pushpull);
struct PageRankStats *runPageRankAlgorithm(void *graph, __u32 datastructure, double epsilon, __u32 iterations, __u32 pushpull);
struct DFSStats *runDepthFirstSearchAlgorithm(void *graph, __u32 datastructure, int root);
struct IncrementalAggregationStats *runIncrementalAggregationAlgorithm(void *graph, __u32 datastructure);
struct BellmanFordStats *runBellmanFordAlgorithm(void *graph, __u32 datastructure, __u32 root, __u32 iterations, __u32 pushpull);
struct SSSPStats *runSSSPAlgorithm(void *graph, __u32 datastructure, __u32 root, __u32 iterations, __u32 pushpull, __u32 delta);
struct SPMVStats *runSPMVAlgorithm(void *graph, __u32 datastructure, __u32 iterations, __u32 pushpull);
struct CCStats *runConnectedComponentsAlgorithm(void *graph, __u32 datastructure, __u32 iterations, __u32 pushpull);
struct TCStats *runTriangleCountAlgorithm(void *graph, __u32 datastructure, __u32 pushpull);


#endif



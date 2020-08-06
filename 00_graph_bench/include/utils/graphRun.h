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
uint32_t generateRandomRootGraphCSR(struct GraphCSR *graph);
uint32_t generateRandomRootGraphGrid(struct GraphGrid *graph);
uint32_t generateRandomRootGraphAdjLinkedList(struct GraphAdjLinkedList *graph);
uint32_t generateRandomRootGraphAdjArrayList(struct GraphAdjArrayList *graph);
uint32_t generateRandomRootGeneral(void *graph, struct Arguments *arguments);

void freeGraphDataStructure(void *graph, uint32_t datastructure);
void freeGraphStatsGeneral(void *stats, uint32_t algorithm);

void writeSerializedGraphDataStructure(struct Arguments *arguments);
void readSerializeGraphDataStructure(struct Arguments *arguments);

void generateGraphPrintMessageWithtime(const char *msg, double time);
void *generateGraphDataStructure(struct Arguments *arguments);

void runGraphAlgorithms(void *graph, struct Arguments *arguments);

struct BFSStats *runBreadthFirstSearchAlgorithm(void *graph, uint32_t datastructure, int root, uint32_t pushpull);
struct PageRankStats *runPageRankAlgorithm(void *graph, uint32_t datastructure, double epsilon, uint32_t iterations, uint32_t pushpull);
struct DFSStats *runDepthFirstSearchAlgorithm(void *graph, uint32_t datastructure, int root);
struct IncrementalAggregationStats *runIncrementalAggregationAlgorithm(void *graph, uint32_t datastructure);
struct BellmanFordStats *runBellmanFordAlgorithm(void *graph, uint32_t datastructure, uint32_t root, uint32_t iterations, uint32_t pushpull);
struct SSSPStats *runSSSPAlgorithm(void *graph, uint32_t datastructure, uint32_t root, uint32_t iterations, uint32_t pushpull, uint32_t delta);
struct SPMVStats *runSPMVAlgorithm(void *graph, uint32_t datastructure, uint32_t iterations, uint32_t pushpull);
struct CCStats *runConnectedComponentsAlgorithm(void *graph, uint32_t datastructure, uint32_t iterations, uint32_t pushpull);
struct TCStats *runTriangleCountAlgorithm(void *graph, uint32_t datastructure, uint32_t pushpull);
struct BetweennessCentralityStats *runBetweennessCentralityAlgorithm(void *graph, uint32_t datastructure, uint32_t pushpull);


#endif



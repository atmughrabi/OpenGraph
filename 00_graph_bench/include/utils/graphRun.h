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
uint32_t generateRandomRootGraphCSR(mt19937state *mt19937var, struct GraphCSR *graph);
uint32_t generateRandomRootGraphGrid(mt19937state *mt19937var, struct GraphGrid *graph);
uint32_t generateRandomRootGraphAdjLinkedList(mt19937state *mt19937var, struct GraphAdjLinkedList *graph);
uint32_t generateRandomRootGraphAdjArrayList(mt19937state *mt19937var, struct GraphAdjArrayList *graph);
uint32_t generateRandomRootGeneral(struct Arguments *arguments, void *graph);

void freeGraphDataStructure(void *graph, uint32_t datastructure);
void freeGraphStatsGeneral(void *stats, uint32_t algorithm);

void writeSerializedGraphDataStructure(struct Arguments *arguments);
void readSerializeGraphDataStructure(struct Arguments *arguments);

void generateGraphPrintMessageWithtime(const char *msg, double time);
void *generateGraphDataStructure(struct Arguments *arguments);

void runGraphAlgorithms(struct Arguments *arguments, void *graph);

struct BFSStats *runBreadthFirstSearchAlgorithm(struct Arguments *arguments, void *graph);
struct PageRankStats *runPageRankAlgorithm(struct Arguments *arguments, void *graph);
struct DFSStats *runDepthFirstSearchAlgorithm(struct Arguments *arguments, void *graph);
struct IncrementalAggregationStats *runIncrementalAggregationAlgorithm(struct Arguments *arguments, void *graph);
struct BellmanFordStats *runBellmanFordAlgorithm(struct Arguments *arguments, void *graph);
struct SSSPStats *runSSSPAlgorithm(struct Arguments *arguments, void *graph);
struct SPMVStats *runSPMVAlgorithm(struct Arguments *arguments, void *graph);
struct CCStats *runConnectedComponentsAlgorithm(struct Arguments *arguments, void *graph);
struct TCStats *runTriangleCountAlgorithm(struct Arguments *arguments, void *graph);
struct BetweennessCentralityStats *runBetweennessCentralityAlgorithm(struct Arguments *arguments, void *graph);


#endif



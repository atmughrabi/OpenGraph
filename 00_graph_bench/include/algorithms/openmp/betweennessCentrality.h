#ifndef BETWEENNESSCENTRALITY_H
#define BETWEENNESSCENTRALITY_H

#include <stdint.h>

#include "graphConfig.h"
#include "edgeList.h"
#include "bitmap.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"


// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************
struct Predecessor
{
    uint32_t *nodes;
    uint32_t  degree;
};

struct BetweennessCentralityStats
{
    struct Predecessor *stack;
    uint32_t *distances;
    uint32_t *realRanks;
    int *parents;
    int *sigma;
    float    *dependency;
    float    *betweennessCentrality;
    struct Predecessor *predecessors;
    uint32_t  iteration;
    uint32_t  processed_nodes;
    uint32_t  num_vertices;
    double time_total;
};

struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphCSR(struct GraphCSR *graph);

//TODO:
struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphGrid(struct GraphGrid *graph);
struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphAdjArrayList(struct GraphAdjArrayList *graph);
struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph);

void freeBetweennessCentralityStats(struct BetweennessCentralityStats *stats);
void clearBetweennessCentralityStats(struct BetweennessCentralityStats *stats);

// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************
uint32_t generateRandomRootBetweennessCentrality(struct GraphCSR *graph);
void copyBitmapToStack(struct Bitmap *q_bitmap, struct Predecessor *stack, uint32_t num_vertices);
struct Predecessor *creatNewPredecessorList(uint32_t *degrees, uint32_t num_vertices);
struct BetweennessCentralityStats *betweennessCentralityBFSPullGraphCSR(uint32_t source, struct GraphCSR *graph, struct BetweennessCentralityStats *stats);
uint32_t betweennessCentralityBottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BetweennessCentralityStats *stats);

// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BetweennessCentralityStats *betweennessCentralityGraphCSR(uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph);
struct BetweennessCentralityStats *betweennessCentralityBrandesGraphCSR(uint32_t iterations, struct GraphCSR *graph);


#endif
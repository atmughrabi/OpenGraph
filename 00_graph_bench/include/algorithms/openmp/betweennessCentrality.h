#ifndef BETWEENNESSCENTRALITY_H
#define BETWEENNESSCENTRALITY_H

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

struct BetweennessCentralityStats
{
    uint32_t *distances;
    int *parents;
    float    *scores;
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


// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************



// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BetweennessCentralityStats *betweennessCentralityGraphCSR(uint32_t pushpull, struct GraphCSR *graph);
struct BetweennessCentralityStats *betweennessCentralityBrandesGraphCSR(struct GraphCSR *graph);


#endif
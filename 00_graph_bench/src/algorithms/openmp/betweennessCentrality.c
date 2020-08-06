// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : betweennessCentrality.c
// Create : 2019-09-28 15:20:58
// Revise : 2019-09-28 15:34:05
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#include "timer.h"
#include "myMalloc.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"
#include "graphConfig.h"

#include "graphCSR.h"
// #include "graphGrid.h"
// #include "graphAdjArrayList.h"
// #include "graphAdjLinkedList.h"

#include "betweennessCentrality.h"


// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t vertex_id;

    struct BetweennessCentralityStats *stats = (struct BetweennessCentralityStats *) my_malloc(sizeof(struct BetweennessCentralityStats));

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->scores  = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = 0;
        stats->scores[vertex_id] = 0.0f;
        if(graph->vertices->out_degree[vertex_id])
            stats->parents[vertex_id] = graph->vertices->out_degree[vertex_id] * (-1);
        else
            stats->parents[vertex_id] = -1;
    }

    return stats;

}


void freeBetweennessCentralityStats(struct BetweennessCentralityStats *stats)
{
    if(stats)
    {
        if(stats->distances)
            free(stats->distances);
        if(stats->parents)
            free(stats->parents);
        if(stats->scores)
            free(stats->scores);
        free(stats);
    }
}

// ********************************************************************************************
// ***************					Auxiliary functions  	  					 **************
// ********************************************************************************************


// ********************************************************************************************
// ***************					CSR DataStructure							 **************
// ********************************************************************************************

struct BetweennessCentralityStats *betweennessCentralityGraphCSR(uint32_t pushpull, struct GraphCSR *graph)
{
    struct BetweennessCentralityStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // Brandes
        stats = betweennessCentralityBrandesGraphCSR(graph);
        break;
    default:// Brandes
        stats = betweennessCentralityBrandesGraphCSR(graph);
        break;
    }


    return stats;
}

struct BetweennessCentralityStats *betweennessCentralityBrandesGraphCSR(struct GraphCSR *graph)
{

    struct BetweennessCentralityStats *stats = newBetweennessCentralityStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Brandes Betweenness Centrality");
    printf(" -----------------------------------------------------\n");

    return stats;

}
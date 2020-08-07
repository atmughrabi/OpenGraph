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
#include <stdint.h>
#include <omp.h>

#include "mt19937.h"
#include "timer.h"
#include "myMalloc.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"
#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "betweennessCentrality.h"


// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct BetweennessCentralityStats *newBetweennessCentralityStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t vertex_id;

    struct BetweennessCentralityStats *stats = (struct BetweennessCentralityStats *) my_malloc(sizeof(struct BetweennessCentralityStats));
    stats->betweennessCentrality = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->stack  = (struct Predecessor *) my_malloc(graph->num_vertices * sizeof(struct Predecessor));
    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->realRanks  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->stack->nodes  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->stack->degree  = 0;
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->sigma = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->dependency  = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->predecessors = creatNewPredecessorList(graph->vertices->in_degree, graph->num_vertices);

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = UINT32_MAX;
        stats->dependency[vertex_id] = 0.0f;
        stats->betweennessCentrality[vertex_id] = 0.0f;
        stats->sigma[vertex_id] = 0;
        stats->realRanks[vertex_id] = vertex_id;
        stats->stack->nodes[vertex_id] = 0;
        if(graph->vertices->out_degree[vertex_id])
            stats->parents[vertex_id] = graph->vertices->out_degree[vertex_id] * (-1);
        else
            stats->parents[vertex_id] = -1;
    }

    return stats;

}


void printRanksBetweennessCentralityStats(struct BetweennessCentralityStats *stats)
{

    uint32_t vertex_id;
    uint32_t v;
    uint32_t topK = 10;
    float bc;

    if(topK > stats->num_vertices)
        topK = stats->num_vertices;

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-33s | \n", "Vertex(id)",  "Betweenness Centrality");
    for(vertex_id = stats->num_vertices - 1; vertex_id > stats->num_vertices - topK ; vertex_id--)
    {
        v = stats->realRanks[vertex_id];
        bc = stats->betweennessCentrality[v];
        printf("| %-15u | %-33f | \n", v,  bc);
    }
    printf(" -----------------------------------------------------\n");
}

void clearBetweennessCentralityStats(struct BetweennessCentralityStats *stats)
{

    uint32_t vertex_id;

    stats->stack->degree  = 0;
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = stats->num_vertices;
    // stats->time_total = 0.0f;

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats)
    for(vertex_id = 0; vertex_id < stats->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = UINT32_MAX;
        stats->dependency[vertex_id] = 0.0f;
        stats->sigma[vertex_id] = 0;
        stats->stack->nodes[vertex_id] = 0;
        stats->predecessors[vertex_id].degree = 0;
        stats->parents[vertex_id] = -1;
    }
}

void freeBetweennessCentralityStats(struct BetweennessCentralityStats *stats)
{
    uint32_t i;

    if(stats)
    {
        if(stats->distances)
            free(stats->distances);
        if(stats->parents)
            free(stats->parents);
        if(stats->dependency)
            free(stats->dependency);
        if(stats->sigma)
            free(stats->sigma);
        if(stats->betweennessCentrality)
            free(stats->betweennessCentrality);
        if(stats->stack)
        {
            if(stats->stack->nodes)
                free(stats->stack->nodes);
        }
        free(stats->stack);
        if(stats->predecessors)
        {
            for(i = 0; i < stats->num_vertices; i++)
            {
                free( stats->predecessors[i].nodes);
            }
            free(stats->predecessors);
        }
        free(stats);
    }
}

// ********************************************************************************************
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************

uint32_t generateRandomRootBetweennessCentrality(struct GraphCSR *graph)
{

    uint32_t root = 0;

    while(1)
    {
        root = generateRandInt(mt19937var);
        if(root < graph->num_vertices)
        {
            if(graph->vertices->out_degree[root] > 0)
                break;
        }
    }

    return root;

}

struct Predecessor *creatNewPredecessorList(uint32_t *degrees, uint32_t num_vertices)
{
    struct Predecessor *predecessors = (struct Predecessor *) my_malloc( num_vertices * sizeof(struct Predecessor));

    uint32_t i;
    for(i = 0; i < num_vertices; i++)
    {
        predecessors[i].nodes = (uint32_t *) my_malloc(degrees[i] * sizeof(uint32_t));
        predecessors[i].degree = 0;
    }

    return predecessors;
}

void copyBitmapToStack(struct Bitmap *q_bitmap, struct Predecessor *stack, uint32_t num_vertices)
{
    uint32_t v;

    for (v = 0; v < num_vertices; ++v)
    {
        if(getBit(q_bitmap, v))
        {
            stack->nodes[stack->degree] = v;
            stack->degree++;
        }
    }
}

// breadth-first-search(graph, source)
//  sharedFrontierQueue ← {source}
//  next ← {}
//  parents ← [-1,-1,. . . -1]
//      while sharedFrontierQueue 6= {} do
//          top-down-step(graph, sharedFrontierQueue, next, parents)
//          sharedFrontierQueue ← next
//          next ← {}
//      end while
//  return parents

struct BetweennessCentralityStats *betweennessCentralityBFSPullGraphCSR(uint32_t source, struct GraphCSR *graph, struct BetweennessCentralityStats *stats)
{



    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-33u | \n", "SOURCE NODE", source);

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t nf = 0; // number of vertices in sharedFrontierQueue

    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;
    stats->sigma[source] = 1;
    copyBitmapToStack(sharedFrontierQueue->q_bitmap_next, stats->stack, stats->num_vertices);

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);

    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {
        nf = betweennessCentralityBottomUpStepGraphCSR(graph, sharedFrontierQueue->q_bitmap, sharedFrontierQueue->q_bitmap_next, stats);
        sharedFrontierQueue->q_bitmap_next->numSetBits = nf;
        copyBitmapToStack(sharedFrontierQueue->q_bitmap_next, stats->stack, stats->num_vertices);
        swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
        clearBitmap(sharedFrontierQueue->q_bitmap_next);
        stats->processed_nodes += nf;


    } // end while

    freeArrayQueue(sharedFrontierQueue);

    return stats;
}


// bottom-up-step(graph, sharedFrontierQueue, next, parents) //pull
//  for v ∈ vertices do
//      if parents[v] = -1 then
//          for u ∈ neighbors[v] do
//              if u ∈ sharedFrontierQueue then
//              parents[v] ← u
//              next ← next ∪ {v}
//              break
//              end if
//          end for
//      end if
//  end for

uint32_t betweennessCentralityBottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BetweennessCentralityStats *stats)
{
    uint32_t v;
    uint32_t u;
    uint32_t j;
    uint32_t edge_idx;
    uint32_t out_degree;
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;

    uint32_t nf = 0; // number of vertices in sharedFrontierQueue


#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif

    #pragma omp parallel for default(none) private(j,u,v,out_degree,edge_idx) shared(stats,bitmapCurr,bitmapNext,graph,vertices,sorted_edges_array) reduction(+:nf) schedule(dynamic, 1024)
    for(v = 0 ; v < graph->num_vertices ; v++)
    {
        out_degree = vertices->out_degree[v];
        if(stats->distances[v] == UINT32_MAX)  // optmization
        {
            edge_idx = vertices->edges_idx[v];

            for(j = edge_idx ; j < (edge_idx + out_degree) ; j++)
            {
                u = sorted_edges_array[j];
                if(getBit(bitmapCurr, u))
                {
                    // stats->parents[v] = u;
                    stats->distances[v] = stats->distances[u] + 1;

                    if(stats->distances[v] == stats->distances[u] + 1)
                    {
                        stats->sigma[v] += stats->sigma[u];
                        stats->predecessors[v].nodes[stats->predecessors[v].degree] = u;
                        stats->predecessors[v].degree++;
                    }

                    setBitAtomic(bitmapNext, v);
                    nf++;
                    // break;
                }
            }

        }

    }
    return nf;
}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct BetweennessCentralityStats *betweennessCentralityGraphCSR(uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph)
{
    struct BetweennessCentralityStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // Brandes
        stats = betweennessCentralityBrandesGraphCSR(iterations, graph);
        break;
    default:// Brandes
        stats = betweennessCentralityBrandesGraphCSR(iterations, graph);
        break;
    }


    return stats;
}

struct BetweennessCentralityStats *betweennessCentralityBrandesGraphCSR(uint32_t iterations, struct GraphCSR *graph)
{

    struct BetweennessCentralityStats *stats = newBetweennessCentralityStatsGraphCSR(graph);

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Brandes Betweenness Centrality");
    printf(" -----------------------------------------------------\n");

    uint32_t iter;
    uint32_t s;
    uint32_t v;
    uint32_t w;
    uint32_t u;
    uint32_t t;

    Start(timer);
    for(iter = 0 ; iter < iterations ; iter++)
    {
        s = generateRandomRootBetweennessCentrality(graph);
        Start(timer_inner);
        clearBetweennessCentralityStats(stats);

        stats = betweennessCentralityBFSPullGraphCSR(s, graph, stats);


        // #pragma omp parallel for
        for (t = stats->stack->degree - 1; t > 0; t--)
        {
            w = stats->stack->nodes[t];

            for(u = 0 ; u < stats->predecessors[w].degree; u++)
            {
                v = stats->predecessors[w].nodes[u];
                if(stats->sigma[w] != 0)
                {
                    stats->dependency[v] +=  (stats->sigma[v] * 1.0 / stats->sigma[w]) * (1 + stats->dependency[w]);
                }
            }

            if (w != s)
            {
                stats->betweennessCentrality[w] += stats->dependency[w] / 2;
            }
        }
        Stop(timer_inner);
        stats->time_total += Seconds(timer_inner);

        printf("| %-15s | %-15f | %-15u | \n", "Iter.Time", Seconds(timer_inner), stats->processed_nodes);
    }
    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-33f | \n", "Avg.Time", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    return stats;

}
// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : SSSP.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:34:11
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <limits.h> //UINT_MAX

#include "timer.h"
#include "myMalloc.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"
#include "sortRun.h"

#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "SSSP.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************
struct SSSPStats *newSSSPStatsGeneral(uint32_t num_vertices, uint32_t delta)
{

    uint32_t v;

    struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

    stats->bucket_counter = 0;
    stats->delta = delta;
    stats->bucket_current = 0;
    stats->processed_nodes = 0;
    stats->buckets_total = 0;
    stats->time_total = 0.0;
    stats->num_vertices = num_vertices;

    stats->distances  = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    stats->parents = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    stats->buckets_map = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < num_vertices; v++)
    {
        stats->buckets_map[v] = UINT_MAX / 2;
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;

}


struct SSSPStats *newSSSPStatsGraphCSR(struct GraphCSR *graph, uint32_t delta)
{

    uint32_t v;

    struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

    stats->bucket_counter = 0;
    stats->delta = delta;
    stats->bucket_current = 0;
    stats->processed_nodes = 0;
    stats->buckets_total = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->buckets_map = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->buckets_map[v] = UINT_MAX / 2;
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;

}

struct SSSPStats *newSSSPStatsGraphGrid(struct GraphGrid *graph, uint32_t delta)
{

    uint32_t v;

    struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

    stats->bucket_counter = 0;
    stats->delta = delta;
    stats->bucket_current = 0;
    stats->processed_nodes = 0;
    stats->buckets_total = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->buckets_map = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->buckets_map[v] = UINT_MAX / 2;
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;
}

struct SSSPStats *newSSSPStatsGraphAdjArrayList(struct GraphAdjArrayList *graph, uint32_t delta)
{

    uint32_t v;

    struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

    stats->bucket_counter = 0;
    stats->delta = delta;
    stats->bucket_current = 0;
    stats->processed_nodes = 0;
    stats->buckets_total = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->buckets_map = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->buckets_map[v] = UINT_MAX / 2;
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;
}

struct SSSPStats *newSSSPStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph, uint32_t delta)
{

    uint32_t v;

    struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

    stats->bucket_counter = 0;
    stats->delta = delta;
    stats->bucket_current = 0;
    stats->processed_nodes = 0;
    stats->buckets_total = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->buckets_map = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->buckets_map[v] = UINT_MAX / 2;
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;
}

void freeSSSPStats(struct SSSPStats *stats)
{

    if(stats)
    {
        if(stats->distances)
            free(stats->distances);
        if(stats->parents)
            free(stats->parents);
        if(stats->buckets_map)
            free(stats->buckets_map);
        free(stats);
    }

}



// ********************************************************************************************
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************

uint32_t SSSPAtomicMin(uint32_t *dist, uint32_t newValue)
{

    uint32_t oldValue;
    uint32_t flag = 0;

    do
    {

        oldValue = *dist;
        if(oldValue > newValue)
        {
            if(__sync_bool_compare_and_swap(dist, oldValue, newValue))
            {
                flag = 1;
            }
        }
        else
        {
            return 0;
        }
    }
    while(!flag);

    return 1;
}

uint32_t SSSPCompareDistanceArrays(struct SSSPStats *stats1, struct SSSPStats *stats2)
{

    uint32_t v = 0;


    for(v = 0 ; v < stats1->num_vertices ; v++)
    {

        if(stats1->distances[v] != stats2->distances[v] && stats2->distances[v] != ( UINT_MAX / 2) )
        {

            // printf("v %u s1 %u s2 %u \n",v,stats1->distances[v],stats2->distances[v]);
            // return 0;
        }
        // else if(stats1->distances[v] != UINT_MAX/2)


    }

    return 1;

}

int SSSPAtomicRelax(uint32_t src, uint32_t dest, float weight, struct SSSPStats *stats)
{
    uint32_t oldParent, newParent;
    uint32_t oldDistanceV = UINT_MAX / 2;
    uint32_t oldDistanceU = UINT_MAX / 2;
    uint32_t newDistance = UINT_MAX / 2;
    uint32_t oldBucket = UINT_MAX / 2;
    uint32_t newBucket = UINT_MAX / 2;
    uint32_t flagu = 0;
    uint32_t flagv = 0;


    do
    {

        flagu = 0;
        flagv = 0;



        oldDistanceV = stats->distances[src];
        oldDistanceU = stats->distances[dest];
        oldParent = stats->parents[dest];
        oldBucket = stats->buckets_map[dest];
        newDistance = oldDistanceV + weight;

        if( oldDistanceU > newDistance )
        {

            newParent = src;
            newDistance = oldDistanceV + weight;
            newBucket = newDistance / stats->delta;

            if(__sync_bool_compare_and_swap(&(stats->distances[src]), oldDistanceV, oldDistanceV))
            {
                flagv = 1;
            }

            if(__sync_bool_compare_and_swap(&(stats->distances[dest]), oldDistanceU, newDistance) && flagv)
            {

                flagu = 1;
            }


        }
        else
        {


            return 0;
        }

    }
    while (!flagu);

    __sync_bool_compare_and_swap(&(stats->parents[dest]), oldParent, newParent);

    if(__sync_bool_compare_and_swap(&(stats->buckets_map[dest]), oldBucket, newBucket))
    {
        if(oldBucket == UINT_MAX / 2)
            #pragma omp atomic update
            stats->buckets_total++;
    }

    if(__sync_bool_compare_and_swap(&(stats->bucket_current), newBucket, stats->bucket_current))
    {
        stats->bucket_counter = 1;
    }



    return 1;

}



int SSSPRelax(uint32_t src, uint32_t dest, float weight, struct SSSPStats *stats)
{


    uint32_t newDistance = stats->distances[src] + weight;
    uint32_t bucket = UINT_MAX / 2;


    if( stats->distances[dest] > newDistance )
    {

        if(stats->buckets_map[dest] == UINT_MAX / 2)
            stats->buckets_total++;

        stats->distances[dest] = newDistance;
        stats->parents[dest] = src;

        bucket = newDistance / stats->delta;

        // if(stats->buckets_map[dest] > bucket)
        stats->buckets_map[dest] = bucket;

        if (bucket ==  stats->bucket_current)
            stats->bucket_counter = 1;


        return 1;

    }


    return 0;
}

void SSSPPrintStats(struct SSSPStats *stats)
{
    uint32_t v;

    for(v = 0; v < stats->num_vertices; v++)
    {

        if(stats->distances[v] != UINT_MAX / 2)
        {

            printf("d %u \n", stats->distances[v]);

        }


    }


}

void SSSPPrintStatsDetails(struct SSSPStats *stats)
{
    uint32_t v;
    uint32_t minDistance = UINT_MAX / 2;
    uint32_t maxDistance = 0;
    uint32_t numberOfDiscoverNodes = 0;


    #pragma omp parallel for reduction(max:maxDistance) reduction(+:numberOfDiscoverNodes) reduction(min:minDistance)
    for(v = 0; v < stats->num_vertices; v++)
    {

        if(stats->distances[v] != UINT_MAX / 2)
        {

            numberOfDiscoverNodes++;

            if(minDistance >  stats->distances[v] && stats->distances[v] != 0)
                minDistance = stats->distances[v];

            if(maxDistance < stats->distances[v])
                maxDistance = stats->distances[v];


        }

    }

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Min Dist", "Max Dist", " Discovered");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15u | \n", minDistance, maxDistance, numberOfDiscoverNodes);
    printf(" -----------------------------------------------------\n");


}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

void SSSPSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus, uint32_t delta)
{

    // The first subset, Ef, contains all edges (vi, vj) such that i < j; the second, Eb, contains edges (vi, vj) such that i > j.

    //calculate the size of each edge array
    uint32_t edgesPlusCounter = 0;
    uint32_t edgesMinusCounter = 0;
    // uint32_t numVerticesPlusCounter = 0;
    // uint32_t numVerticesMinusCounter = 0;
    uint32_t e;
    float weight;
    // uint32_t src;
    // uint32_t dest;

    #pragma omp parallel for private(e,weight) shared(graph,delta) reduction(+:edgesPlusCounter,edgesMinusCounter)
    for(e = 0 ; e < graph->num_edges ; e++)
    {

        // src  = graph->sorted_edges_array[e].src;
        // dest = graph->sorted_edges_array[e].dest;


        weight =  1;
#if WEIGHTED
        weight = graph->sorted_edges_array->edges_array_weight[e];
#endif



        if(weight > delta)
        {
            edgesPlusCounter++;

        }
        else if (weight <= delta)
        {
            edgesMinusCounter++;

        }
    }

    *graphPlus = graphCSRNew(graph->num_vertices, edgesPlusCounter, 1);
    *graphMinus =  graphCSRNew(graph->num_vertices, edgesMinusCounter, 1);

    struct EdgeList *edgesPlus = newEdgeList(edgesPlusCounter);
    struct EdgeList *edgesMinus = newEdgeList(edgesMinusCounter);

#if DIRECTED
    struct EdgeList *edgesPlusInverse = newEdgeList(edgesPlusCounter);
    struct EdgeList *edgesMinusInverse = newEdgeList(edgesMinusCounter);
#endif

    uint32_t edgesPlus_idx = 0;
    uint32_t edgesMinus_idx = 0;

    #pragma omp parallel for private(e,weight) shared(edgesMinus_idx,edgesPlus_idx, delta,edgesPlus,edgesMinus,graph)
    for(e = 0 ; e < graph->num_edges ; e++)
    {

        weight =  1;
#if WEIGHTED
        weight = graph->sorted_edges_array->edges_array_weight[e];
#endif
        uint32_t index = 0;

        if(weight > delta)
        {
            index = __sync_fetch_and_add(&edgesPlus_idx, 1);

            edgesPlus->edges_array_dest[index] = graph->sorted_edges_array->edges_array_dest[e];
            edgesPlus->edges_array_src[index] = graph->sorted_edges_array->edges_array_src[e];
#if WEIGHTED
            edgesPlus->edges_array_weight[index] = graph->sorted_edges_array->edges_array_weight[e];
#endif

#if DIRECTED
            edgesPlusInverse->edges_array_dest[index] = graph->sorted_edges_array->edges_array_src[e];
            edgesPlusInverse->edges_array_src[index] = graph->sorted_edges_array->edges_array_dest[e];
#if WEIGHTED
            edgesPlusInverse->edges_array_weight[index] = graph->sorted_edges_array->edges_array_weight[e];
#endif

#endif

        }
        else if (weight <= delta)
        {
            index = __sync_fetch_and_add(&edgesMinus_idx, 1);

            edgesMinus->edges_array_dest[index] = graph->sorted_edges_array->edges_array_dest[e];
            edgesMinus->edges_array_src[index] = graph->sorted_edges_array->edges_array_src[e];
#if WEIGHTED
            edgesMinus->edges_array_weight[index] = graph->sorted_edges_array->edges_array_weight[e];
#endif

#if DIRECTED
            edgesMinusInverse->edges_array_dest[index] = graph->sorted_edges_array->edges_array_src[e];
            edgesMinusInverse->edges_array_src[index] = graph->sorted_edges_array->edges_array_dest[e];
#if WEIGHTED
            edgesMinusInverse->edges_array_weight[index] = graph->sorted_edges_array->edges_array_weight[e];
#endif

#endif

        }
    }


    edgesPlus = sortRunAlgorithms(edgesPlus, 0);
    edgesMinus = sortRunAlgorithms(edgesMinus, 0);

#if DIRECTED
    edgesPlusInverse = sortRunAlgorithms(edgesPlusInverse, 0);
    edgesMinusInverse = sortRunAlgorithms(edgesMinusInverse, 0);
#endif

    graphCSRAssignEdgeList ((*graphPlus), edgesPlus, 0);
    graphCSRAssignEdgeList ((*graphMinus), edgesMinus, 0);

#if DIRECTED
    graphCSRAssignEdgeList ((*graphPlus), edgesPlusInverse, 1);
    graphCSRAssignEdgeList ((*graphMinus), edgesMinusInverse, 1);
#endif


}

int test_13(int x)
{
    int *i;
    i = &x;
    *i = x + 1;
    return x;
}

struct SSSPStats *SSSPGraphCSR(uint32_t source,  uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph, uint32_t delta)
{

    struct SSSPStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // push
        stats = SSSPDataDrivenPushGraphCSR(source, iterations, graph, delta);
        // SSSPDataDrivenPullGraphCSR(source, iterations, graph, delta); BUGGY
        break;
    case 1: // push
        stats = SSSPDataDrivenSplitPushGraphCSR(source, iterations, graph, delta);
        break;
    default:// push
        stats = SSSPDataDrivenPushGraphCSR(source, iterations, graph, delta);
        break;
    }

    return stats;

}

// struct SSSPStats *SSSPDataDrivenPullGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta)
// {

//     uint32_t v;
//     uint32_t iter = 0;
//     iterations = graph->num_vertices - 1;


//     struct SSSPStats *stats = (struct SSSPStats *) malloc(sizeof(struct SSSPStats));

//     stats->bucket_counter = 0;
//     stats->delta = delta;
//     stats->bucket_current = 0;
//     stats->processed_nodes = 0;
//     stats->buckets_total = 0;
//     stats->time_total = 0.0;
//     stats->num_vertices = graph->num_vertices;

//     struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
//     struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

//     struct Bitmap *bitmapSetCurr = newBitmap(graph->num_vertices);

//     uint32_t activeVertices = 0;


//     stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
//     stats->parents = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
//     stats->buckets_map = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));


//     struct GraphCSR *graphHeavy = NULL;
//     struct GraphCSR *graphLight = NULL;

//     printf(" -----------------------------------------------------\n");
//     printf("| %-51s | \n", "Starting Delta-Stepping Algorithm Pull DD (Source)");
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51u | \n", source);
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51s | \n", "Start Split Heavy/Light");
//     printf(" -----------------------------------------------------\n");
//     Start(timer_inner);
//     SSSPSpiltGraphCSR(graph, &graphHeavy, &graphLight, stats->delta);
//     Stop(timer_inner);
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51s | \n", "Graph Light Edges (Number)");
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51u | \n", graphLight->num_edges );
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51s | \n", "Graph Heavy Edges (Number)");
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51u | \n", graphHeavy->num_edges);
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51s | \n", "END Split Heavy/Light");
//     printf(" -----------------------------------------------------\n");
//     printf("| %-51f | \n",  Seconds(timer_inner));
//     printf(" -----------------------------------------------------\n");


//     printf(" -----------------------------------------------------\n");
//     printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active vertices", "Time (Seconds)");
//     printf(" -----------------------------------------------------\n");

//     if(source > graph->num_vertices)
//     {
//         printf(" -----------------------------------------------------\n");
//         printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
//         printf(" -----------------------------------------------------\n");
//         return NULL;
//     }

//     Start(timer);

//     Start(timer_inner);
//     //order vertices according to degree
//     #pragma omp parallel for
//     for(v = 0; v < graph->num_vertices; v++)
//     {

//         stats->buckets_map[v] = UINT_MAX / 2;
//         stats->distances[v] = UINT_MAX / 2;
//         stats->parents[v] = UINT_MAX;

//     }

//     stats->parents[source] = source;
//     stats->distances[source] = 0;

//     stats->buckets_map[source] = 0; // maps to bucket zero
//     stats->bucket_counter = 1;
//     stats->buckets_total = 1;
//     stats->bucket_current = 0;

//     activeVertices = 1;

//     uint32_t degree = graph->vertices[source].out_degree;
//     uint32_t edge_idx = graph->vertices[source].edges_idx;

//     for(v = edge_idx ; v < (edge_idx + degree) ; v++)
//     {

//         uint32_t t = graph->sorted_edges_array[v].dest;
//         stats->buckets_map[t] = stats->bucket_current;
//         stats->parents[t] = source;
//         setBitAtomic(bitmapSetCurr, t);
//         // activeVertices++;
//         stats->buckets_total++;

//         // printf("tb : %u v: %u cb %u \n",stats->buckets_total,  t, stats->bucket_current);

//     }

//     Stop(timer_inner);

//     printf("| %-15s | %-15u | %-15f | \n", "Init", stats->buckets_total,  Seconds(timer_inner));
//     printf(" -----------------------------------------------------\n");


//     while (stats->buckets_total)
//     {
//         // Start(timer_inner);
//         stats->processed_nodes += activeVertices;
//         activeVertices = 0;
//         stats->bucket_counter = 1;
//         clearBitmap(bitmapSetCurr);

//         while(stats->bucket_counter)
//         {
//             Start(timer_inner);
//             stats->bucket_counter = 0;
//             // uint32_t buckets_total_local =
//             // process light edges
//             #pragma omp parallel for private(v) shared(bitmapSetCurr, graphLight, stats) reduction(+ : activeVertices)
//             for(v = 0; v < graphLight->num_vertices; v++)
//             {

//                 uint32_t minDistance = UINT_MAX / 2;
//                 uint32_t degree;
//                 uint32_t j, u, w;
//                 uint32_t edge_idx;
//                 uint32_t src = UINT_MAX;
//                 struct Edge *edge = (struct Edge *) my_malloc(sizeof(struct Edge));

//                 if(__sync_bool_compare_and_swap(&(stats->buckets_map[v]), stats->bucket_current, (UINT_MAX / 2)))
//                 {
//                     // pop vertex from bucket list
//                     setBitAtomic(bitmapSetCurr, v);

//                     #pragma omp atomic update
//                     stats->buckets_total--;


//                     // printf("light tb : %u v: %u cb %u \n",stats->buckets_total,  v, stats->bucket_current);

//                     degree = graphLight->inverse_vertices[v].out_degree;
//                     edge_idx = graphLight->inverse_vertices[v].edges_idx;

//                     for(j = edge_idx ; j < (edge_idx + degree) ; j++)
//                     {
//                         u = graphLight->inverse_sorted_edges_array[j].dest;
//                         w = graphLight->inverse_sorted_edges_array[j].weight;

//                         if (minDistance > (stats->distances[u] + w))
//                         {
//                             minDistance = (stats->distances[u] + w);
//                             src = u;

//                             dest = v;
//                             weight = w;
//                             src = u;
//                         }
//                     }

//                     if(src != UINT_MAX)
//                         if(SSSPAtomicRelax(edge, stats))
//                         {


//                             // printf("relax tb : %u v: %u u: %u dis %u \n",stats->buckets_total,  dest, src,stats->distances[dest]);

//                             degree = graph->vertices[v].out_degree;
//                             edge_idx = graph->vertices[v].edges_idx;

//                             for(j = edge_idx ; j < (edge_idx + degree) ; j++)
//                             {
//                                 u = graph->sorted_edges_array[j].dest;
//                                 w = graph->sorted_edges_array[j].weight;

//                                 uint32_t oldBucket = stats->buckets_map[u];
//                                 uint32_t newBucket = stats->bucket_current;

//                                 if(__sync_bool_compare_and_swap(&(stats->buckets_map[u]), oldBucket, newBucket))
//                                 {
//                                     if(oldBucket == UINT_MAX / 2)
//                                         #pragma omp atomic update
//                                         stats->buckets_total++;

//                                     stats->bucket_counter = 1;
//                                 }
//                             }
//                             activeVertices++;
//                         }

//                 }
//             }

//             Stop(timer_inner);

//             if(activeVertices)
//                 printf("| L%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));
//         }

//         Start(timer_inner);

//         #pragma omp parallel for private(v) shared(bitmapSetCurr, graphHeavy, stats) reduction(+ : activeVertices)
//         for(v = 0; v < graphHeavy->num_vertices; v++)
//         {

//             uint32_t minDistance = UINT_MAX / 2;
//             uint32_t degree;
//             uint32_t j, u, w;
//             uint32_t edge_idx;
//             uint32_t src = UINT_MAX;
//             struct Edge *edge = (struct Edge *) my_malloc(sizeof(struct Edge));

//             if(getBit(bitmapSetCurr, v))
//             {

//                 // printf("heavy tb : %u v: %u cb %u \n",stats->buckets_total,  v, stats->bucket_current);

//                 degree = graphHeavy->inverse_vertices[v].out_degree;
//                 edge_idx = graphHeavy->inverse_vertices[v].edges_idx;

//                 for(j = edge_idx ; j < (edge_idx + degree) ; j++)
//                 {
//                     u = graphHeavy->inverse_sorted_edges_array[j].dest;
//                     w = graphHeavy->inverse_sorted_edges_array[j].weight;

//                     if (minDistance > (stats->distances[u] + w))
//                     {
//                         minDistance = (stats->distances[u] + w);
//                         src = u;

//                         dest = v;
//                         weight = w;
//                         src = u;
//                     }
//                 }


//                 if(src != UINT_MAX)
//                     if(SSSPAtomicRelax(edge, stats))
//                     {

//                         degree = graph->vertices[v].out_degree;
//                         edge_idx = graph->vertices[v].edges_idx;

//                         // printf("relax tb : %u v: %u u: %u dis %u \n",stats->buckets_total,  dest, src,stats->distances[dest]);


//                         // for(j = edge_idx ; j < (edge_idx + degree) ; j++){
//                         //      u = graph->sorted_edges_array[j].dest;
//                         //      w = graph->sorted_edges_array[j].weight;

//                         //      uint32_t oldBucket = stats->buckets_map[u];
//                         //      uint32_t newBucket;

//                         //      if(stats->bucket_counter)
//                         //       newBucket = stats->bucket_current;
//                         //      else
//                         //       newBucket = stats->bucket_current+1;

//                         //      if(__sync_bool_compare_and_swap(&(stats->buckets_map[u]), oldBucket, newBucket)){
//                         //          if(oldBucket == UINT_MAX/2)
//                         //  #pragma omp atomic update
//                         //  stats->buckets_total++;

//                         //      }
//                         //  }
//                         activeVertices++;
//                     }

//             }
//         }

//         iter++;
//         stats->bucket_current++;
//         // clearBitmap(bitmapSetCurr);
//         Stop(timer_inner);
//         if(activeVertices)
//             printf("| H%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));


//     }



//     Stop(timer);
//     stats->time_total += Seconds(timer);
//     printf(" -----------------------------------------------------\n");
//     printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
//     printf(" -----------------------------------------------------\n");
//     SSSPPrintStatsDetails(stats);

//     // free resources
//     free(timer);
//     free(timer_inner);
//     freeBitmap(bitmapSetCurr);
//     graphCSRFree(graphHeavy);
//     graphCSRFree(graphLight);


//     // SSSPPrintStats(stats);
//     return stats;
// }


struct SSSPStats *SSSPDataDrivenPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta)
{

    uint32_t v;
    uint32_t iter = 0;

    iterations = graph->num_vertices - 1;

    struct SSSPStats *stats = newSSSPStatsGraphCSR(graph, delta);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Delta-Stepping Algorithm Push DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct Bitmap *bitmapSetCurr = newBitmap(graph->num_vertices);

    uint32_t activeVertices = 0;



    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active vertices", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");



    Start(timer);

    Start(timer_inner);
    //order vertices according to degree


    stats->parents[source] = source;
    stats->distances[source] = 0;

    stats->buckets_map[source] = 0; // maps to bucket zero
    stats->bucket_counter = 1;
    stats->buckets_total = 1;
    stats->bucket_current = 0;

    activeVertices = 1;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", stats->buckets_total,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");


    while (stats->buckets_total)
    {
        // Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;
        stats->bucket_counter = 1;
        clearBitmap(bitmapSetCurr);
        Start(timer_inner);

        while(stats->bucket_counter)
        {
            Start(timer_inner);
            stats->bucket_counter = 0;
            // uint32_t buckets_total_local =
            // process light edges
            #pragma omp parallel for private(v) shared(bitmapSetCurr, graph, stats) reduction(+ : activeVertices)
            for(v = 0; v < graph->num_vertices; v++)
            {

                if(__sync_bool_compare_and_swap(&(stats->buckets_map[v]), stats->bucket_current, (UINT_MAX / 2)))
                {
                    // if(stats->buckets_map[v] == stats->bucket_current) {

                    // pop vertex from bucket list
                    setBitAtomic(bitmapSetCurr, v);

                    #pragma omp atomic update
                    stats->buckets_total--;

                    // stats->buckets_map[v] = UINT_MAX/2;

                    uint32_t degree = graph->vertices->out_degree[v];
                    uint32_t edge_idx = graph->vertices->edges_idx[v];
                    uint32_t j;
                    for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                    {
                        uint32_t src = graph->sorted_edges_array->edges_array_src[j];
                        uint32_t dest = graph->sorted_edges_array->edges_array_dest[j];
                        float weight = 1;
#if WEIGHTED
                        weight = graph->sorted_edges_array->edges_array_weight[j];
#endif

                        if(numThreads == 1)
                            activeVertices += SSSPRelax(src, dest, weight, stats);
                        else
                            activeVertices += SSSPAtomicRelax(src, dest, weight, stats);
                    }

                }
            }

            Stop(timer_inner);

            if(activeVertices)
                printf("| L%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));
        }

        iter++;
        stats->bucket_current++;
        Stop(timer_inner);
        if(activeVertices)
            printf("| H%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));
    }



    Stop(timer);
    stats->time_total += Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    SSSPPrintStatsDetails(stats);

    // free resources
    free(timer);
    free(timer_inner);
    freeBitmap(bitmapSetCurr);


    // SSSPPrintStats(stats);
    return stats;
}


struct SSSPStats *SSSPDataDrivenSplitPushGraphCSR(uint32_t source,  uint32_t iterations, struct GraphCSR *graph, uint32_t delta)
{

    uint32_t v;
    uint32_t iter = 0;

    iterations = graph->num_vertices - 1;

    struct SSSPStats *stats = newSSSPStatsGraphCSR(graph, delta);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Delta-Stepping Algorithm Push DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct Bitmap *bitmapSetCurr = newBitmap(graph->num_vertices);

    uint32_t activeVertices = 0;




    struct GraphCSR *graphHeavy = NULL;
    struct GraphCSR *graphLight = NULL;


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Start Split Heavy/Light");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    SSSPSpiltGraphCSR(graph, &graphHeavy, &graphLight, stats->delta);
    Stop(timer_inner);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Graph Light Edges (Number)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", graphLight->num_edges );
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Graph Heavy Edges (Number)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", graphHeavy->num_edges);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "END Split Heavy/Light");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n",  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");


    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active vertices", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");



    Start(timer);

    Start(timer_inner);
    //order vertices according to degree


    stats->parents[source] = source;
    stats->distances[source] = 0;

    stats->buckets_map[source] = 0; // maps to bucket zero
    stats->bucket_counter = 1;
    stats->buckets_total = 1;
    stats->bucket_current = 0;

    activeVertices = 1;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", stats->buckets_total,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");


    while (stats->buckets_total)
    {
        // Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;
        stats->bucket_counter = 1;
        clearBitmap(bitmapSetCurr);

        while(stats->bucket_counter)
        {
            Start(timer_inner);
            stats->bucket_counter = 0;
            // uint32_t buckets_total_local =
            // process light edges
            #pragma omp parallel for private(v) shared(bitmapSetCurr, graphLight, stats) reduction(+ : activeVertices)
            for(v = 0; v < graphLight->num_vertices; v++)
            {

                if(__sync_bool_compare_and_swap(&(stats->buckets_map[v]), stats->bucket_current, (UINT_MAX / 2)))
                {
                    // if(stats->buckets_map[v] == stats->bucket_current) {

                    // pop vertex from bucket list
                    setBitAtomic(bitmapSetCurr, v);

                    #pragma omp atomic update
                    stats->buckets_total--;

                    // stats->buckets_map[v] = UINT_MAX/2;

                    uint32_t degree = graphLight->vertices->out_degree[v];
                    uint32_t edge_idx = graphLight->vertices->edges_idx[v];
                    uint32_t j;
                    for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                    {
                        uint32_t src = graphLight->sorted_edges_array->edges_array_src[j];
                        uint32_t dest = graphLight->sorted_edges_array->edges_array_dest[j];
                        float weight = 1;
#if WEIGHTED
                        weight = graphLight->sorted_edges_array->edges_array_weight[j];
#endif

                        if(numThreads == 1)
                            activeVertices += SSSPRelax(src, dest, weight, stats);
                        else
                            activeVertices += SSSPAtomicRelax(src, dest, weight, stats);
                    }

                }
            }

            Stop(timer_inner);

            if(activeVertices)
                printf("| L%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));
        }

        Start(timer_inner);

        #pragma omp parallel for private(v) shared(bitmapSetCurr, graphHeavy, stats) reduction(+ : activeVertices)
        for(v = 0; v < graphHeavy->num_vertices; v++)
        {
            if(getBit(bitmapSetCurr, v))
            {

                uint32_t degree = graphHeavy->vertices->out_degree[v];
                uint32_t edge_idx = graphHeavy->vertices->edges_idx[v];
                uint32_t j;


                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    uint32_t src = graphHeavy->sorted_edges_array->edges_array_src[j];
                    uint32_t dest = graphHeavy->sorted_edges_array->edges_array_dest[j];
                    float weight = 1;
#if WEIGHTED
                    weight = graphHeavy->sorted_edges_array->edges_array_weight[j];
#endif

                    if(numThreads == 1)
                        activeVertices += SSSPRelax(src, dest, weight, stats);
                    else
                        activeVertices += SSSPAtomicRelax(src, dest, weight, stats);
                }
            }
        }

        iter++;
        stats->bucket_current++;
        Stop(timer_inner);
        if(activeVertices)
            printf("| H%-14u | %-15u | %-15f |\n", iter, stats->buckets_total, Seconds(timer_inner));
    }



    Stop(timer);
    stats->time_total += Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    SSSPPrintStatsDetails(stats);

    // free resources
    free(timer);
    free(timer_inner);
    freeBitmap(bitmapSetCurr);
    graphCSRFree(graphHeavy);
    graphCSRFree(graphLight);


    // SSSPPrintStats(stats);
    return stats;
}





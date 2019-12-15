// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : bellmanFord.c
// Create : 2019-09-28 15:19:50
// Revise : 2019-09-28 15:34:29
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <limits.h> //UINT_MAX

#include "myMalloc.h"
#include "timer.h"
#include "mt19937.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"

#include "graphConfig.h"
#include "sortRun.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "bellmanFord.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct BellmanFordStats *newBellmanFordStatsGraphCSR(struct GraphCSR *graph)
{
    __u32 v;
    struct BellmanFordStats *stats = (struct BellmanFordStats *) my_malloc(sizeof(struct BellmanFordStats));
    stats->processed_nodes = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;
    stats->distances  = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));
    stats->parents = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {

        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;

    }

    return stats;
}
struct BellmanFordStats *newBellmanFordStatsGraphGrid(struct GraphGrid *graph)
{

    __u32 v;

    struct BellmanFordStats *stats = (struct BellmanFordStats *) my_malloc(sizeof(struct BellmanFordStats));
    stats->processed_nodes = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;
    stats->distances  = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));
    stats->parents = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {

        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;

    }

    return stats;


}
struct BellmanFordStats *newBellmanFordStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    __u32 v;
    struct BellmanFordStats *stats = (struct BellmanFordStats *) my_malloc(sizeof(struct BellmanFordStats));
    stats->processed_nodes = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;
    stats->distances  = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));
    stats->parents = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;
    }

    return stats;

}
struct BellmanFordStats *newBellmanFordStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    __u32 v;
    struct BellmanFordStats *stats = (struct BellmanFordStats *) my_malloc(sizeof(struct BellmanFordStats));
    stats->processed_nodes = 0;
    stats->time_total = 0.0;
    stats->num_vertices = graph->num_vertices;
    stats->distances  = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));
    stats->parents = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {

        stats->distances[v] = UINT_MAX / 2;
        stats->parents[v] = UINT_MAX;

    }

    return stats;


}

void freeBellmanFordStats(struct BellmanFordStats *stats)
{
    if(stats)
    {
        if(stats->distances)
            free(stats->distances);
        if(stats->parents)
            free(stats->parents);
        free(stats);
    }
}

// ********************************************************************************************
// ***************                  Auxiliary functions                          **************
// ********************************************************************************************

__u32 bellmanFordAtomicMin(__u32 *dist, __u32 newValue)
{

    __u32 oldValue;
    __u32 flag = 0;

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

__u32 bellmanFordCompareDistanceArrays(struct BellmanFordStats *stats1, struct BellmanFordStats *stats2)
{

    __u32 v = 0;


    for(v = 0 ; v < stats1->num_vertices ; v++)
    {

        if(stats1->distances[v] != stats2->distances[v])
        {

            return 0;
        }
        // else if(stats1->distances[v] != UINT_MAX/2)


    }

    return 1;

}

int bellmanFordAtomicRelax(__u32 src, __u32 dest, __u32 weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext)
{
    // __u32 oldParent, newParent;
    __u32 oldDistanceV = UINT_MAX / 2;
    __u32 oldDistanceU = UINT_MAX / 2;
    __u32 newDistance = UINT_MAX / 2;
    __u32 flagu = 0;
    __u32 flagv = 0;
    // __u32 flagp = 0;
    __u32 activeVertices = 0;

    do
    {

        flagu = 0;
        flagv = 0;
        // flagp = 0;


        oldDistanceV = stats->distances[src];
        oldDistanceU = stats->distances[dest];
        // oldParent = stats->parents[dest];
        newDistance = oldDistanceV + weight;

        if( oldDistanceU > newDistance )
        {

            // newParent = src;
            newDistance = oldDistanceV + weight;

            if(__sync_bool_compare_and_swap(&(stats->distances[src]), oldDistanceV, oldDistanceV))
            {
                flagv = 1;
            }

            if(__sync_bool_compare_and_swap(&(stats->distances[dest]), oldDistanceU, newDistance) && flagv)
            {
                flagu = 1;
                setBitAtomic(bitmapNext, dest);
                activeVertices++;
                stats->parents[dest] = src;
            }

            // if(__sync_bool_compare_and_swap(&(stats->parents[dest]), oldParent, newParent) && flagv && flagu)
            // {
            //     flagp = 1;
            // }

            // if(!getBit(bitmapNext, dest) && flagv && flagu && flagp)
            // {
            //     setBitAtomic(bitmapNext, dest);
            //     activeVertices++;
            // }

        }
        else
        {
            return activeVertices;
        }

    }
    while (!flagu || !flagv );


    return activeVertices;

}



int bellmanFordRelax(__u32 src, __u32 dest, __u32 weight, struct BellmanFordStats *stats, struct Bitmap *bitmapNext)
{

    __u32 activeVertices = 0;
    __u32 newDistance = stats->distances[src] + weight;

    if( stats->distances[dest] > newDistance )
    {


        stats->distances[dest] = newDistance;
        stats->parents[dest] = src;

        if(!getBit(bitmapNext, dest))
        {
            activeVertices++;
            setBit(bitmapNext, dest);
        }
    }


    return activeVertices;

}

void bellmanFordPrintStats(struct BellmanFordStats *stats)
{
    __u32 v;
    __u32 sum = 0;
    for(v = 0; v < stats->num_vertices; v++)
    {

        if(stats->distances[v] != UINT_MAX / 2)
        {
            sum += stats->distances[v];
            printf("p %u d %u \n", stats->parents[v], stats->distances[v]);

        }


    }

    printf("sum %u \n", sum);
}

void bellmanFordPrintStatsDetails(struct BellmanFordStats *stats)
{
    __u32 v;
    __u32 minDistance = UINT_MAX / 2;
    __u32 maxDistance = 0;
    __u32 numberOfDiscoverNodes = 0;


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


// -- To shuffle an array a of n elements (indices 0..n-1):
// for i from 0 to n−2 do
//      j ← random integer such that i ≤ j < n
//      exchange a[i] and a[j]


// used with Bannister, M. J.; Eppstein, D. (2012). Randomized speedup of the Bellman–Ford algorithm

void durstenfeldShuffle(__u32 *vertices, __u32 size)
{

    __u32 v;
    for(v = 0; v < size; v++)
    {

        __u32 idx = (generateRandInt(mt19937var) % (size - 1));
        __u32 temp = vertices[v];
        vertices[v] = vertices[idx];
        vertices[idx] = temp;

    }

}


// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphGrid(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphGrid *graph)
{

    struct BellmanFordStats *stats;

    switch (pushpull)
    {
    case 0: // pull
        stats = bellmanFordPullRowGraphGrid(source, iterations, graph);
        break;
    case 1: // push
        stats = bellmanFordPushColumnGraphGrid(source, iterations, graph);
        break;
    default:// push
        stats = bellmanFordPushColumnGraphGrid(source, iterations, graph);
        break;
    }

    return stats;

}
struct BellmanFordStats *bellmanFordPullRowGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph)
{

    __u32 iter = 0;
    __u32 totalPartitions  = graph->grid->num_partitions;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm ROW-WISE DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }





    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;





    Start(timer);

    Start(timer_inner);
    //order vertices according to degree


    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;

        __u32 i;
        #pragma omp parallel for private(i) reduction(+ : activeVertices) schedule (dynamic,numThreads)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            __u32 j;
            // #pragma omp parallel for private(j) reduction(+ : activeVertices) schedule (dynamic,8)
            for (j = 0; j < totalPartitions; ++j)
            {
                __u32 k;

                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {

                    __u32 src = partition->edgeList->edges_array_src[k];
                    __u32 dest = partition->edgeList->edges_array_dest[k];
                    __u32 weight = partition->edgeList->edges_array_weight[k];

                    if(getBit(bitmapCurr, src))
                    {
                        if(numThreads == 1)
                            activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                        else
                            activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                    }
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total += Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);


    // bellmanFordPrintStats(stats);
    return stats;


}
struct BellmanFordStats *bellmanFordPushColumnGraphGrid(__u32 source,  __u32 iterations, struct GraphGrid *graph)
{

    __u32 iter = 0;
    __u32 totalPartitions  = graph->grid->num_partitions;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm COL-WISE DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;




    Start(timer);
    Start(timer_inner);


    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;

        __u32 j;
        #pragma omp parallel for private(j) reduction(+ : activeVertices) schedule (dynamic,numThreads)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
        {
            __u32 i;
            // #pragma omp parallel for private(i) reduction(+ : activeVertices) schedule (dynamic,numThreads)
            for (i = 0; i < totalPartitions; ++i)
            {
                __u32 k;

                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {

                    __u32 src = partition->edgeList->edges_array_src[k];
                    __u32 dest = partition->edgeList->edges_array_dest[k];
                    __u32 weight = partition->edgeList->edges_array_weight[k];

                    if(getBit(bitmapCurr, src))
                    {
                        // if(numThreads == 1)
                        activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                        // else
                        // activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                    }
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);

    // bellmanFordPrintStats(stats);
    return stats;


}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

void bellmanFordSpiltGraphCSR(struct GraphCSR *graph, struct GraphCSR **graphPlus, struct GraphCSR **graphMinus)
{

    // The first subset, Ef, contains all edges (vi, vj) such that i < j; the second, Eb, contains edges (vi, vj) such that i > j.

    //calculate the size of each edge array
    __u32 edgesPlusCounter = 0;
    __u32 edgesMinusCounter = 0;
    __u32 e;
    __u32 src;
    __u32 dest;

    #pragma omp parallel for private(e,src,dest) shared(graph) reduction(+:edgesPlusCounter,edgesMinusCounter)
    for(e = 0 ; e < graph->num_edges ; e++)
    {

        src  = graph->sorted_edges_array->edges_array_src[e];
        dest = graph->sorted_edges_array->edges_array_dest[e];
        if(src <= dest)
        {
            edgesPlusCounter++;
        }
        else if (src > dest)
        {
            edgesMinusCounter++;
        }
    }

    *graphPlus = graphCSRNew(graph->num_vertices, edgesPlusCounter, 1);
    *graphMinus =  graphCSRNew(graph->num_vertices, edgesMinusCounter, 1);

    struct EdgeList *edgesPlus = newEdgeList(edgesPlusCounter);
    struct EdgeList *edgesMinus = newEdgeList(edgesMinusCounter);


    edgesPlus->num_vertices = graph->num_vertices;
    edgesMinus->num_vertices = graph->num_vertices;

    __u32 edgesPlus_idx = 0;
    __u32 edgesMinus_idx = 0;

    #pragma omp parallel for private(e,src,dest) shared(edgesMinus_idx,edgesPlus_idx, edgesPlus,edgesMinus,graph)
    for(e = 0 ; e < graph->num_edges ; e++)
    {
        __u32 localEdgesPlus_idx = 0;
        __u32 localEdgesMinus_idx = 0;

        src  = graph->sorted_edges_array->edges_array_src[e];
        dest = graph->sorted_edges_array->edges_array_dest[e];
        if(src <= dest)
        {
            localEdgesPlus_idx = __sync_fetch_and_add(&edgesPlus_idx, 1);

            edgesPlus->edges_array_src[localEdgesPlus_idx] = graph->sorted_edges_array->edges_array_src[e];
            edgesPlus->edges_array_dest[localEdgesPlus_idx] = graph->sorted_edges_array->edges_array_dest[e];
            edgesPlus->edges_array_weight[localEdgesPlus_idx] = graph->sorted_edges_array->edges_array_weight[e];
        }
        else if (src > dest)
        {
            localEdgesMinus_idx = __sync_fetch_and_add(&edgesMinus_idx, 1);

            edgesMinus->edges_array_src[localEdgesMinus_idx] = graph->sorted_edges_array->edges_array_src[e];
            edgesMinus->edges_array_dest[localEdgesMinus_idx] = graph->sorted_edges_array->edges_array_dest[e];
            edgesMinus->edges_array_weight[localEdgesMinus_idx] = graph->sorted_edges_array->edges_array_weight[e];

        }
    }



    edgesPlus = sortRunAlgorithms(edgesPlus, 0);
    edgesMinus = sortRunAlgorithms(edgesMinus, 0);

    graphCSRAssignEdgeList ((*graphPlus), edgesPlus, 0);

    printf("\n");

    graphCSRAssignEdgeList ((*graphMinus), edgesMinus, 0);


}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

void printDistances(struct BellmanFordStats *stats)
{

    __u32 vertex_id;
    __u32 sum = 0;
    for(vertex_id = 0; vertex_id < stats->num_vertices ; vertex_id++)
    {
        sum += stats->distances[vertex_id];
        printf("v: %u d: %u \n", vertex_id, sum);
    }

}

struct BellmanFordStats *bellmanFordGraphCSR(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphCSR *graph)
{

    struct BellmanFordStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = bellmanFordDataDrivenPullGraphCSR(source, iterations, graph);
        break;
    case 1: // push
        stats = bellmanFordDataDrivenPushGraphCSR(source, iterations, graph);
        break;
    case 2: // randomized push
        stats = bellmanFordRandomizedDataDrivenPushGraphCSR(source, iterations, graph);
        break;
    default:// push
        stats = bellmanFordDataDrivenPushGraphCSR(source, iterations, graph);
        break;
    }

    return stats;
}

struct BellmanFordStats *bellmanFordDataDrivenPullGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph)
{


    __u32 v;
    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Pull DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;

    struct Vertex *vertices = NULL;
    struct EdgeList  *sorted_edges_array = NULL;

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array;
#endif


    Start(timer);

    Start(timer_inner);

    setBit(bitmapNext, source);
    bitmapNext->numSetBits++;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    __u32 degree = graph->vertices->out_degree[source];
    __u32 edge_idx = graph->vertices->edges_idx[source];

    for(v = edge_idx ; v < (edge_idx + degree) ; v++)
    {

        __u32 t = graph->sorted_edges_array->edges_array_dest[v];
        stats->parents[t] = source;
        bitmapNext->numSetBits++;
        setBit(bitmapNext, t);
        activeVertices++;

    }

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(v) shared(vertices,sorted_edges_array,graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            __u32 minDistance = UINT_MAX / 2;
            __u32 minParent = UINT_MAX;
            __u32 degree;
            __u32 j, u, w;
            __u32 edge_idx;

            if(getBit(bitmapCurr, v))
            {

                degree = vertices->out_degree[v];
                edge_idx = vertices->edges_idx[v];
                // printf("degree %u source %u \n",degree,v );
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = sorted_edges_array->edges_array_dest[j];
                    w = sorted_edges_array->edges_array_weight[j];


                    if (minDistance > (stats->distances[u] + w))
                    {
                        minDistance = (stats->distances[u] + w);
                        minParent = u;
                    }
                }

                if(bellmanFordAtomicMin(&(stats->distances[v]), minDistance))
                {
                    stats->parents[v] = minParent;

                    degree = graph->vertices->out_degree[v];
                    edge_idx = graph->vertices->edges_idx[v];

                    for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                    {
                        u = graph->sorted_edges_array->edges_array_dest[j];

                        if(!getBit(bitmapNext, u))
                        {
                            activeVertices++;
                            setBitAtomic(bitmapNext, u);
                        }
                    }
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);



        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);

    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);



    // bellmanFordPrintStats(stats);
    return stats;


}



struct BellmanFordStats *bellmanFordDataDrivenPushGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph)
{

    __u32 v;


    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Push DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;




    Start(timer);

    Start(timer_inner);
    //order vertices according to degree

    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(v) shared(graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            if(getBit(bitmapCurr, v))
            {

                __u32 degree = graph->vertices->out_degree[v];
                __u32 edge_idx = graph->vertices->edges_idx[v];
                __u32 j;
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    __u32 src = graph->sorted_edges_array->edges_array_src[j];
                    __u32 dest = graph->sorted_edges_array->edges_array_dest[j];
                    __u32 weight = graph->sorted_edges_array->edges_array_weight[j];

                    if(numThreads == 1)
                        activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                    else
                        activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);


    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);



    // bellmanFordPrintStats(stats);
    return stats;
}


// Randomized Speedup of the Bellman–Ford Algorithm

// number the vertices randomly such that all permutations with s first are equally likely
//  C ← {s}
//      while C 6= ∅ do
//          for each vertex u in numerical order do
//              if u ∈ C or D[v] has changed since start of iteration then
//                  for each edge uv in graph G+ do
//                      relax(u, v)
// for each vertex u in reverse numerical order do
//  if u ∈ C or D[v] has changed since start of iteration then
//      for each edge uv in graph G− do
//          relax(u, v)
// C ← {vertices v for which D[v] changed}

struct BellmanFordStats *bellmanFordRandomizedDataDrivenPushGraphCSR(__u32 source,  __u32 iterations, struct GraphCSR *graph)
{

    __u32 v;

    __u32 n;
    __u32 *vertices;
    __u32 *degrees;
    __u32 iter = 0;
    struct GraphCSR *graphPlus = NULL;
    struct GraphCSR *graphMinus = NULL;

    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Push DD");
    printf("| %-51s | \n", "Randomized G+/G- optimization (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;


    vertices = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));
    degrees = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Start Split G+/G-");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    bellmanFordSpiltGraphCSR(graph, &graphPlus, &graphMinus);
    Stop(timer_inner);
    printf("| %-51f | \n",  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");




    Start(timer);



    Start(timer_inner);

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vertices[v] = v;
        degrees[v] = graph->vertices->out_degree[v];
    }

    //randomize iteratiing accross verticess

    durstenfeldShuffle(vertices, graph->num_vertices);

    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(v,n) shared(vertices,graphPlus,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(n = 0; n < graphPlus->num_vertices; n++)
        {

            v = vertices[n];

            if(getBit(bitmapCurr, v))
            {

                __u32 degree = graphPlus->vertices->out_degree[v];
                __u32 edge_idx = graphPlus->vertices->edges_idx[v];
                __u32 j;
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {

                    __u32 src = graphPlus->sorted_edges_array->edges_array_src[j];
                    __u32 dest = graphPlus->sorted_edges_array->edges_array_dest[j];
                    __u32 weight = graphPlus->sorted_edges_array->edges_array_weight[j];

                    if(numThreads == 1)
                        activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                    else
                        activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                }
            }
        }

        #pragma omp parallel for private(v,n) shared(vertices,graphMinus,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(n = 0; n < graphMinus->num_vertices; n++)
        {

            v = vertices[n];

            if(getBit(bitmapCurr, v))
            {

                __u32 degree = graphMinus->vertices->out_degree[v];
                __u32 edge_idx = graphMinus->vertices->edges_idx[v];
                __u32 j;
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {

                    __u32 src = graphMinus->sorted_edges_array->edges_array_src[j];
                    __u32 dest = graphMinus->sorted_edges_array->edges_array_dest[j];
                    __u32 weight = graphMinus->sorted_edges_array->edges_array_weight[j];

                    if(numThreads == 1)
                        activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                    else
                        activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    free(vertices);
    free(degrees);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);

    // graphCSRFree(graphPlus);
    // graphCSRFree(graphMinus);


    // bellmanFordPrintStats(stats);
    return stats;
}

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjArrayList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph)
{

    struct BellmanFordStats *stats;

    switch (pushpull)
    {
    case 0: // push
        stats = bellmanFordDataDrivenPullGraphAdjArrayList(source, iterations, graph);
        break;
    case 1: // pull
        stats = bellmanFordDataDrivenPushGraphAdjArrayList(source, iterations, graph);
        break;
    default:// push
        stats = bellmanFordDataDrivenPushGraphAdjArrayList(source, iterations, graph);
        break;
    }
   
    return stats;
}


struct BellmanFordStats *bellmanFordDataDrivenPullGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    __u32 degree;
    __u32 v;

    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct EdgeList *nodes;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphAdjArrayList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Pull DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }



    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;





    Start(timer);

    Start(timer_inner);

    setBit(bitmapNext, source);
    bitmapNext->numSetBits++;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    nodes = graph->vertices[source].outNodes;
    degree = graph->vertices[source].out_degree;

    for(v = 0 ; v < (degree) ; v++)
    {

        __u32 t = nodes->edges_array_dest[v];
        stats->parents[t] = source;
        bitmapNext->numSetBits++;
        setBit(bitmapNext, t);
        activeVertices++;

    }

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;

        #pragma omp parallel for private(nodes,v) shared(graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            __u32 minDistance = UINT_MAX / 2;
            __u32 degree;
            __u32 j, u, w;

            __u32 minParent = UINT_MAX;

            if(getBit(bitmapCurr, v))
            {
#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
                nodes = graph->vertices[v].inNodes;
                degree = graph->vertices[v].in_degree;
#else
                nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;
#endif
                // printf("degree %u source %u \n",degree,v );
                for(j = 0 ; j < (degree) ; j++)
                {
                    u = nodes->edges_array_dest[j];
                    w = nodes->edges_array_weight[j];
                    // printf("w %u \n",w );
                    if (minDistance > (stats->distances[u] + w))
                    {
                        minDistance = (stats->distances[u] + w);
                        minParent = u;
                    }
                }

                if(bellmanFordAtomicMin(&(stats->distances[v]), minDistance))
                {
                    // stats->distances[v] = minDistance;
                    stats->parents[v] = minParent;
                    nodes = graph->vertices[v].outNodes;
                    degree = graph->vertices[v].out_degree;

                    for(j = 0 ; j < (degree) ; j++)
                    {
                        u = nodes->edges_array_dest[j];
                        w = nodes->edges_array_weight[j];

                        if(!getBit(bitmapNext, u))
                        {
                            activeVertices++;
                            setBitAtomic(bitmapNext, u);
                        }
                    }
                }
            }
        }

        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);



    // bellmanFordPrintStats(stats);
    return stats;

}
struct BellmanFordStats *bellmanFordDataDrivenPushGraphAdjArrayList(__u32 source,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    __u32 v;

    struct EdgeList *nodes;
    __u32 degree;
    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphAdjArrayList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Push DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }



    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;



    Start(timer);

    Start(timer_inner);


    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(nodes,degree,v) shared(graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            if(getBit(bitmapCurr, v))
            {

                degree = graph->vertices[v].out_degree;
                nodes = graph->vertices[v].outNodes;
                __u32 j;
                for(j = 0 ; j < (degree) ; j++)
                {

                    __u32 src = nodes->edges_array_src[j];
                    __u32 dest = nodes->edges_array_dest[j];
                    __u32 weight = nodes->edges_array_weight[j];

                    if(numThreads == 1)
                        activeVertices += bellmanFordRelax(src, dest, weight, stats, bitmapNext);
                    else
                        activeVertices += bellmanFordAtomicRelax(src, dest, weight, stats, bitmapNext);
                }

            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total += Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);


    // bellmanFordPrintStats(stats);
    return stats;

}

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct BellmanFordStats *bellmanFordGraphAdjLinkedList(__u32 source,  __u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph)
{

    struct BellmanFordStats *stats;

    switch (pushpull)
    {
    case 0: // pull
        stats = bellmanFordPullGraphAdjLinkedList(source, iterations, graph);
        break;
    case 1: // push
        stats = bellmanFordPushGraphAdjLinkedList(source, iterations, graph);
        break;
    default:// push
        stats = bellmanFordPushGraphAdjLinkedList(source, iterations, graph);
        break;
    }
    
    return stats;

}

struct BellmanFordStats *bellmanFordPullGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    __u32 degree;
    __u32 v;
    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphAdjLinkedList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Pull DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct AdjLinkedListNode *nodes;

    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;


    Start(timer);

    Start(timer_inner);

    setBit(bitmapNext, source);
    bitmapNext->numSetBits++;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    nodes = graph->vertices[source].outNodes;
    degree = graph->vertices[source].out_degree;

    for(v = 0 ; v < (degree) ; v++)
    {

        __u32 t = nodes->dest;
        nodes = nodes->next;
        stats->parents[t] = source;
        bitmapNext->numSetBits++;
        setBit(bitmapNext, t);
        activeVertices++;

    }

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(nodes,v) shared(graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            __u32 minDistance = UINT_MAX / 2;
            __u32 degree;
            __u32 j, u, w;


            if(getBit(bitmapCurr, v))
            {

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
                nodes = graph->vertices[v].inNodes;
                degree = graph->vertices[v].in_degree;
#else
                nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;
#endif

                for(j = 0 ; j < (degree) ; j++)
                {
                    u = nodes->dest;
                    w = nodes->weight;
                    nodes = nodes->next;

                    if (minDistance > (stats->distances[u] + w))
                    {
                        minDistance = (stats->distances[u] + w);
                    }
                }

                if(bellmanFordAtomicMin(&(stats->distances[v]), minDistance))
                {
                    stats->parents[v] = minDistance;

                    nodes = graph->vertices[v].outNodes;
                    degree = graph->vertices[v].out_degree;

                    for(j = 0 ; j < (degree) ; j++)
                    {
                        u = nodes->dest;
                        w = nodes->weight;
                        nodes = nodes->next;

                        if(!getBit(bitmapNext, u))
                        {
                            activeVertices++;
                            setBitAtomic(bitmapNext, u);
                        }
                    }
                }
            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);



        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);



    // bellmanFordPrintStats(stats);
    return stats;

}
struct BellmanFordStats *bellmanFordPushGraphAdjLinkedList(__u32 source,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    __u32 v;

    struct AdjLinkedListNode *nodes;
    __u32 degree;
    __u32 iter = 0;
    iterations = graph->num_vertices - 1;
    struct BellmanFordStats *stats = newBellmanFordStatsGraphAdjLinkedList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Bellman-Ford Algorithm Push DD (Source)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Active Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }



    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);
    int activeVertices = 0;





    Start(timer);

    Start(timer_inner);

    setBit(bitmapNext, source);
    bitmapNext->numSetBits = 1;
    stats->parents[source] = source;
    stats->distances[source] = 0;

    swapBitmaps(&bitmapCurr, &bitmapNext);
    clearBitmap(bitmapNext);
    activeVertices++;

    Stop(timer_inner);

    printf("| %-15s | %-15u | %-15f | \n", "Init", activeVertices,  Seconds(timer_inner));
    printf(" -----------------------------------------------------\n");

    for(iter = 0; iter < iterations; iter++)
    {
        Start(timer_inner);
        stats->processed_nodes += activeVertices;
        activeVertices = 0;


        #pragma omp parallel for private(nodes,degree,v) shared(graph,stats,bitmapNext,bitmapCurr) reduction(+ : activeVertices) schedule (dynamic,128)
        for(v = 0; v < graph->num_vertices; v++)
        {

            if(getBit(bitmapCurr, v))
            {

                degree = graph->vertices[v].out_degree;
                nodes = graph->vertices[v].outNodes;
                __u32 j;
                for(j = 0 ; j < (degree) ; j++)
                {
                    __u32   u = nodes->dest;
                    __u32   w = nodes->weight;
                    nodes = nodes->next;

                    if(numThreads == 1)
                        activeVertices += bellmanFordRelax(v, u, w, stats, bitmapNext);
                    else
                        activeVertices += bellmanFordAtomicRelax(v, u, w, stats, bitmapNext);
                }

            }
        }


        swapBitmaps(&bitmapCurr, &bitmapNext);
        clearBitmap(bitmapNext);

        Stop(timer_inner);

        printf("| %-15u | %-15u | %-15f | \n", iter, activeVertices, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }


    Stop(timer);
    stats->time_total = Seconds(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    bellmanFordPrintStatsDetails(stats);



    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);



    // bellmanFordPrintStats(stats);
    return stats;

}

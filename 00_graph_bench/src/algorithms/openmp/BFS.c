// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : BFS.c
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
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "BFS.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************

struct BFSStats *newBFSStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t vertex_id;

    struct BFSStats *stats = (struct BFSStats *) my_malloc(sizeof(struct BFSStats));

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = 0;
        if(graph->vertices->out_degree[vertex_id])
            stats->parents[vertex_id] = graph->vertices->out_degree[vertex_id] * (-1);
        else
            stats->parents[vertex_id] = -1;
    }

    return stats;

}

struct BFSStats *newBFSStatsGraphGrid(struct GraphGrid *graph)
{

    uint32_t vertex_id;

    struct BFSStats *stats = (struct BFSStats *) my_malloc(sizeof(struct BFSStats));

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;

    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = 0;
        stats->parents[vertex_id] = -1;
    }

    return stats;
}

struct BFSStats *newBFSStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    uint32_t vertex_id;

    struct BFSStats *stats = (struct BFSStats *) my_malloc(sizeof(struct BFSStats));

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = 0;
        if(graph->vertices[vertex_id].out_degree)
            stats->parents[vertex_id] = graph->vertices[vertex_id].out_degree * (-1);
        else
            stats->parents[vertex_id] = -1;
    }

    return stats;
}

struct BFSStats *newBFSStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    uint32_t vertex_id;

    struct BFSStats *stats = (struct BFSStats *) my_malloc(sizeof(struct BFSStats));

    stats->distances  = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->parents = (int *) my_malloc(graph->num_vertices * sizeof(int));
    stats->processed_nodes = 0;
    stats->iteration = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;

    // optimization for BFS implentaion instead of -1 we use -out degree to for hybrid approach counter
    #pragma omp parallel for default(none) private(vertex_id) shared(stats,graph)
    for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
    {
        stats->distances[vertex_id] = 0;
        if(graph->vertices[vertex_id].out_degree)
            stats->parents[vertex_id] = graph->vertices[vertex_id].out_degree * (-1);
        else
            stats->parents[vertex_id] = -1;
    }

    return stats;
}

void freeBFSStats(struct BFSStats *stats)
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
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************
struct BFSStats *breadthFirstSearchGraphCSR(uint32_t source, uint32_t pushpull, struct GraphCSR *graph)
{

    struct BFSStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = breadthFirstSearchPullGraphCSR(source, graph);
        break;
    case 1: // push
        stats = breadthFirstSearchPushGraphCSR(source, graph);
        break;
    case 2: // pull/push
        stats = breadthFirstSearchDirectionOptimizedGraphCSR(source, graph);
        break;
    case 3: // push-bitmap queue instead of array queue
        stats = breadthFirstSearchPushBitmapGraphCSR(source, graph);
        break;
    case 4: // pull/push-bitmap queue instead of array queue
        stats = breadthFirstSearchPushDirectionOptimizedBitmapGraphCSR(source, graph);
        break;
    default:// push
        stats = breadthFirstSearchDirectionOptimizedGraphCSR(source, graph);
        break;
    }


    return stats;

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

struct BFSStats *breadthFirstSearchPullGraphCSR(uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PULL/BU (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t nf = 0; // number of vertices in sharedFrontierQueue

    Start(timer_inner);
    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);

    printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {

        Start(timer_inner);
        nf = bottomUpStepGraphCSR(graph, sharedFrontierQueue->q_bitmap, sharedFrontierQueue->q_bitmap_next, stats);
        sharedFrontierQueue->q_bitmap_next->numSetBits = nf;
        swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
        clearBitmap(sharedFrontierQueue->q_bitmap_next);
        Stop(timer_inner);

        //stats
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += nf;
        printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);

    return stats;
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

struct BFSStats *breadthFirstSearchPushGraphCSR(uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/TD (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t P = numThreads;

    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }

    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        Start(timer_inner);
        topDownStepGraphCSR(graph, sharedFrontierQueue, localFrontierQueues, stats);
        slideWindowArrayQueue(sharedFrontierQueue);
        Stop(timer_inner);

        //stats collection
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;
        printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);


    return stats;
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

struct BFSStats *breadthFirstSearchDirectionOptimizedGraphCSR(uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/PULL(SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);
    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    uint32_t P = numThreads;
    uint32_t mu = graph->num_edges; // number of edges to check from sharedFrontierQueue
    uint32_t mf = graph->vertices->out_degree[source]; // number of edges from unexplored verticies
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    uint32_t nf_prev = 0; // number of vertices in sharedFrontierQueue
    uint32_t n = graph->num_vertices; // number of nodes
    uint32_t alpha = 15;
    uint32_t beta = 18;


    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }



    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        if(mf > (mu / alpha))
        {

            Start(timer_inner);
            arrayQueueToBitmap(sharedFrontierQueue, bitmapCurr);
            nf = sizeArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);
            printf("| E  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            do
            {
                Start(timer_inner);
                nf_prev = nf;
                nf = bottomUpStepGraphCSR(graph, bitmapCurr, bitmapNext, stats);
                swapBitmaps(&bitmapCurr, &bitmapNext);
                clearBitmap(bitmapNext);
                Stop(timer_inner);

                //stats collection
                stats->time_total +=  Seconds(timer_inner);
                stats->processed_nodes += nf;
                printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

            }
            while(( nf > nf_prev) ||  // growing;
                    ( nf > (n / beta)));

            Start(timer_inner);
            bitmapToArrayQueue(bitmapCurr, sharedFrontierQueue, localFrontierQueues);
            Stop(timer_inner);
            printf("| C  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            mf = 1;

        }
        else
        {

            Start(timer_inner);
            mu -= mf;
            mf = topDownStepGraphCSR(graph, sharedFrontierQueue, localFrontierQueues, stats);
            slideWindowArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);

            //stats collection
            stats->time_total +=  Seconds(timer_inner);
            stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;
            printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

        }



    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);
    free(timer);
    free(timer_inner);


    return stats;
}


// top-down-step(graph, sharedFrontierQueue, next, parents)
//  for v ∈ sharedFrontierQueue do
//      for u ∈ neighbors[v] do
//          if parents[u] = -1 then
//              parents[u] ← v
//              next ← next ∪ {u}
//          end if
//      end for
//  end for

uint32_t topDownStepGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{



    uint32_t v;
    uint32_t u;
    uint32_t i;
    uint32_t j;
    uint32_t edge_idx;
    uint32_t mf = 0;


    #pragma omp parallel default (none) private(u,v,j,i,edge_idx) shared(stats,localFrontierQueues,graph,sharedFrontierQueue,mf)
    {
        uint32_t t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


        #pragma omp for reduction(+:mf) schedule(auto)
        for(i = sharedFrontierQueue->head ; i < sharedFrontierQueue->tail; i++)
        {
            v = sharedFrontierQueue->queue[i];
            edge_idx = graph->vertices->edges_idx[v];

            for(j = edge_idx ; j < (edge_idx + graph->vertices->out_degree[v]) ; j++)
            {

                u = graph->sorted_edges_array->edges_array_dest[j];
                int u_parent = stats->parents[u];
                if(u_parent < 0 )
                {
                    if(__sync_bool_compare_and_swap(&stats->parents[u], u_parent, v))
                    {
                        enArrayQueue(localFrontierQueue, u);
                        mf +=  -(u_parent);
                        stats->distances[u] = stats->distances[v] + 1;
                    }
                }
            }

        }

        flushArrayQueueToShared(localFrontierQueue, sharedFrontierQueue);
    }

    return mf;
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

uint32_t bottomUpStepGraphCSR(struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats)
{


    uint32_t v;
    uint32_t u;
    uint32_t j;
    uint32_t edge_idx;
    uint32_t out_degree;
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;

    // uint32_t processed_nodes = bitmapCurr->numSetBits;
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    // stats->processed_nodes += processed_nodes;

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
        if(stats->parents[v] < 0)  // optmization
        {
            edge_idx = vertices->edges_idx[v];

            for(j = edge_idx ; j < (edge_idx + out_degree) ; j++)
            {
                u = sorted_edges_array[j];
                if(getBit(bitmapCurr, u))
                {
                    stats->parents[v] = u;
                    stats->distances[v] = stats->distances[u] + 1;
                    setBitAtomic(bitmapNext, v);
                    nf++;
                    break;
                }
            }

        }

    }
    return nf;
}


// ********************************************************************************************
// ***************      CSR DataStructure/Bitmap Frontiers                       **************
// ********************************************************************************************

// / breadth-first-search(graph, source)
//  sharedFrontierQueue ← {source}
//  next ← {}
//  parents ← [-1,-1,. . . -1]
//      while sharedFrontierQueue 6= {} do
//          top-down-step(graph, sharedFrontierQueue, next, parents)
//          sharedFrontierQueue ← next
//          next ← {}
//      end while
//  return parents

struct BFSStats *breadthFirstSearchPushBitmapGraphCSR(uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/Bitmap (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    Start(timer_inner);
    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);

    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {

        Start(timer_inner);
        topDownStepUsingBitmapsGraphCSR(graph, sharedFrontierQueue, stats);

        sharedFrontierQueue->q_bitmap_next->numSetBits = getNumOfSetBits(sharedFrontierQueue->q_bitmap_next);
        swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
        clearBitmap(sharedFrontierQueue->q_bitmap_next);
        Stop(timer_inner);


        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += sharedFrontierQueue->q_bitmap->numSetBits;
        printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->q_bitmap->numSetBits, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);

    return stats;
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

struct BFSStats *breadthFirstSearchPushDirectionOptimizedBitmapGraphCSR(uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/PULL Bitmap (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t mu = graph->num_edges; // number of edges to check from sharedFrontierQueue
    uint32_t mf = graph->vertices->out_degree[source]; // number of edges from unexplored verticies
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    uint32_t nf_prev = 0; // number of vertices in sharedFrontierQueue
    uint32_t n = graph->num_vertices; // number of nodes
    uint32_t alpha = 15;
    uint32_t beta = 18;


    Start(timer_inner);
    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {

        if(mf > (mu / alpha))
        {

            nf = sharedFrontierQueue->q_bitmap->numSetBits;
            printf("| E  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            do
            {

                Start(timer_inner);
                nf_prev = nf;
                nf = bottomUpStepGraphCSR(graph, sharedFrontierQueue->q_bitmap, sharedFrontierQueue->q_bitmap_next, stats);

                sharedFrontierQueue->q_bitmap_next->numSetBits = nf;
                swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
                clearBitmap(sharedFrontierQueue->q_bitmap_next);
                Stop(timer_inner);

                //stats
                stats->time_total +=  Seconds(timer_inner);
                stats->processed_nodes += nf;
                printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

            }
            while(( nf > nf_prev) ||  // growing;
                    ( nf > (n / beta)));

            printf("| C  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            mf = 1;

        }
        else
        {

            mu -= mf;

            Start(timer_inner);
            mf = topDownStepUsingBitmapsGraphCSR(graph, sharedFrontierQueue, stats);

            sharedFrontierQueue->q_bitmap_next->numSetBits = getNumOfSetBits(sharedFrontierQueue->q_bitmap_next);
            swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
            clearBitmap(sharedFrontierQueue->q_bitmap_next);
            Stop(timer_inner);


            stats->time_total +=  Seconds(timer_inner);
            stats->processed_nodes += sharedFrontierQueue->q_bitmap->numSetBits;
            printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->q_bitmap->numSetBits, Seconds(timer_inner));

        }



    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);

    return stats;
}


uint32_t topDownStepUsingBitmapsGraphCSR(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct BFSStats *stats)
{



    uint32_t v;
    uint32_t u;
    uint32_t i;
    uint32_t j;
    uint32_t edge_idx;
    uint32_t mf = 0;

    #pragma omp parallel default (none) private(u,v,j,i,edge_idx) shared(stats,graph,sharedFrontierQueue,mf)
    {


        #pragma omp for reduction(+:mf)
        for(i = 0 ; i < (sharedFrontierQueue->q_bitmap->size); i++)
        {
            if(getBit(sharedFrontierQueue->q_bitmap, i))
            {
                // processed_nodes++;
                v = i;
                edge_idx = graph->vertices->edges_idx[v];

                for(j = edge_idx ; j < (edge_idx + graph->vertices->out_degree[v]) ; j++)
                {


                    u = graph->sorted_edges_array->edges_array_dest[j];
                    int u_parent = stats->parents[u];

                    if(u_parent < 0 )
                    {
                        if(__sync_bool_compare_and_swap(&stats->parents[u], u_parent, v))
                        {
                            mf +=  -(u_parent);
                            stats->distances[u] = stats->distances[v] + 1;
                            setBitAtomic(sharedFrontierQueue->q_bitmap_next, u);


                        }

                    }
                }

            }

        }

    }



    return mf;
}



// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphGrid(uint32_t source, uint32_t pushpull, struct GraphGrid *graph)
{

    struct BFSStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = breadthFirstSearchRowGraphGrid(source, graph);
        break;
    case 1: // push
        stats = breadthFirstSearchRowGraphGridBitmap(source, graph);
        break;
    case 2: // pull
        stats = breadthFirstSearchColumnGraphGrid(source, graph);
        break;
    case 3: // push
        stats = breadthFirstSearchColumnGraphGridBitmap(source, graph);
        break;
    default:// push
        stats = breadthFirstSearchRowGraphGrid(source, graph);
        break;
    }

    return stats;

}

// function STREAMVERTICES(Fv,F)
//  Sum = 0
//      for each vertex do
//          if F(vertex) then
//              Sum += Fv(edge)
//          end if
//      end for
//  return Sum
// end function

// function STREAMEDGES(Fe,F)
//  Sum = 0
//      for each active block do >> block with active edges
//          for each edge ∈ block do
//              if F(edge.source) then
//                  Sum += Fe(edge)
//              end if
//          end for
//      end for
//  return Sum
// end function
//we assume that the edges are not sorted in each partition

struct BFSStats *breadthFirstSearchRowGraphGrid(uint32_t source, struct GraphGrid *graph)
{
    struct BFSStats *stats = newBFSStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS-Row (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_iteration = (struct Timer *) malloc(sizeof(struct Timer));
    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);


    uint32_t P = numThreads;



    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    #pragma omp parallel for
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);
    }


    graphGridReset(graph);

    uint32_t processed_nodes = 0;

    Start(timer_iteration);
    enArrayQueue(sharedFrontierQueue, source);
    arrayQueueGenerateBitmap(sharedFrontierQueue);
    stats->parents[source] = source;
    // graphGridSetActivePartitions(graph->grid, source);
    graphGridSetActivePartitionsMap(graph->grid, source);
    Stop(timer_iteration);


    printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, ++processed_nodes, Seconds(timer_iteration));

    stats->time_total += Seconds(timer_iteration);
    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        Start(timer_iteration);
        breadthFirstSearchStreamEdgesRowGraphGrid(graph, sharedFrontierQueue, localFrontierQueues, stats);
        Stop(timer_iteration);


        processed_nodes = sharedFrontierQueue->tail_next - sharedFrontierQueue->tail;
        slideWindowArrayQueue(sharedFrontierQueue);
        arrayQueueGenerateBitmap(sharedFrontierQueue);
        breadthFirstSearchSetActivePartitions(graph, sharedFrontierQueue);

        stats->time_total += Seconds(timer_iteration);
        printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));
    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", sharedFrontierQueue->tail_next, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "**", sharedFrontierQueue->tail_next, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    freeArrayQueue(sharedFrontierQueue);
    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }

    //   #pragma omp parallel for
    //   for(i=0 ; i < P*P ; i++){
    //  freeArrayQueue(localFrontierQueuesL2[i]);
    // }

    // free(localFrontierQueuesL2);
    free(localFrontierQueues);
    free(timer_iteration);
    free(timer);

    return stats;
}

struct BFSStats *breadthFirstSearchColumnGraphGrid(uint32_t source, struct GraphGrid *graph)
{
    struct BFSStats *stats = newBFSStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS-Column (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_iteration = (struct Timer *) malloc(sizeof(struct Timer));
    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);


    uint32_t P = numThreads;



    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    #pragma omp parallel for
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);
    }


    graphGridReset(graph);

    uint32_t processed_nodes = 0;

    Start(timer_iteration);
    enArrayQueue(sharedFrontierQueue, source);
    arrayQueueGenerateBitmap(sharedFrontierQueue);
    stats->parents[source] = source;
    // graphGridSetActivePartitions(graph->grid, source);
    graphGridSetActivePartitionsMap(graph->grid, source);
    Stop(timer_iteration);


    printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, ++processed_nodes, Seconds(timer_iteration));

    stats->time_total += Seconds(timer_iteration);
    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        Start(timer_iteration);
        breadthFirstSearchStreamEdgesColumnGraphGrid(graph, sharedFrontierQueue, localFrontierQueues, stats);
        Stop(timer_iteration);


        processed_nodes = sharedFrontierQueue->tail_next - sharedFrontierQueue->tail;
        slideWindowArrayQueue(sharedFrontierQueue);
        arrayQueueGenerateBitmap(sharedFrontierQueue);
        breadthFirstSearchSetActivePartitions(graph, sharedFrontierQueue);

        stats->time_total += Seconds(timer_iteration);
        printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));
    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", sharedFrontierQueue->tail_next, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "**", sharedFrontierQueue->tail_next, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    freeArrayQueue(sharedFrontierQueue);
    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }

    //   #pragma omp parallel for
    //   for(i=0 ; i < P*P ; i++){
    //  freeArrayQueue(localFrontierQueuesL2[i]);
    // }

    // free(localFrontierQueuesL2);
    free(localFrontierQueues);
    free(timer_iteration);
    free(timer);

    return stats;
}

// function STREAMEDGES(Fe,F)
//  Sum = 0
//      for each active block do >> block with active edges
//          for each edge ∈ block do
//              if F(edge.source) then
//                  Sum += Fe(edge)
//              end if
//          end for
//      end for
//  return Sum
// end function
//we assume that the edges are not sorted in each partition
void breadthFirstSearchStreamEdgesRowGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{
    // struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
    uint32_t totalPartitions = 0;
    totalPartitions = graph->grid->num_partitions; // PxP




    uint32_t i;

    for (i = 0; i < totalPartitions; ++i)
    {
        uint32_t j;
        #pragma omp parallel for default(none) shared(i,stats,totalPartitions,localFrontierQueues ,sharedFrontierQueue, graph)
        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t t_id = omp_get_thread_num();
            // uint32_t A = 0;
            struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


            if(getBit(graph->grid->activePartitionsMap, (i * totalPartitions) + j))
            {
                // #pragma  omp task untied
                // {

                breadthFirstSearchPartitionGraphGrid(graph, &(graph->grid->partitions[(i * totalPartitions) + j]), sharedFrontierQueue, localFrontierQueue, stats);
                flushArrayQueueToShared(localFrontierQueue, sharedFrontierQueue);
                // }

            }
        }
    }
    // flushArrayQueueToShared(localFrontierQueue,sharedFrontierQueue);
    // }
}


void breadthFirstSearchStreamEdgesColumnGraphGrid(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{
    // struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
    uint32_t totalPartitions = 0;
    totalPartitions = graph->grid->num_partitions; // PxP


    #pragma omp parallel default(none) shared(stats,totalPartitions,localFrontierQueues ,sharedFrontierQueue, graph)
    // #pragma  omp single nowait
    {

        uint32_t t_id = omp_get_thread_num();
        // uint32_t A = 0;
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


        uint32_t j;
        #pragma omp for
        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t i;
            for (i = 0; i < totalPartitions; ++i)
            {


                if(getBit(graph->grid->activePartitionsMap, (i * totalPartitions) + j))
                {
                    // #pragma  omp task untied
                    // {

                    breadthFirstSearchPartitionGraphGrid(graph, &(graph->grid->partitions[(i * totalPartitions) + j]), sharedFrontierQueue, localFrontierQueue, stats);
                    flushArrayQueueToShared(localFrontierQueue, sharedFrontierQueue);
                    // }

                }
            }
        }

    }

    // flushArrayQueueToShared(localFrontierQueue,sharedFrontierQueue);
    // }
}

void breadthFirstSearchPartitionGraphGrid(struct GraphGrid *graph, struct Partition *partition, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue *localFrontierQueue, struct BFSStats *stats)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;


    // #pragma omp parallel default(none) private(i,src,dest) shared(localFrontierQueuesL2,graph,partition,sharedFrontierQueue,localFrontierQueue)
    //    {

    //        uint32_t t_id = omp_get_thread_num();
    //        struct ArrayQueue* localFrontierQueueL2 = localFrontierQueuesL2[t_id];

    //  #pragma omp for schedule(dynamic, 1024)
    for (i = 0; i < partition->num_edges; ++i)
    {

        src  = partition->edgeList->edges_array_src[i];
        dest = partition->edgeList->edges_array_dest[i];
        int v_dest = stats->parents[dest];
        if(isEnArrayQueued(sharedFrontierQueue, src) && (v_dest < 0))
        {
            // if(__sync_bool_compare_and_swap(&stats->parents[dest], v_dest, src))
            // {
            stats->parents[dest] = src;
            stats->distances[dest] = stats->distances[src] + 1;
            enArrayQueue(localFrontierQueue, dest);
            // }
        }
    }

    //      flushArrayQueueToShared(localFrontierQueueL2,localFrontierQueue);
    //      // slideWindowArrayQueue(localFrontierQueue);
    //      localFrontierQueue->tail = localFrontierQueue->tail_next; // to apply to condition to the next flush
    // }


}

void breadthFirstSearchSetActivePartitions(struct GraphGrid *graph, struct ArrayQueue *sharedFrontierQueue)
{

    uint32_t i;
    uint32_t v;

    // graphGridResetActivePartitions(graph->grid);
    graphGridResetActivePartitionsMap(graph->grid);

    #pragma omp parallel for default(none) shared(graph,sharedFrontierQueue) private(i,v) schedule(dynamic,1024)
    for(i = sharedFrontierQueue->head ; i < sharedFrontierQueue->tail; i++)
    {
        v = sharedFrontierQueue->queue[i];
        // graphGridSetActivePartitions(graph->grid, v);
        // if(getBit(graph->grid->activePartitionsMap,i))
        graphGridSetActivePartitionsMap(graph->grid, v);
    }
}


// ********************************************************************************************
// ***************                  GRID DataStructure/Bitmap Frontiers          **************
// ********************************************************************************************

// function STREAMVERTICES(Fv,F)
//  Sum = 0
//      for each vertex do
//          if F(vertex) then
//              Sum += Fv(edge)
//          end if
//      end for
//  return Sum
// end function

// function STREAMEDGES(Fe,F)
//  Sum = 0
//      for each active block do >> block with active edges
//          for each edge ∈ block do
//              if F(edge.source) then
//                  Sum += Fe(edge)
//              end if
//          end for
//      end for
//  return Sum
// end function
//we assume that the edges are not sorted in each partition

struct BFSStats *breadthFirstSearchRowGraphGridBitmap(uint32_t source, struct GraphGrid *graph)
{

    struct BFSStats *stats = newBFSStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS-Row Bitmap (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_iteration = (struct Timer *) malloc(sizeof(struct Timer));
    struct Bitmap *FrontierBitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *FrontierBitmapNext = newBitmap(graph->num_vertices);



    graphGridReset(graph);
    uint32_t processed_nodes = 0;
    uint32_t total_processed_nodes = 0;

    Start(timer_iteration);
    setBit(FrontierBitmapNext, source);
    stats->parents[source] = source;
    processed_nodes = getNumOfSetBits(FrontierBitmapNext);
    swapBitmaps (&FrontierBitmapCurr, &FrontierBitmapNext);
    clearBitmap(FrontierBitmapNext);
    // printf("%u %u\n",getNumOfSetBits(FrontierBitmapCurr),getNumOfSetBits(FrontierBitmapNext) );
    breadthFirstSearchSetActivePartitionsBitmap(graph, FrontierBitmapCurr);


    Stop(timer_iteration);


    total_processed_nodes += processed_nodes;
    printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));

    stats->time_total += Seconds(timer_iteration);
    Start(timer);

    while(processed_nodes)  // start while
    {

        Start(timer_iteration);
        breadthFirstSearchStreamEdgesRowGraphGridBitmap(graph, FrontierBitmapCurr, FrontierBitmapNext, stats);
        Stop(timer_iteration);

        processed_nodes = getNumOfSetBits(FrontierBitmapNext);
        swapBitmaps (&FrontierBitmapCurr, &FrontierBitmapNext);
        clearBitmap(FrontierBitmapNext);
        breadthFirstSearchSetActivePartitionsBitmap(graph, FrontierBitmapCurr);
        total_processed_nodes += processed_nodes;
        stats->time_total += Seconds(timer_iteration);
        printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));
    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", total_processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "**", total_processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeBitmap(FrontierBitmapCurr);
    freeBitmap(FrontierBitmapNext);
    free(timer_iteration);
    free(timer);

    return stats;
}

struct BFSStats *breadthFirstSearchColumnGraphGridBitmap(uint32_t source, struct GraphGrid *graph)
{

    struct BFSStats *stats = newBFSStatsGraphGrid(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS-Column Bitmap (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_iteration = (struct Timer *) malloc(sizeof(struct Timer));
    struct Bitmap *FrontierBitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *FrontierBitmapNext = newBitmap(graph->num_vertices);



    graphGridReset(graph);
    uint32_t processed_nodes = 0;
    uint32_t total_processed_nodes = 0;

    Start(timer_iteration);
    setBit(FrontierBitmapNext, source);
    stats->parents[source] = source;
    processed_nodes = getNumOfSetBits(FrontierBitmapNext);
    swapBitmaps (&FrontierBitmapCurr, &FrontierBitmapNext);
    clearBitmap(FrontierBitmapNext);
    // printf("%u %u\n",getNumOfSetBits(FrontierBitmapCurr),getNumOfSetBits(FrontierBitmapNext) );
    breadthFirstSearchSetActivePartitionsBitmap(graph, FrontierBitmapCurr);


    Stop(timer_iteration);


    total_processed_nodes += processed_nodes;
    printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));

    stats->time_total += Seconds(timer_iteration);
    Start(timer);

    while(processed_nodes)  // start while
    {

        Start(timer_iteration);
        breadthFirstSearchStreamEdgesColumnGraphGridBitmap(graph, FrontierBitmapCurr, FrontierBitmapNext, stats);
        Stop(timer_iteration);

        processed_nodes = getNumOfSetBits(FrontierBitmapNext);
        swapBitmaps (&FrontierBitmapCurr, &FrontierBitmapNext);
        clearBitmap(FrontierBitmapNext);
        breadthFirstSearchSetActivePartitionsBitmap(graph, FrontierBitmapCurr);
        total_processed_nodes += processed_nodes;
        stats->time_total += Seconds(timer_iteration);
        printf("| %-15u | %-15u | %-15f | \n", stats->iteration++, processed_nodes, Seconds(timer_iteration));
    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", total_processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "**", total_processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeBitmap(FrontierBitmapCurr);
    freeBitmap(FrontierBitmapNext);
    free(timer_iteration);
    free(timer);

    return stats;
}

// function STREAMEDGES(Fe,F)
//  Sum = 0
//      for each active block do >> block with active edges
//          for each edge ∈ block do
//              if F(edge.source) then
//                  Sum += Fe(edge)
//              end if
//          end for
//      end for
//  return Sum
// end function
//we assume that the edges are not sorted in each partition
void breadthFirstSearchStreamEdgesRowGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats)
{
    // struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
    uint32_t totalPartitions = 0;
    totalPartitions = graph->grid->num_partitions; // PxP

    uint32_t i;
    for (i = 0; i < totalPartitions; ++i)
    {
        uint32_t j;
        #pragma omp parallel for default(none) shared(i,stats,totalPartitions,FrontierBitmapCurr ,FrontierBitmapNext, graph)
        for (j = 0; j < totalPartitions; ++j)
        {

            if(getBit(graph->grid->activePartitionsMap, (i * totalPartitions) + j) && graph->grid->partitions[(i * totalPartitions) + j].num_edges)
            {
                breadthFirstSearchPartitionGraphGridBitmap(graph, &(graph->grid->partitions[(i * totalPartitions) + j]), FrontierBitmapCurr, FrontierBitmapNext, stats);
            }
        }
    }

}

void breadthFirstSearchStreamEdgesColumnGraphGridBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats)
{
    // struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
    uint32_t totalPartitions = 0;
    totalPartitions = graph->grid->num_partitions; // PxP


    #pragma omp parallel default(none) shared(stats,totalPartitions,FrontierBitmapCurr ,FrontierBitmapNext, graph)
    // #pragma  omp single nowait
    {
        uint32_t j;

        // #pragma omp for schedule(dynamic, 256)
        #pragma omp for
        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t i;
            for (i = 0; i < totalPartitions; ++i)
            {
                if(getBit(graph->grid->activePartitionsMap, (i * totalPartitions) + j) && graph->grid->partitions[(i * totalPartitions) + j].num_edges)
                {
                    breadthFirstSearchPartitionGraphGridBitmap(graph, &(graph->grid->partitions[(i * totalPartitions) + j]), FrontierBitmapCurr, FrontierBitmapNext, stats);
                }
            }
        }
    }
}


void breadthFirstSearchPartitionGraphGridBitmap(struct GraphGrid *graph, struct Partition *partition, struct Bitmap *FrontierBitmapCurr, struct Bitmap *FrontierBitmapNext, struct BFSStats *stats)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;


    for (i = 0; i < partition->num_edges; ++i)
    {

        src  = partition->edgeList->edges_array_src[i];
        dest = partition->edgeList->edges_array_dest[i];
        int v_dest = stats->parents[dest];
        if((v_dest < 0))
        {
            if(getBit(FrontierBitmapCurr, src))
            {
                // if(__sync_bool_compare_and_swap(&stats->parents[dest], v_dest, src))
                // {
                stats->parents[dest] = src;
                stats->distances[dest] = stats->distances[src] + 1;
                setBitAtomic(FrontierBitmapNext, dest);
                // }
            }
        }
    }


}

void breadthFirstSearchSetActivePartitionsBitmap(struct GraphGrid *graph, struct Bitmap *FrontierBitmap)
{

    uint32_t i;

    graphGridResetActivePartitionsMap(graph->grid);

    #pragma omp parallel for default(none) shared(graph,FrontierBitmap) private(i) schedule(dynamic,1024)
    for(i = 0 ; i < FrontierBitmap->size; i++)
    {
        if(getBit(FrontierBitmap, i))
            graphGridSetActivePartitionsMap(graph->grid, i);
    }
}


// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphAdjArrayList(uint32_t source, uint32_t pushpull, struct GraphAdjArrayList *graph)
{

    struct BFSStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = breadthFirstSearchPullGraphAdjArrayList(source, graph);
        break;
    case 1: // push
        stats = breadthFirstSearchPushGraphAdjArrayList(source, graph);
        break;
    case 2: // pull/push
        stats = breadthFirstSearchDirectionOptimizedGraphAdjArrayList(source, graph);
        break;
    default:// push
        stats = breadthFirstSearchDirectionOptimizedGraphAdjArrayList(source, graph);
        break;
    }

    return stats;

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


struct BFSStats *breadthFirstSearchPullGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjArrayList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PULL/BU (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t nf = 0; // number of vertices in sharedFrontierQueue



    Start(timer_inner);
    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);

    printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {

        Start(timer_inner);
        nf = bottomUpStepGraphAdjArrayList(graph, sharedFrontierQueue->q_bitmap, sharedFrontierQueue->q_bitmap_next, stats);
        sharedFrontierQueue->q_bitmap_next->numSetBits = nf;
        swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
        clearBitmap(sharedFrontierQueue->q_bitmap_next);
        Stop(timer_inner);

        //stats
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += nf;
        printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);

    return stats;
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


struct BFSStats *breadthFirstSearchPushGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjArrayList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/TD (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t P = numThreads;

    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }


    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        Start(timer_inner);
        topDownStepGraphAdjArrayList(graph, sharedFrontierQueue, localFrontierQueues, stats);
        slideWindowArrayQueue(sharedFrontierQueue);
        Stop(timer_inner);

        //stats collection
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;
        printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);


    return stats;
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


struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjArrayList(uint32_t source, struct GraphAdjArrayList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjArrayList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);
    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    uint32_t P = numThreads;
    uint32_t mu = graph->num_edges; // number of edges to check from sharedFrontierQueue
    uint32_t mf = graph->vertices[source].out_degree; // number of edges from unexplored verticies
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    uint32_t nf_prev = 0; // number of vertices in sharedFrontierQueue
    uint32_t n = graph->num_vertices; // number of nodes
    uint32_t alpha = 15;
    uint32_t beta = 18;


    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }

    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        if(mf > (mu / alpha))
        {

            Start(timer_inner);
            arrayQueueToBitmap(sharedFrontierQueue, bitmapCurr);
            nf = sizeArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);
            printf("| E  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            do
            {
                Start(timer_inner);
                nf_prev = nf;
                nf = bottomUpStepGraphAdjArrayList(graph, bitmapCurr, bitmapNext, stats);
                swapBitmaps(&bitmapCurr, &bitmapNext);
                clearBitmap(bitmapNext);
                Stop(timer_inner);

                //stats collection
                stats->time_total +=  Seconds(timer_inner);
                stats->processed_nodes += nf;
                printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

            }
            while(( nf > nf_prev) ||  // growing;
                    ( nf > (n / beta)));

            Start(timer_inner);
            bitmapToArrayQueue(bitmapCurr, sharedFrontierQueue, localFrontierQueues);
            Stop(timer_inner);
            printf("| C  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            mf = 1;

        }
        else
        {

            Start(timer_inner);
            mu -= mf;
            mf = topDownStepGraphAdjArrayList(graph, sharedFrontierQueue, localFrontierQueues, stats);
            slideWindowArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);
            //stats collection
            stats->time_total +=  Seconds(timer_inner);
            stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;;
            printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

        }



    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);
    free(timer);
    free(timer_inner);

    return stats;
}


// top-down-step(graph, sharedFrontierQueue, next, parents)
//  for v ∈ sharedFrontierQueue do
//      for u ∈ neighbors[v] do
//          if parents[u] = -1 then
//              parents[u] ← v
//              next ← next ∪ {u}
//          end if
//      end for
//  end for

uint32_t topDownStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{



    uint32_t v;
    uint32_t u;
    uint32_t i;
    uint32_t j;
    uint32_t mf = 0;

    uint32_t out_degree;
    struct EdgeList *outNodes;

    #pragma omp parallel default (none) private(out_degree,outNodes,u,v,j,i) shared(stats,localFrontierQueues,graph,sharedFrontierQueue,mf)
    {
        uint32_t t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


        #pragma omp for reduction(+:mf) schedule(auto)
        for(i = sharedFrontierQueue->head ; i < sharedFrontierQueue->tail; i++)
        {
            v = sharedFrontierQueue->queue[i];
            // v = deArrayQueue(sharedFrontierQueue);
            outNodes = graph->vertices[v].outNodes;
            out_degree = graph->vertices[v].out_degree;

            for(j = 0 ; j < out_degree ; j++)
            {

                u = outNodes->edges_array_dest[j];
                int u_parent = stats->parents[u];
                if(u_parent < 0 )
                {
                    if(__sync_bool_compare_and_swap(&stats->parents[u], u_parent, v))
                    {
                        enArrayQueue(localFrontierQueue, u);
                        stats->distances[u] = stats->distances[v] + 1;
                        mf +=  -(u_parent);
                    }
                }
            }

        }

        flushArrayQueueToShared(localFrontierQueue, sharedFrontierQueue);
    }

    return mf;
}

// bottom-up-step(graph, sharedFrontierQueue, next, parents)
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

uint32_t bottomUpStepGraphAdjArrayList(struct GraphAdjArrayList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats)
{


    uint32_t v;
    uint32_t u;
    uint32_t j;

    // uint32_t processed_nodes = bitmapCurr->numSetBits;
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    // stats->processed_nodes += processed_nodes;


    uint32_t degree;
    struct EdgeList *Nodes;


    #pragma omp parallel for default(none) private(Nodes,j,u,v,degree) shared(stats,bitmapCurr,bitmapNext,graph) reduction(+:nf) schedule(dynamic, 1024)
    for(v = 0 ; v < graph->num_vertices ; v++)
    {
        if(stats->parents[v] < 0)  // optmization
        {

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[v].inNodes;
            degree = graph->vertices[v].in_degree;
#else
            Nodes = graph->vertices[v].outNodes;
            degree = graph->vertices[v].out_degree;
#endif

            for(j = 0 ; j < (degree) ; j++)
            {
                u = Nodes->edges_array_dest[j];
                if(getBit(bitmapCurr, u))
                {
                    stats->parents[v] = u;
                    setBitAtomic(bitmapNext, v);
                    stats->distances[v] = stats->distances[u] + 1;
                    nf++;
                    break;
                }
            }

        }

    }

    return nf;
}


// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct BFSStats *breadthFirstSearchGraphAdjLinkedList(uint32_t source, uint32_t pushpull, struct GraphAdjLinkedList *graph)
{

    struct BFSStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = breadthFirstSearchPullGraphAdjLinkedList(source, graph);
        break;
    case 1: // push
        stats = breadthFirstSearchPushGraphAdjLinkedList(source, graph);
        break;
    case 2: // pull/push
        stats = breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(source, graph);
        break;
    default:// push
        stats = breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(source, graph);
        break;
    }

    return stats;

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


struct BFSStats *breadthFirstSearchPullGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjLinkedList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PULL/BU (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }


    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t nf = 0; // number of vertices in sharedFrontierQueue



    Start(timer_inner);
    setBit(sharedFrontierQueue->q_bitmap_next, source);
    sharedFrontierQueue->q_bitmap_next->numSetBits = 1;
    stats->parents[source] = source;

    swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
    clearBitmap(sharedFrontierQueue->q_bitmap_next);
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);

    printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while (sharedFrontierQueue->q_bitmap->numSetBits)
    {

        Start(timer_inner);
        nf = bottomUpStepGraphAdjLinkedList(graph, sharedFrontierQueue->q_bitmap, sharedFrontierQueue->q_bitmap_next, stats);
        sharedFrontierQueue->q_bitmap_next->numSetBits = nf;
        swapBitmaps(&sharedFrontierQueue->q_bitmap, &sharedFrontierQueue->q_bitmap_next);
        clearBitmap(sharedFrontierQueue->q_bitmap_next);
        Stop(timer_inner);

        //stats
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += nf;
        printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);

    return stats;
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


struct BFSStats *breadthFirstSearchPushGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjLinkedList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PUSH/TD (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);

    uint32_t P = numThreads;
    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }

    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;


    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        Start(timer_inner);
        topDownStepGraphAdjLinkedList(graph, sharedFrontierQueue, localFrontierQueues, stats);
        slideWindowArrayQueue(sharedFrontierQueue);
        Stop(timer_inner);

        //stats collection
        stats->time_total +=  Seconds(timer_inner);
        stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;
        printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    free(timer);
    free(timer_inner);


    return stats;
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


struct BFSStats *breadthFirstSearchDirectionOptimizedGraphAdjLinkedList(uint32_t source, struct GraphAdjLinkedList *graph)
{

    struct BFSStats *stats = newBFSStatsGraphAdjLinkedList(graph);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS PULL/PUSH (SOURCE NODE)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", source);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iteration", "Nodes", "Time (Seconds)");
    printf(" -----------------------------------------------------\n");

    if(source > graph->num_vertices)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "ERROR!! CHECK SOURCE RANGE");
        printf(" -----------------------------------------------------\n");
        return stats;
    }

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct ArrayQueue *sharedFrontierQueue = newArrayQueue(graph->num_vertices);
    struct Bitmap *bitmapCurr = newBitmap(graph->num_vertices);
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    uint32_t P = numThreads;
    uint32_t mu = graph->num_edges; // number of edges to check from sharedFrontierQueue
    uint32_t mf = graph->vertices[source].out_degree; // number of edges from unexplored verticies
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    uint32_t nf_prev = 0; // number of vertices in sharedFrontierQueue
    uint32_t n = graph->num_vertices; // number of nodes
    uint32_t alpha = 15;
    uint32_t beta = 18;


    struct ArrayQueue **localFrontierQueues = (struct ArrayQueue **) my_malloc( P * sizeof(struct ArrayQueue *));


    uint32_t i;
    for(i = 0 ; i < P ; i++)
    {
        localFrontierQueues[i] = newArrayQueue(graph->num_vertices);

    }



    Start(timer_inner);
    enArrayQueue(sharedFrontierQueue, source);
    // setBit(sharedFrontierQueue->q_bitmap,source);
    stats->parents[source] = source;
    Stop(timer_inner);
    stats->time_total +=  Seconds(timer_inner);
    // graph->vertices[source].visited = 1;

    printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, ++stats->processed_nodes, Seconds(timer_inner));

    Start(timer);
    while(!isEmptyArrayQueue(sharedFrontierQueue))  // start while
    {

        if(mf > (mu / alpha))
        {

            Start(timer_inner);
            arrayQueueToBitmap(sharedFrontierQueue, bitmapCurr);
            nf = sizeArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);
            printf("| E  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            do
            {
                Start(timer_inner);
                nf_prev = nf;
                nf = bottomUpStepGraphAdjLinkedList(graph, bitmapCurr, bitmapNext, stats);
                swapBitmaps(&bitmapCurr, &bitmapNext);
                clearBitmap(bitmapNext);
                Stop(timer_inner);

                //stats collection
                stats->time_total +=  Seconds(timer_inner);
                stats->processed_nodes += nf;
                printf("| BU %-12u | %-15u | %-15f | \n", stats->iteration++, nf, Seconds(timer_inner));

            }
            while(( nf > nf_prev) ||  // growing;
                    ( nf > (n / beta)));

            Start(timer_inner);
            bitmapToArrayQueue(bitmapCurr, sharedFrontierQueue, localFrontierQueues);
            Stop(timer_inner);
            printf("| C  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            mf = 1;

        }
        else
        {

            Start(timer_inner);
            mu -= mf;
            mf = topDownStepGraphAdjLinkedList(graph, sharedFrontierQueue, localFrontierQueues, stats);
            slideWindowArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);
            //stats collection
            stats->time_total +=  Seconds(timer_inner);
            stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;;
            printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

        }



    } // end while
    Stop(timer);
    stats->time_total =  Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    for(i = 0 ; i < P ; i++)
    {
        freeArrayQueue(localFrontierQueues[i]);
    }
    free(localFrontierQueues);
    freeArrayQueue(sharedFrontierQueue);
    freeBitmap(bitmapNext);
    freeBitmap(bitmapCurr);
    free(timer);
    free(timer_inner);

    return stats;
}


// top-down-step(graph, sharedFrontierQueue, next, parents)
//  for v ∈ sharedFrontierQueue do
//      for u ∈ neighbors[v] do
//          if parents[u] = -1 then
//              parents[u] ← v
//              next ← next ∪ {u}
//          end if
//      end for
//  end for

uint32_t topDownStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{



    uint32_t v;
    uint32_t u;
    uint32_t i;
    uint32_t j;
    uint32_t mf = 0;

    uint32_t out_degree;
    struct AdjLinkedListNode *outNodes;

    #pragma omp parallel default (none) private(out_degree,outNodes,u,v,j,i) shared(stats,localFrontierQueues,graph,sharedFrontierQueue,mf)
    {
        uint32_t t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


        #pragma omp for reduction(+:mf) schedule(auto)
        for(i = sharedFrontierQueue->head ; i < sharedFrontierQueue->tail; i++)
        {
            v = sharedFrontierQueue->queue[i];
            // v = deArrayQueue(sharedFrontierQueue);
            outNodes = graph->vertices[v].outNodes;
            out_degree = graph->vertices[v].out_degree;

            for(j = 0 ; j < out_degree ; j++)
            {

                u = outNodes->dest;
                outNodes = outNodes->next; // travers pointer

                int u_parent = stats->parents[u];
                if(u_parent < 0 )
                {
                    if(__sync_bool_compare_and_swap(&stats->parents[u], u_parent, v))
                    {
                        enArrayQueue(localFrontierQueue, u);
                        stats->distances[u] = stats->distances[v] + 1;
                        mf +=  -(u_parent);
                    }
                }
            }

        }

        flushArrayQueueToShared(localFrontierQueue, sharedFrontierQueue);
    }

    return mf;
}

// bottom-up-step(graph, sharedFrontierQueue, next, parents)
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

uint32_t bottomUpStepGraphAdjLinkedList(struct GraphAdjLinkedList *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats)
{


    uint32_t v;
    uint32_t u;
    uint32_t j;

    // uint32_t processed_nodes = bitmapCurr->numSetBits;
    uint32_t nf = 0; // number of vertices in sharedFrontierQueue
    // stats->processed_nodes += processed_nodes;


    uint32_t degree;
    struct AdjLinkedListNode *Nodes;


    #pragma omp parallel for default(none) private(Nodes,j,u,v,degree) shared(stats,bitmapCurr,bitmapNext,graph) reduction(+:nf) schedule(dynamic, 1024)
    for(v = 0 ; v < graph->num_vertices ; v++)
    {
        if(stats->parents[v] < 0)  // optmization
        {

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[v].inNodes;
            degree = graph->vertices[v].in_degree;
#else
            Nodes = graph->vertices[v].outNodes;
            degree = graph->vertices[v].out_degree;
#endif

            for(j = 0 ; j < (degree) ; j++)
            {
                u = Nodes->dest;
                Nodes = Nodes->next;
                if(getBit(bitmapCurr, u))
                {
                    stats->parents[v] = u;
                    setBitAtomic(bitmapNext, v);
                    stats->distances[v] = stats->distances[u] + 1;
                    nf++;
                    break;
                }
            }

        }

    }

    return nf;
}

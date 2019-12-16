// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : epochReorder.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:35:52
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <limits.h>

#include "graphCSR.h"
#include "mt19937.h"
#include "reorder.h"
#include "arrayQueue.h"
#include "myMalloc.h"
#include "epochReorder.h"
#include "bitmap.h"
#include "timer.h"

#include "pageRank.h"
#include "BFS.h"


uint32_t epochAtomicMin(uint32_t *dist, uint32_t newValue)
{

    uint32_t oldValue;
    uint32_t flag = 0;

    do
    {

        oldValue = *dist;
        if(oldValue > newValue)
        {
            if(__sync_bool_compare_and_swap(dist, oldValue, 0))
            {
                flag = 1;
                // printf("| %-15s | %-30u | \n","Hard hardThreshold", newValue);
                // printf("| %-15s | %-30u | \n","Hard counter", oldValue);
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



struct EpochReorder *newEpochReoder( uint32_t softThreshold, uint32_t hardThreshold, uint32_t numCounters, uint32_t numVertices)
{

    uint32_t v = 0;
    uint32_t n = 0;

    struct EpochReorder *epochReorder = (struct EpochReorder *) my_malloc(sizeof(struct EpochReorder));


    epochReorder->rrIndex = 0;
    epochReorder->softcounter = 0;
    epochReorder->hardcounter = 0;
    epochReorder->softThreshold = softThreshold;
    epochReorder->hardThreshold = hardThreshold;
    epochReorder->numCounters = numCounters;
    epochReorder->numVertices = numVertices;

    epochReorder->recencyBits = newBitmap(numVertices);


    epochReorder->frequency = (uint32_t **) my_malloc(sizeof(uint32_t *) * numCounters);
    epochReorder->reuse = (uint32_t **) my_malloc(sizeof(uint32_t *) * numCounters);

    epochReorder->base_reuse = (uint32_t *) my_malloc(sizeof(uint32_t) * numVertices);

    for(v = 0; v < numCounters; v++)
    {
        epochReorder->frequency[v] = (uint32_t *) my_malloc(sizeof(uint32_t) * numVertices);
        epochReorder->reuse[v] = (uint32_t *) my_malloc(sizeof(uint32_t) * numVertices);
    }



    #pragma omp parallel for
    for(v = 0; v < numCounters; v++)
    {
        for(n = 0; n < numVertices; n++)
        {
            epochReorder->frequency[v][n] = 0;
            epochReorder->reuse[v][n] = 0;
        }
    }

    #pragma omp parallel for
    for(n = 0; n < numVertices; n++)
    {
        epochReorder->base_reuse[n] = UINT_MAX;
    }

    return epochReorder;

}

uint32_t *epochReorderPageRank(struct GraphCSR *graph)
{


    uint32_t numCounters = 20;
    uint32_t hardThreshold = 4096;
    uint32_t softThreshold = 8192;
    double epsilon = 1e-6;
    uint32_t iterations = 2;
    uint32_t *labels;

    struct EpochReorder *epochReorder = newEpochReoder(softThreshold, hardThreshold, numCounters, graph->num_vertices);

    epochReorderPageRankPullGraphCSR(epochReorder, epsilon, iterations, graph);

    labels = epochReorderCreateLabels(epochReorder);

    freeEpochReorder(epochReorder);

    return labels;

}



float *epochReorderPageRankPullGraphCSR(struct EpochReorder *epochReorder, double epsilon,  uint32_t iterations, struct GraphCSR *graph)
{

    uint32_t iter;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double error_total = 0;
    // float init_pr = 1.0f / (float)graph->num_vertices;
    float base_pr = (1.0f - Damp);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif


    float *pageRanks = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Epoch Pull (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext,pageRanks,base_pr)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanks[v] = base_pr;
        pageRanksNext[v] = 0;
    }

    for(iter = 0; iter < iterations; iter++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            epochReorderIncrementCounters(epochReorder, v);

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = sorted_edges_array[j];
                nodeIncomingPR += riDividedOnDiClause[u]; // pageRanks[v]/graph->vertices[v].out_degree;
                atomicEpochReorderIncrementCounters( epochReorder, u);
            }

            // epochReorderIncrementCounters(epochReorder,v);
            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(epsilon, pageRanks,pageRanksNext,base_pr) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  pageRanks[v];
            float nextPageRank =  base_pr + (Damp * pageRanksNext[v]);
            pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", iter, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanks[v] = pageRanks[v] / graph->num_vertices;
        sum += pageRanks[v];
    }

    Stop(timer);


    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", iter, sum, error_total, Seconds(timer));
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, iter, Seconds(timer));
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    return pageRanks;
}


uint32_t *epochReorderRecordBFS(struct GraphCSR *graph)
{


    uint32_t *labels;
    uint32_t *labelsInverse;
    uint32_t *degrees;
    uint32_t root = generateRandInt(mt19937var);
    int v;


    labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    labelsInverse = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    degrees = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));



    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        labelsInverse[v] = v;
        degrees[v] = graph->vertices->out_degree[v];

        // printf("%u %u \n",labelsInverse[v],degrees[v] );
    }

    labelsInverse = radixSortEdgesByDegree(degrees, labelsInverse, graph->num_vertices);

    //decending order mapping
    // #pragma omp parallel for
    // for(v = 0; v < graph->num_vertices; v++){
    //   labels[labelsInverse[v]] = graph->num_vertices -1 - v;

    //   // printf("%u %u \n",labelsInverse[v],degrees[v] );
    // }

    uint32_t numCounters = 32;
    uint32_t hardThreshold = degrees[graph->num_vertices - 1 ];
    uint32_t softThreshold = hardThreshold / 4;
    uint32_t t1 = 0;
    uint32_t t2 = 1;
    uint32_t nextTerm = 1;

    struct EpochReorder *epochReorder = newEpochReoder(softThreshold, hardThreshold, numCounters, graph->num_vertices);

    // #pragma omp parallel for
    // for(v = graph->num_vertices - 1 ; v > graph->num_vertices - 10; v--)
    for(v = 0 ; v < graph->num_vertices ; v += t1 / 2)
    {
        root = labelsInverse[v];

        epochReorderBreadthFirstSearchGraphCSR( epochReorder, root, graph);
        // hardThreshold = stats->processed_nodes;
        // printf(" -----------------------------------------------------\n");
        // printf("| %-15s | %-30f | \n","SUM total", (stats->processed_nodes*100.0f)/graph->num_vertices);
        // printf(" -----------------------------------------------------\n");
        // printf(" -----------------------------------------------------\n");
        // printf("| %-15s | %-30u | \n","EPOCH ", epochReorder->rrIndex);
        // printf(" -----------------------------------------------------\n");

        nextTerm = t1 + t2;
        t1 = t2;
        t2 = nextTerm;
    }


    // printEpochs(epochReorder);

    free(labels);

    labels = epochReorderCreateLabels(epochReorder);



    freeEpochReorder(epochReorder);
    free(labelsInverse);
    free(degrees);


    return labels;


}



// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************



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

void epochReorderBreadthFirstSearchGraphCSR(struct EpochReorder *epochReorder, uint32_t source, struct GraphCSR *graph)
{

    struct BFSStats *stats = newBFSStatsGraphCSR(graph);
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

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Breadth First Search (SOURCE NODE)");
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
        return;
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
                nf = epochReorderBottomUpStepGraphCSR(epochReorder, graph, bitmapCurr, bitmapNext, stats);
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
            bitmapToArrayQueue( bitmapCurr, sharedFrontierQueue, localFrontierQueues);
            Stop(timer_inner);
            printf("| C  %-12s | %-15s | %-15f | \n", " ", " ", Seconds(timer_inner));

            mf = 1;

        }
        else
        {

            Start(timer_inner);
            mu -= mf;
            mf = epochReorderTopDownStepGraphCSR(epochReorder, graph, sharedFrontierQueue, localFrontierQueues, stats);
            slideWindowArrayQueue(sharedFrontierQueue);
            Stop(timer_inner);

            //stats collection
            stats->time_total +=  Seconds(timer_inner);
            stats->processed_nodes += sharedFrontierQueue->tail - sharedFrontierQueue->head;;
            printf("| TD %-12u | %-15u | %-15f | \n", stats->iteration++, sharedFrontierQueue->tail - sharedFrontierQueue->head, Seconds(timer_inner));

        }



    } // end while
    Stop(timer);
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "No OverHead", stats->processed_nodes, stats->time_total);
    printf(" -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15u | %-15f | \n", "total", stats->processed_nodes, Seconds(timer));
    printf(" -----------------------------------------------------\n");

    // graphCSRReset(graph); // no need to recet once processed_nodes = num vertices we traveresed the whole graph
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

uint32_t epochReorderTopDownStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue, struct ArrayQueue **localFrontierQueues, struct BFSStats *stats)
{



    uint32_t v;
    uint32_t u;
    uint32_t i;
    uint32_t j;
    uint32_t edge_idx;
    uint32_t mf = 0;


    #pragma omp parallel default (none) private(u,v,j,i,edge_idx) shared(stats,epochReorder,localFrontierQueues,graph,sharedFrontierQueue,mf)
    {
        uint32_t t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];


        #pragma omp for reduction(+:mf) schedule(auto)
        for(i = sharedFrontierQueue->head ; i < sharedFrontierQueue->tail; i++)
        {
            v = sharedFrontierQueue->queue[i];
            edge_idx = graph->vertices->edges_idx[v];
            // atomicEpochReorderIncrementCounters( epochReorder, v);
            for(j = edge_idx ; j < (edge_idx + graph->vertices->out_degree[v]) ; j++)
            {
                u = graph->sorted_edges_array->edges_array_dest[j];
                int u_parent = stats->parents[u];
                if(u_parent < 0 )
                {
                    atomicEpochReorderIncrementCounters( epochReorder, u);
                    if(__sync_bool_compare_and_swap(&(stats->parents[u]), u_parent, v))
                    {
                        atomicEpochReorderIncrementCounters( epochReorder, u);
                        enArrayQueue(localFrontierQueue, u);
                        mf +=  -(u_parent);
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

uint32_t epochReorderBottomUpStepGraphCSR(struct EpochReorder *epochReorder, struct GraphCSR *graph, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, struct BFSStats *stats)
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

    #pragma omp parallel for default(none) private(j,u,v,out_degree,edge_idx) shared(stats,epochReorder,bitmapCurr,bitmapNext,graph,vertices,sorted_edges_array) reduction(+:nf) schedule(dynamic, 1024)
    for(v = 0 ; v < graph->num_vertices ; v++)
    {
        out_degree = vertices->out_degree[v];
        if(stats->parents[v] < 0)  // optmization
        {
            edge_idx = vertices->edges_idx[v];

            for(j = edge_idx ; j < (edge_idx + out_degree) ; j++)
            {
                u = sorted_edges_array[j];
                atomicEpochReorderIncrementCounters( epochReorder, u);
                if(getBit(bitmapCurr, u))
                {
                    atomicEpochReorderIncrementCounters( epochReorder, v);
                    stats->parents[v] = u;
                    setBitAtomic(bitmapNext, v);
                    nf++;
                    break;
                }
            }

        }

    }
    return nf;
}



uint32_t *epochReorderCreateLabels(struct EpochReorder *epochReorder)
{

    // uint32_t *labels;
    uint32_t *labelsInverse = NULL;
    uint32_t *histMaps = NULL;
    uint32_t *histValues = NULL;
    uint32_t *histDegree = NULL;
    uint32_t v = 0;
    uint32_t h = 0;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Epoch Create Labels");
    printf(" -----------------------------------------------------\n");

    Start(timer);



    // labels = (uint32_t *) my_malloc(epochReorder->numVertices * sizeof(uint32_t));
    labelsInverse = (uint32_t *) my_malloc(epochReorder->numVertices * sizeof(uint32_t));
    histMaps = (uint32_t *) my_malloc(epochReorder->numVertices * sizeof(uint32_t));
    histValues = (uint32_t *) my_malloc(epochReorder->numVertices * sizeof(uint32_t));
    histDegree = (uint32_t *) my_malloc((epochReorder->numCounters) * sizeof(uint32_t));


    #pragma omp parallel for
    for(v = 0; v < epochReorder->numVertices; v++)
    {
        labelsInverse[v] = v;
        histMaps[v] = 0;
        histValues[v] = 0;
    }

    #pragma omp parallel for
    for(v = 0; v < (epochReorder->numCounters); v++)
    {
        histDegree[v] = 0;
    }


    //find max histogram values
    #pragma omp parallel for private(h) shared(histValues, histMaps, epochReorder)
    for(v = 0; v < epochReorder->numVertices; v++)
    {
        uint32_t maxValue = 0;
        uint32_t maxIndex = UINT_MAX;
        for(h = 0; h < epochReorder->numCounters; h++ )
        {
            if(epochReorder->reuse[h][v] > maxValue)
            {
                maxValue = epochReorder->reuse[h][v];
                maxIndex = h;

            }
        }

        if(maxIndex == UINT_MAX)
            maxIndex = 0;
        // histDegree[maxIndex]++;
        // else
        // {
        // uint32_t random_hist =  generateRandInt(mt19937var) % (epochReorder->numCounters - 1);
        // histDegree[random_hist]++;

        // }

        histMaps[v] = maxIndex;
        histValues[v] = maxValue;
    }

    labelsInverse = radixSortEdgesByEpochs(histValues, histMaps, labelsInverse, epochReorder->numVertices);

    // #pragma omp parallel for
    //  for(v = 0; v < epochReorder->numVertices; v++){
    //     labels[labelsInverse[v]] = v;
    //  }

    // uint32_t Accume = histMaps[0];
    // for(v = 0; v < (epochReorder->numVertices); v++)
    // {

    //     // if(histValues[v])
    //     printf("Rank[%u] EPOCH[%u] FREQ[%u] V[%u] \n", v, histMaps[v], histValues[v], labelsInverse[v]);
    // }


    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Epoch Create Labels Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(histMaps);
    free(histValues);
    free(histDegree);

    return labelsInverse;

}

void printEpochs(struct EpochReorder *epochReorder)
{

    uint32_t v = 0;
    uint32_t h = 0;

    for(v = 0; v < 10; v++)
    {
        printf("v[%u] ", v);
        for(h = 0; h < epochReorder->numCounters; h++ )
        {
            printf("%u[%u] ", h, epochReorder->frequency[h][v]);
        }
        printf("\n");
    }

}

void epochReorderIncrementCounters(struct EpochReorder *epochReorder, uint32_t v)
{

    uint32_t histogramIndex = 0;


    if(epochReorder->hardcounter > epochReorder->hardThreshold)
    {
        epochReorder->hardcounter = 0;
        epochReorder->rrIndex++;
        epochReorder->rrIndex  %= epochReorder->numCounters;

    }


    histogramIndex = epochReorder->rrIndex;
    epochReorder->frequency[histogramIndex][v]++;
    epochReorder->hardcounter++;
    epochReorder->softcounter++;

    if(epochReorder->base_reuse[v] != UINT_MAX)
    {
        epochReorder->reuse[histogramIndex][v] += (epochReorder->softcounter - epochReorder->base_reuse[v]);
                epochReorder->base_reuse[v] = epochReorder->softcounter;
    }
    else
    {
        epochReorder->base_reuse[v] = epochReorder->softcounter;
    }



}

void atomicEpochReorderIncrementCounters(struct EpochReorder *epochReorder, uint32_t v)
{


    uint32_t histogramIndex = 0;

    if(epochAtomicMin(&(epochReorder->hardcounter), epochReorder->hardThreshold))
    {

        epochReorder->rrIndex++;
        epochReorder->rrIndex  %= epochReorder->numCounters;

    }
    else
    {
        #pragma omp atomic update
        epochReorder->hardcounter++;

        #pragma omp atomic update
        epochReorder->softcounter++;
    }

    // if(epochAtomicMin(&(epochReorder->softcounter), epochReorder->softThreshold))
    // {
    //     clearBitmap(epochReorder->recencyBits);
    // }
    // else
    // {
    //     #pragma omp atomic update
    //     epochReorder->softcounter++;
    // }

    #pragma omp atomic read
    histogramIndex = epochReorder->rrIndex;

    // if(!getBit(epochReorder->recencyBits, v))
    // {
    //     setBitAtomic(epochReorder->recencyBits, v);
    #pragma omp atomic update
    epochReorder->frequency[histogramIndex][v]++;



    if(!__sync_bool_compare_and_swap(&(epochReorder->base_reuse[v]), UINT_MAX, epochReorder->base_reuse[v]))
    {

        #pragma omp atomic update
        epochReorder->reuse[histogramIndex][v] += (epochReorder->softcounter - epochReorder->base_reuse[v]);

        #pragma omp atomic write
        epochReorder->base_reuse[v] = epochReorder->softcounter;

    }
    else
    {

        epochReorder->base_reuse[v] = epochReorder->softcounter;

    }




}

void freeEpochReorder(struct EpochReorder *epochReorder)
{

    uint32_t v;

    if(epochReorder)
    {
        freeBitmap(epochReorder->recencyBits);
        for(v = 0; v < epochReorder->numCounters; v++)
        {
            free(epochReorder->frequency[v]);
            free(epochReorder->reuse[v]);
        }
        free( epochReorder->frequency);
        free( epochReorder->reuse);
        free( epochReorder->base_reuse);
        free( epochReorder);
    }
}




void radixSortCountSortEdgesByEpochs (uint32_t **histValues, uint32_t **histValuesTemp, uint32_t **histMaps, uint32_t **histMapsTemp, uint32_t **labels, uint32_t **labelsTemp, uint32_t radix, uint32_t buckets, uint32_t *buckets_count, uint32_t num_vertices)
{

    uint32_t *tempPointer1 = NULL;
    uint32_t *tempPointer2 = NULL;
    uint32_t *tempPointer3 = NULL;
    uint32_t t = 0;
    uint32_t o = 0;
    uint32_t u = 0;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t base = 0;

    #pragma omp parallel default(none) shared(histValues, histValuesTemp, histMaps, histMapsTemp,radix,labels,labelsTemp,buckets,buckets_count, num_vertices) firstprivate(t_id, P, offset_end,offset_start,base,i,j,t,u,o)
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id * (num_vertices / P);


        if(t_id == (P - 1))
        {
            offset_end = offset_start + (num_vertices / P) + (num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (num_vertices / P);
        }


        //HISTOGRAM-KEYS
        for(i = 0; i < buckets; i++)
        {
            buckets_count[(t_id * buckets) + i] = 0;
        }


        for (i = offset_start; i < offset_end; i++)
        {
            u = (*histMaps)[i];
            t = (u >> (radix * 8)) & 0xff;
            buckets_count[(t_id * buckets) + t]++;
        }


        #pragma omp barrier


        // SCAN BUCKETS
        if(t_id == 0)
        {
            for(i = 0; i < buckets; i++)
            {
                for(j = 0 ; j < P; j++)
                {
                    t = buckets_count[(j * buckets) + i];
                    buckets_count[(j * buckets) + i] = base;
                    base += t;
                }
            }
        }


        #pragma omp barrier

        //RANK-AND-PERMUTE
        for (i = offset_start; i < offset_end; i++)         /* radix sort */
        {
            u = (*histMaps)[i];
            t = (u >> (radix * 8)) & 0xff;
            o = buckets_count[(t_id * buckets) + t];
            (*histMapsTemp)[o] = (*histMaps)[i];
            (*histValuesTemp)[o] = (*histValues)[i];
            (*labelsTemp)[o] = (*labels)[i];
            buckets_count[(t_id * buckets) + t]++;

        }

    }

    tempPointer1 = *labels;
    *labels = *labelsTemp;
    *labelsTemp = tempPointer1;


    tempPointer2 = *histValues;
    *histValues = *histValuesTemp;
    *histValuesTemp = tempPointer2;

    tempPointer3 = *histMaps;
    *histMaps = *histMapsTemp;
    *histMapsTemp = tempPointer3;

}


uint32_t *radixSortEdgesByEpochs (uint32_t *histValues, uint32_t *histMaps, uint32_t *labels, uint32_t num_vertices)
{


    // printf("*** START Radix Sort Edges By Source *** \n");

    // struct Graph* graph = graphNew(edgeList->num_vertices, edgeList->num_edges, inverse);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    uint32_t radix = 4;  // 32/8 8 bit radix needs 4 iterations
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t buckets = 256; // 2^radix = 256 buckets
    uint32_t *buckets_count = NULL;

    // omp_set_num_threads(P);

    uint32_t j = 0; //1,2,3 iteration
    uint32_t v = 0;


    uint32_t *histValuesTemp = NULL;
    uint32_t *histMapsTemp = NULL;
    uint32_t *labelsTemp = NULL;

    buckets_count = (uint32_t *) my_malloc(P * buckets * sizeof(uint32_t));
    histValuesTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    histMapsTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    labelsTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < num_vertices; v++)
    {
        histValuesTemp[v] = 0;
        histMapsTemp[v] = UINT_MAX / 2;
    }

    for(j = 0 ; j < radix ; j++)
    {
        radixSortCountSortEdgesByEpochs (&histMaps, &histMapsTemp, &histValues, &histValuesTemp, &labels, &labelsTemp, j, buckets, buckets_count, num_vertices);
    }

    for(j = 0 ; j < radix ; j++)
    {
        radixSortCountSortEdgesByEpochs (&histValues, &histValuesTemp, &histMaps, &histMapsTemp, &labels, &labelsTemp, j, buckets, buckets_count, num_vertices);
    }

    free(buckets_count);
    free(histValuesTemp);
    free(histMapsTemp);
    free(labelsTemp);

    return labels;

}
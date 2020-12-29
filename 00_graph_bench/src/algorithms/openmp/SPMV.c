// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : SPMV.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:34:11
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "timer.h"
#include "myMalloc.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"
#include "reorder.h"
#include "graphConfig.h"

#include "fixedPoint.h"
#include "quantization.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "SPMV.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************


struct SPMVStats *newSPMVStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t v;

    struct SPMVStats *stats = (struct SPMVStats *) my_malloc(sizeof(struct SPMVStats));

    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->vector_output = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->vector_input = (float *) my_malloc(graph->num_vertices * sizeof(float));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->vector_output[v] =  0.0f;
        stats->vector_input[v] =  0.0f;
    }

    return stats;

}
struct SPMVStats *newSPMVStatsGraphGrid(struct GraphGrid *graph)
{

    uint32_t v;

    struct SPMVStats *stats = (struct SPMVStats *) my_malloc(sizeof(struct SPMVStats));

    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->vector_output = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->vector_input = (float *) my_malloc(graph->num_vertices * sizeof(float));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->vector_output[v] =  0.0f;
        stats->vector_input[v] =  0.0f;
    }

    return stats;

}
struct SPMVStats *newSPMVStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    uint32_t v;

    struct SPMVStats *stats = (struct SPMVStats *) my_malloc(sizeof(struct SPMVStats));

    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->vector_output = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->vector_input = (float *) my_malloc(graph->num_vertices * sizeof(float));


    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->vector_output[v] =  0.0f;
        stats->vector_input[v] =  0.0f;
    }

    return stats;

}
struct SPMVStats *newSPMVStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    uint32_t v;

    struct SPMVStats *stats = (struct SPMVStats *) my_malloc(sizeof(struct SPMVStats));

    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->vector_output = (float *) my_malloc(graph->num_vertices * sizeof(float));
    stats->vector_input = (float *) my_malloc(graph->num_vertices * sizeof(float));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->vector_output[v] =  0.0f;
        stats->vector_input[v] =  0.0f;
    }

    return stats;

}

void freeSPMVStats(struct SPMVStats *stats)
{

    if(stats)
    {
        if(stats->vector_output)
            free(stats->vector_output);
        if(stats->vector_input)
            free(stats->vector_input);

        free(stats);
    }

}

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphGrid( struct Arguments *arguments, struct GraphGrid *graph)
{

    struct SPMVStats *stats = NULL;

    switch (arguments->pushpull)
    {
    case 0: // push
        stats = SPMVPullRowGraphGrid( arguments, graph);
        break;
    case 1: // pull
        stats = SPMVPushColumnGraphGrid( arguments, graph);
        break;
    case 2: // pull
        stats = SPMVPullRowFixedPointGraphGrid( arguments, graph);
        break;
    case 3: // push
        stats = SPMVPushColumnFixedPointGraphGrid( arguments, graph);
        break;
    default:// pull
        stats = SPMVPullRowGraphGrid( arguments, graph);
        break;
    }

    return stats;

}
struct SPMVStats *SPMVPullRowGraphGrid( struct Arguments *arguments, struct GraphGrid *graph)
{

    uint32_t v;
    double sum = 0.0;

    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct SPMVStats *stats = newSPMVStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-Row");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->grid->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->grid->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        uint32_t i;
        // #pragma omp parallel for private(i) schedule (dynamic,arguments->algo_numThreads)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t j;
            #pragma omp parallel for private(j) schedule (dynamic,arguments->algo_numThreads)
            for (j = 0; j < totalPartitions; ++j)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                float weight = 0.0001f;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

#if WEIGHTED
                    weight = partition->edgeList->edges_array_weight[k];
#endif

                    // #pragma omp atomic update
                    // __sync_fetch_and_add(&stats->vector_output[dest],(weight * stats->vector_input[src]));
                    // addAtomicFloat(&stats->vector_output[dest], (weight * stats->vector_input[src])

                    // #pragma omp atomic update
                    stats->vector_output[dest] +=  (weight * stats->vector_input[src]);
                }
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;


}
struct SPMVStats *SPMVPushColumnGraphGrid( struct Arguments *arguments, struct GraphGrid *graph)
{
    uint32_t v;
    double sum = 0.0;

    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct SPMVStats *stats = newSPMVStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-Column");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->grid->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->grid->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        uint32_t j;
        #pragma omp parallel for private(j) schedule (dynamic,arguments->algo_numThreads)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
        {
            uint32_t i;
            // #pragma omp parallel for private(i) schedule (dynamic,arguments->algo_numThreads)
            for (i = 0; i < totalPartitions; ++i)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                float weight = 0.0001f;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

#if WEIGHTED
                    weight = partition->edgeList->edges_array_weight[k];
#endif

                    // #pragma omp atomic update
                    stats->vector_output[dest] +=  (weight * stats->vector_input[src]);
                }
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;

}
struct SPMVStats *SPMVPullRowFixedPointGraphGrid( struct Arguments *arguments, struct GraphGrid *graph)
{

    uint32_t v;
    double sum = 0.0;

    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct SPMVStats *stats = newSPMVStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-Row Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->grid->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->grid->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }


    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        uint32_t i;
        // #pragma omp parallel for private(i) schedule (dynamic,arguments->algo_numThreads)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t j;
            #pragma omp parallel for private(j) schedule (dynamic,arguments->algo_numThreads)
            for (j = 0; j < totalPartitions; ++j)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                uint64_t weight = DoubleToFixed64(0.0001f);
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

#if WEIGHTED
                    weight = DoubleToFixed64(partition->edgeList->edges_array_weight[k]);
#endif
                    // #pragma omp atomic update
                    vector_output[dest] += MULFixed64V1(weight, vector_input[src]);
                }
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;


}
struct SPMVStats *SPMVPushColumnFixedPointGraphGrid( struct Arguments *arguments, struct GraphGrid *graph)
{

    uint32_t v;
    double sum = 0.0;

    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct SPMVStats *stats = newSPMVStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-Column Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->grid->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->grid->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }


    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        uint32_t j;
        #pragma omp parallel for private(j) schedule (dynamic,arguments->algo_numThreads)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
        {
            uint32_t i;
            // #pragma omp parallel for private(i) schedule (dynamic,arguments->algo_numThreads)
            for (i = 0; i < totalPartitions; ++i)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                uint64_t weight = DoubleToFixed64(0.0001f);
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

#if WEIGHTED
                    weight = DoubleToFixed64(partition->edgeList->edges_array_weight[k]);
#endif

                    // #pragma omp atomic update
                    vector_output[dest] += MULFixed64V1(weight, vector_input[src]);
                }
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;

}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphCSR( struct Arguments *arguments, struct GraphCSR *graph)
{

    struct SPMVStats *stats = NULL;

    switch (arguments->pushpull)
    {

    case 0: // pull
        stats = SPMVPullGraphCSR( arguments, graph);
        break;
    case 1: // push
        stats = SPMVPushGraphCSR( arguments, graph);
        break;
    case 2: // pull
        stats = SPMVPullFixedPointGraphCSR( arguments, graph);
        break;
    case 3: // push
        stats = SPMVPushFixedPointGraphCSR( arguments, graph);
        break;
    default:// pull
        stats = SPMVPullGraphCSR( arguments, graph);
        break;
    }

    return stats;

}
struct SPMVStats *SPMVPullGraphCSR( struct Arguments *arguments, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    double sum = 0.0;

    struct SPMVStats *stats = newSPMVStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
#if WEIGHTED
    float *edges_array_weight = NULL;
#endif

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->inverse_sorted_edges_array->edges_array_weight;
#endif
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->sorted_edges_array->edges_array_weight;
#endif
#endif

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->vertices->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src ;
            uint32_t dest = v;
            float weight = 0.0001f;
            degree = vertices->out_degree[dest];
            edge_idx = vertices->edges_idx[dest];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                src = EXTRACT_VALUE(sorted_edges_array[j]);
#if WEIGHTED
                weight = edges_array_weight[j];
#endif
                stats->vector_output[dest] +=  (weight * stats->vector_input[src]);
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    return stats;

}
struct SPMVStats *SPMVPushGraphCSR( struct Arguments *arguments, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    double sum = 0.0;

    struct SPMVStats *stats = newSPMVStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
#if WEIGHTED
    float *edges_array_weight = NULL;
#endif

    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->sorted_edges_array->edges_array_weight;
#endif


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->vertices->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            float weight = 0.0001f;
            degree = vertices->out_degree[src];
            edge_idx = vertices->edges_idx[src];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                dest = EXTRACT_VALUE(sorted_edges_array[j]);
#if WEIGHTED
                weight = edges_array_weight[j];
#endif

                #pragma omp atomic update
                stats->vector_output[dest] += (weight * stats->vector_input[src]);
            }

        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;

}

struct SPMVStats *SPMVPullFixedPointGraphCSR( struct Arguments *arguments, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t w;
    double sum = 0.0;

    struct SPMVStats *stats = newSPMVStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint32_t *vector_input = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    uint32_t *edges_array_weight_fixedPoint = (uint32_t *) my_malloc(graph->num_edges * sizeof(uint32_t));


    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
#if WEIGHTED
    float *edges_array_weight = NULL;
#endif

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->inverse_sorted_edges_array->edges_array_weight;
#endif
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->sorted_edges_array->edges_array_weight;
#endif
#endif

    #pragma omp parallel for
    for (w = 0; w < graph->num_edges ; ++w)
    {
#if WEIGHTED
        edges_array_weight_fixedPoint[w] = FloatToFixed32(edges_array_weight[w]);
#else
        edges_array_weight_fixedPoint[w] = FloatToFixed32(0.0001f);
#endif
    }


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->vertices->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;


    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = FloatToFixed32(stats->vector_input[v]);
    }

    Start(timer);

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src;
            uint32_t dest = v;
            float weight = FloatToFixed32(0.0001f);
            degree = vertices->out_degree[dest];
            edge_idx = vertices->edges_idx[dest];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                src = EXTRACT_VALUE(sorted_edges_array[j]);
#if WEIGHTED
                weight = edges_array_weight_fixedPoint[j];
#endif
                vector_output[dest] += MULFixed32V1(weight, vector_input[src]); // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed32ToFloat(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);
    free(edges_array_weight_fixedPoint);

    return stats;

}
struct SPMVStats *SPMVPushFixedPointGraphCSR( struct Arguments *arguments, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    double sum = 0.0;

    struct SPMVStats *stats = newSPMVStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
#if WEIGHTED
    float *edges_array_weight = NULL;
#endif

    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#if WEIGHTED
    edges_array_weight = graph->sorted_edges_array->edges_array_weight;
#endif


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices->out_degree[v])
            stats->vector_input[v] =  (1.0f / graph->vertices->out_degree[v]);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            uint64_t weight = DoubleToFixed64(0.0001f);
            degree = vertices->out_degree[src];
            edge_idx = vertices->edges_idx[src];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                dest = EXTRACT_VALUE(sorted_edges_array[j]);
#if WEIGHTED
                weight = DoubleToFixed64(edges_array_weight[j]);
#endif

                #pragma omp atomic update
                vector_output[dest] += MULFixed64V1(weight, vector_input[src]);
            }

        }
        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;

}

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjArrayList( struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    struct SPMVStats *stats = NULL;

    switch (arguments->pushpull)
    {

    case 0: // pull
        stats = SPMVPullGraphAdjArrayList( arguments, graph);
        break;
    case 1: // push
        stats = SPMVPushGraphAdjArrayList( arguments, graph);
        break;
    case 2: // pull
        stats = SPMVPullFixedPointGraphAdjArrayList( arguments, graph);
        break;
    case 3: // push
        stats = SPMVPushFixedPointGraphAdjArrayList( arguments, graph);
        break;
    default:// push
        stats = SPMVPullGraphAdjArrayList( arguments, graph);
        break;
    }


    return stats;

}
struct SPMVStats *SPMVPullGraphAdjArrayList( struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct EdgeList *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src;
            uint32_t dest = v;
            float weight = 0.0001f;

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[dest].inNodes;
            degree = graph->vertices[dest].in_degree;
#else
            Nodes = graph->vertices[dest].outNodes;
            degree = graph->vertices[dest].out_degree;
#endif

            for(j = 0 ; j < (degree) ; j++)
            {
                src = Nodes->edges_array_dest[j];
#if WEIGHTED
                weight = Nodes->edges_array_weight[j];
#endif
                stats->vector_output[dest] +=  (weight * stats->vector_input[src]); // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;


}
struct SPMVStats *SPMVPushGraphAdjArrayList( struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct EdgeList *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            float weight = 0.0001f;

            Nodes = graph->vertices[src].outNodes;
            degree = graph->vertices[src].out_degree;

            for(j = 0 ; j <  (degree) ; j++)
            {
                dest =  Nodes->edges_array_dest[j];
#if WEIGHTED
                weight = Nodes->edges_array_weight[j];
#endif

                #pragma omp atomic update
                stats->vector_output[dest] += (weight * stats->vector_input[src]);
            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;

}

struct SPMVStats *SPMVPullFixedPointGraphAdjArrayList( struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct EdgeList *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src;
            uint32_t dest = v;
            uint64_t weight = DoubleToFixed64(0.0001f);

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[dest].inNodes;
            degree = graph->vertices[dest].in_degree;
#else
            Nodes = graph->vertices[dest].outNodes;
            degree = graph->vertices[dest].out_degree;
#endif

            for(j = 0 ; j < (degree) ; j++)
            {
                src = Nodes->edges_array_dest[j];


#if WEIGHTED
                weight = DoubleToFixed64(Nodes->edges_array_weight[j]);
#endif

                vector_output[dest] +=  MULFixed64V1(weight, vector_input[src]); // stats->pageRanks[v]/graph->vertices[v].out_degree;

            }

        }



        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);

    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }


    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;

}
struct SPMVStats *SPMVPushFixedPointGraphAdjArrayList( struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct EdgeList *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            uint64_t weight = DoubleToFixed64(0.0001f);

            Nodes = graph->vertices[src].outNodes;
            degree = graph->vertices[src].out_degree;

            for(j = 0 ; j <  (degree) ; j++)
            {
                dest =  Nodes->edges_array_dest[j];
#if WEIGHTED
                weight = DoubleToFixed64(Nodes->edges_array_weight[j]);
#endif
                #pragma omp atomic update
                vector_output[dest] += MULFixed64V1(weight, vector_input[src]);
            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }


    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;

}

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct SPMVStats *SPMVGraphAdjLinkedList( struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    struct SPMVStats *stats = NULL;

    switch (arguments->pushpull)
    {

    case 0: // pull
        stats = SPMVPullGraphAdjLinkedList( arguments, graph);
        break;
    case 1: // push
        stats = SPMVPushGraphAdjLinkedList( arguments, graph);
        break;
    case 2: // pull
        stats = SPMVPullFixedPointGraphAdjLinkedList( arguments, graph);
        break;
    case 3: // push
        stats = SPMVPushFixedPointGraphAdjLinkedList( arguments, graph);
        break;
    default:// push
        stats = SPMVPullGraphAdjLinkedList( arguments, graph);
        break;
    }


    return stats;

}
struct SPMVStats *SPMVPullGraphAdjLinkedList( struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct AdjLinkedListNode *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src;
            uint32_t dest = v;
            float weight = 0.0001f;

#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[dest].inNodes;
            degree = graph->vertices[dest].in_degree;
#else
            Nodes = graph->vertices[dest].outNodes;
            degree = graph->vertices[dest].out_degree;
#endif
            for(j = 0 ; j < (degree) ; j++)
            {
                src =  Nodes->dest;
#if WEIGHTED
                weight = Nodes->weight;
#endif
                Nodes = Nodes->next;

                stats->vector_output[dest] +=  (weight * stats->vector_input[src]); // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }


        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;
}
struct SPMVStats *SPMVPushGraphAdjLinkedList( struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct AdjLinkedListNode *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            float weight = 0.0001f;

            Nodes = graph->vertices[src].outNodes;
            degree = graph->vertices[src].out_degree;

            for(j = 0 ; j <  (degree) ; j++)
            {

                dest =  Nodes->dest;
#if WEIGHTED
                weight = Nodes->weight;
#endif
                Nodes = Nodes->next;

                #pragma omp atomic update
                stats->vector_output[dest] += (weight * stats->vector_input[src]);
            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    return stats;
}
struct SPMVStats *SPMVPullFixedPointGraphAdjLinkedList( struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct AdjLinkedListNode *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PULL Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src;
            uint32_t dest = v;
            uint64_t weight = DoubleToFixed64(0.0001f);
#if DIRECTED // will look at the other neighbours if directed by using inverese edge list
            Nodes = graph->vertices[dest].inNodes;
            degree = graph->vertices[dest].in_degree;
#else
            Nodes = graph->vertices[dest].outNodes;
            degree = graph->vertices[dest].out_degree;
#endif
            for(j = 0 ; j < (degree) ; j++)
            {
                src =  Nodes->dest;

#if WEIGHTED
                weight = DoubleToFixed64(Nodes->weight);
#endif
                Nodes = Nodes->next;

                vector_output[dest] +=  MULFixed64V1(weight, vector_input[src]); // stats->pageRanks[v]/graph->vertices[v].out_degree;

            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }


    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;
}
struct SPMVStats *SPMVPushFixedPointGraphAdjLinkedList( struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    uint32_t degree;
    double sum = 0.0;
    struct AdjLinkedListNode *Nodes;

    struct SPMVStats *stats = newSPMVStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    uint64_t *vector_input = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *vector_output = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting SPMV-PUSH Fixed-Point");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    //assume any vector input for benchamrking purpose.
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(graph->vertices[v].out_degree)
            stats->vector_input[v] =  (1.0f / graph->vertices[v].out_degree);
        else
            stats->vector_input[v] = 0.001f;
    }

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        vector_output[v] = 0;
        vector_input[v] = DoubleToFixed64(stats->vector_input[v]);
    }

    Start(timer);
    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;
            uint64_t weight = DoubleToFixed64(0.0001f);

            Nodes = graph->vertices[src].outNodes;
            degree = graph->vertices[src].out_degree;

            for(j = 0 ; j <  (degree) ; j++)
            {
                dest =  Nodes->dest;

#if WEIGHTED
                weight = DoubleToFixed64(Nodes->weight);
#endif
                Nodes = Nodes->next;

                #pragma omp atomic update
                vector_output[dest] += MULFixed64V1(weight, vector_input[src]);
            }

        }


        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));

    }// end iteration loop

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->vector_output[v] = Fixed64ToDouble(vector_output[v]);
    }


    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {

        sum += ((int)(stats->vector_output[v] * 100 + .5) / 100.0);
    }


    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Sum", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15lf | %-15f | \n", stats->iterations, sum, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(vector_output);
    free(vector_input);

    return stats;
}

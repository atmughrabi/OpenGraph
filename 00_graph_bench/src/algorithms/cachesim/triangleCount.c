// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : triangleCount.c
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

#include "graphConfig.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "triangleCount.h"


struct TCStats *newTCStatsGraphCSR(struct GraphCSR *graph)
{
    uint32_t v;
    struct TCStats *stats = (struct TCStats *) my_malloc(sizeof(struct TCStats));

    stats->total_counts = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->counts = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->counts[v] =  0;
    }

    return stats;

}
struct TCStats *newTCStatsGraphGrid(struct GraphGrid *graph)
{
    uint32_t v;
    struct TCStats *stats = (struct TCStats *) my_malloc(sizeof(struct TCStats));

    stats->total_counts = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->counts = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->counts[v] =  0;
    }
    return stats;
}
struct TCStats *newTCStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{
    uint32_t v;
    struct TCStats *stats = (struct TCStats *) my_malloc(sizeof(struct TCStats));

    stats->total_counts = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->counts = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->counts[v] =  0;
    }
    return stats;

}
struct TCStats *newTCStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    struct TCStats *stats = (struct TCStats *) my_malloc(sizeof(struct TCStats));

    stats->total_counts = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->counts = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->counts[v] =  0;
    }
    return stats;

}
void freeTCStats(struct TCStats *stats)
{

    if(stats)
    {
        if(stats->counts)
            free(stats->counts);

        free(stats);
    }

}

// ********************************************************************************************
// ***************                  Helper Functions                             **************
// ********************************************************************************************

uint32_t minTwoNodes(uint32_t node_v, uint32_t node_u, uint32_t degree_v, uint32_t degree_u)
{

    if(degree_v < degree_u)
        return node_v;
    else
        return node_u;

}

uint32_t maxTwoNodes(uint32_t node_v, uint32_t node_u, uint32_t degree_v, uint32_t degree_u)
{

    if(degree_u > degree_v)
        return node_u;
    else
        return node_v;

}

uint32_t countIntersectionsBinarySearch(uint32_t u, uint32_t v, struct GraphCSR *graph)
{

    uint32_t count = 0;

    uint32_t degree_iter = graph->vertices->out_degree[v];
    uint32_t edge_idx_iter = graph->vertices->edges_idx[v];

    uint32_t degree_comp = graph->vertices->out_degree[u];
    uint32_t edge_idx_comp = graph->vertices->edges_idx[u];

    uint32_t iter;

    for(iter = edge_idx_iter ; iter < (edge_idx_iter + degree_iter); iter++ )
    {

        uint32_t u_iter = graph->sorted_edges_array->edges_array_dest[iter];
        if(u_iter > v)
            break;

        uint32_t bottom = 0;
        uint32_t top = degree_comp;
        uint32_t mid = (top + bottom) >> 1;
        uint32_t v_comp = graph->sorted_edges_array->edges_array_dest[edge_idx_comp + mid];

        while( bottom < (top - 1))
        {

            if(u_iter < v_comp)
            {
                top = mid;

            }
            else if ( u_iter > v_comp)
            {
                bottom = mid;

            }
            else
            {
                count++;
                break;
            }


            mid = (top + bottom) >> 1;
            v_comp = graph->sorted_edges_array->edges_array_dest[edge_idx_comp + mid];
            u_iter = graph->sorted_edges_array->edges_array_dest[iter];

        }

        if((top - 1) == 0 && u_iter == v_comp)
            count++;


    }
    return count;
}


// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct TCStats *triangleCountGraphCSR(uint32_t pushpull, struct GraphCSR *graph)
{
    struct TCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // basic slow
        stats = triangleCountBasicGraphCSR(graph);
        break;
    case 1: // pull
        stats = triangleCountPullGraphCSR(graph);
        break;
    case 2: // push
        stats = triangleCountPushGraphCSR(graph);
        break;
    case 3: // With binary intersection
        stats = triangleCountBinaryIntersectionGraphCSR(graph);
        break;
    default:// pull
        stats = triangleCountPullGraphCSR(graph);
        break;
    }

    return stats;

}
struct TCStats *triangleCountBasicGraphCSR(struct GraphCSR *graph)
{

    uint32_t u;
    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count-basic");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    Start(timer);
    #pragma omp parallel for shared(stats) schedule(dynamic, 128)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t degree_u = graph->vertices->out_degree[u];
        uint32_t edge_idx_u = graph->vertices->edges_idx[u];
        uint32_t v;

        for(v = edge_idx_u; v < (edge_idx_u + degree_u) ; v++)
        {
            uint32_t node_v = graph->sorted_edges_array->edges_array_dest[v];
            uint32_t degree_v = graph->vertices->out_degree[node_v];
            uint32_t edge_idx_v = graph->vertices->edges_idx[node_v];
            uint32_t w;

            uint32_t degree_iter = graph->vertices->out_degree[u];
            uint32_t edge_idx_iter = graph->vertices->edges_idx[u];
            uint32_t iter;

            for(w = edge_idx_v; w < (edge_idx_v + degree_v) ; w++)
            {
                uint32_t node_w = graph->sorted_edges_array->edges_array_dest[w];
                uint32_t node_iter = graph->sorted_edges_array->edges_array_dest[edge_idx_iter];

                for(iter = edge_idx_iter; iter < (edge_idx_iter + degree_iter) ; iter++)
                {
                    node_iter = graph->sorted_edges_array->edges_array_dest[iter];

                    if(node_iter == node_w)
                        stats->counts[u]++;
                }
            }
        }
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    #pragma omp parallel for default(none) reduction (+ : counts) private(u) shared(stats)
    for(u = 0; u < stats->num_vertices; u++)
    {
        counts += stats->counts[u];
    }

    stats->total_counts = counts / 6;

    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    return stats;

}

struct TCStats *triangleCountPullGraphCSR(struct GraphCSR *graph)
{

    uint32_t u;
    uint64_t counts = 0;
    uint64_t steps = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count-PULL");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    Start(timer);
    #pragma omp parallel for shared(stats) reduction(+:counts) schedule(dynamic, 128)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t degree_u = graph->vertices->out_degree[u];
        uint32_t edge_idx_u = graph->vertices->edges_idx[u];
        uint32_t v;

        steps++;
        for(v = edge_idx_u; v < (edge_idx_u + degree_u) ; v++)
        {
            uint32_t node_v = graph->sorted_edges_array->edges_array_dest[v];
            uint32_t degree_v = graph->vertices->out_degree[node_v];

            if(node_v > u)
                break;

            uint32_t edge_idx_v = graph->vertices->edges_idx[node_v];
            uint32_t w;


            uint32_t degree_iter = graph->vertices->out_degree[u];
            uint32_t edge_idx_iter = graph->vertices->edges_idx[u];
            uint32_t iter;

            for(w = edge_idx_v; w < (edge_idx_v + degree_v) ; w++)
            {

                uint32_t node_w = graph->sorted_edges_array->edges_array_dest[w];
                if(node_w > node_v)
                    break;

                uint32_t node_iter = graph->sorted_edges_array->edges_array_dest[edge_idx_iter];



                for(iter = edge_idx_iter; iter < (edge_idx_iter + degree_iter) ; iter++)
                {
                    node_iter = graph->sorted_edges_array->edges_array_dest[iter];

                    if(node_iter >= node_w)
                        break;
                }


                if(node_w == node_iter)
                {
                    counts++;
                }
            }
        }
    }
    Stop(timer);
    stats->time_total = Seconds(timer);

    stats->total_counts = counts;

    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    return stats;

}
struct TCStats *triangleCountPushGraphCSR(struct GraphCSR *graph)
{

    uint32_t u;
    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count-PUSH");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    Start(timer);
    #pragma omp parallel for shared(stats) schedule(dynamic, 128)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t degree_u = graph->vertices->out_degree[u];
        uint32_t edge_idx_u = graph->vertices->edges_idx[u];
        uint32_t v;

        for(v = edge_idx_u; v < (edge_idx_u + degree_u) ; v++)
        {
            uint32_t node_v = graph->sorted_edges_array->edges_array_dest[v];

            if(node_v > u)
                break;

            uint32_t degree_v = graph->vertices->out_degree[node_v];
            uint32_t edge_idx_v = graph->vertices->edges_idx[node_v];
            uint32_t w;

            uint32_t degree_iter = graph->vertices->out_degree[u];
            uint32_t edge_idx_iter = graph->vertices->edges_idx[u];
            uint32_t iter;

            for(w = edge_idx_v; w < (edge_idx_v + degree_v) ; w++)
            {

                uint32_t node_w = graph->sorted_edges_array->edges_array_dest[w];

                if(node_w > node_v)
                    break;

                uint32_t node_iter = graph->sorted_edges_array->edges_array_dest[edge_idx_iter];

                for(iter = edge_idx_iter; iter < (edge_idx_iter + degree_iter) ; iter++)
                {
                    node_iter = graph->sorted_edges_array->edges_array_dest[iter];

                    if(node_iter >= node_w)
                        break;
                }

                if(node_w == node_iter)
                {
                    #pragma omp atomic update
                    stats->counts[node_w]++;
                }
            }
        }
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    #pragma omp parallel for default(none) reduction (+ : counts) private(u) shared(stats)
    for(u = 0; u < stats->num_vertices; u++)
    {
        counts += stats->counts[u];
    }

    stats->total_counts = counts;

    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    return stats;

}


struct TCStats *triangleCountBinaryIntersectionGraphCSR(struct GraphCSR *graph)
{

    uint32_t u;
    uint64_t counts = 0;
    uint64_t steps = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Binary-Intersection");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    Start(timer);
    #pragma omp parallel for shared(stats) reduction(+:counts) schedule(dynamic, 128)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t degree_u = graph->vertices->out_degree[u];
        uint32_t edge_idx_u = graph->vertices->edges_idx[u];
        uint32_t v;

        steps++;
        for(v = edge_idx_u; v < (edge_idx_u + degree_u) ; v++)
        {
            uint32_t node_v = graph->sorted_edges_array->edges_array_dest[v];

            if(node_v > u)
                break;
            counts += countIntersectionsBinarySearch(u, node_v, graph);
        }
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    stats->total_counts = counts;

    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    return stats;

}

// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct TCStats *triangleCountGraphGrid(uint32_t pushpull, struct GraphGrid *graph)
{
    struct TCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = triangleCountRowGraphGrid(graph);
        break;
    case 1: // push
        stats = triangleCountColumnGraphGrid(graph);
        break;
    default:// pull
        stats = triangleCountRowGraphGrid(graph);
        break;
    }

    return stats;

}
struct TCStats *triangleCountRowGraphGrid(struct GraphGrid *graph)
{

    
    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;

}
struct TCStats *triangleCountColumnGraphGrid(struct GraphGrid *graph)
{

    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;

}

// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjArrayList(uint32_t pushpull, struct GraphAdjArrayList *graph)
{
    struct TCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = triangleCountPullGraphAdjArrayList(graph);
        break;
    case 1: // push
        stats = triangleCountPullGraphAdjArrayList(graph);
        break;
    default:// pull
        stats = triangleCountPullGraphAdjArrayList(graph);
        break;
    }

    return stats;
}
struct TCStats *triangleCountPullGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;

}
struct TCStats *triangleCountPushGraphAdjArrayList(struct GraphAdjArrayList *graph)
{
    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;
}

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct TCStats *triangleCountGraphAdjLinkedList(uint32_t pushpull, struct GraphAdjLinkedList *graph)
{
    struct TCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = triangleCountPullGraphAdjLinkedList(graph);
        break;
    case 1: // push
        stats = triangleCountPushGraphAdjLinkedList(graph);
        break;
    default:// pull
        stats = triangleCountPullGraphAdjLinkedList(graph);
        break;
    }

    return stats;
}
struct TCStats *triangleCountPullGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{
    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;

}
struct TCStats *triangleCountPushGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    uint64_t counts = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Triangle Count To Be Implemented");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Triangle Counts", "Time (S)");
    printf(" -----------------------------------------------------\n");

    struct TCStats *stats = newTCStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    Stop(timer);
    stats->time_total = Seconds(timer);
    stats->total_counts = counts;
    printf("| %-21lu | %-27f | \n", stats->total_counts, stats->time_total);
    printf(" -----------------------------------------------------\n");
    return stats;

}
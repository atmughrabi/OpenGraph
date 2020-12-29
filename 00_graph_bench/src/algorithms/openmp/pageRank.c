// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : pageRank.c
// Create : 2019-09-28 14:41:30
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
#include "worklist.h"

#include "graphConfig.h"

#include "fixedPoint.h"
#include "quantization.h"
#include "reorder.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "pageRank.h"

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************


struct PageRankStats *newPageRankStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));;
    stats->pageRanks = (float *) my_malloc(graph->num_vertices * sizeof(float));;


    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->pageRanks[v] =  stats->base_pr;
        stats->realRanks[v] =  v;
    }


    return stats;

}
struct PageRankStats *newPageRankStatsGraphGrid(struct GraphGrid *graph)
{

    uint32_t v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));;
    stats->pageRanks = (float *) my_malloc(graph->num_vertices * sizeof(float));;


    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->pageRanks[v] =  stats->base_pr;
        stats->realRanks[v] =  v;
    }


    return stats;


}
struct PageRankStats *newPageRankStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    uint32_t v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));;
    stats->pageRanks = (float *) my_malloc(graph->num_vertices * sizeof(float));;


    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->pageRanks[v] =  stats->base_pr;
        stats->realRanks[v] =  v;
    }


    return stats;


}
struct PageRankStats *newPageRankStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    uint32_t v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));;
    stats->pageRanks = (float *) my_malloc(graph->num_vertices * sizeof(float));;


    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->pageRanks[v] =  stats->base_pr;
        stats->realRanks[v] =  v;
    }


    return stats;


}

void freePageRankStats(struct PageRankStats *stats)
{
    if(stats)
    {
        if(stats->realRanks)
            free(stats->realRanks);
        if(stats->pageRanks)
            free(stats->pageRanks);
        free(stats);
    }
}



// ********************************************************************************************
// ***************          Auxilary functions                                   **************
// ********************************************************************************************

void addAtomicFloat(float *num, float value)
{

    float newV, oldV;
    uint32_t *lnewV;
    uint32_t *loldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
        loldV = (uint32_t *)&oldV;
        lnewV = (uint32_t *)&newV;
    }
    while(!__sync_bool_compare_and_swap((uint32_t *)num, *(loldV), *(lnewV)));

}


void addAtomicDouble(double *num, double value)
{

    double newV, oldV;
    uint64_t *lnewV;
    uint64_t *loldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
        loldV = (uint64_t *)&oldV;
        lnewV = (uint64_t *)&newV;
    }
    while(!__sync_bool_compare_and_swap((uint64_t *)num, *(loldV), *(lnewV)));

}

void setAtomic(uint64_t *num, uint64_t value)
{


    uint64_t newV, oldV;

    do
    {
        oldV = *num;
        newV = value;
    }
    while(!__sync_bool_compare_and_swap(num, oldV, newV));

}

void addAtomicFixedPoint(uint64_t *num, uint64_t value)
{

    uint64_t newV, oldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
    }
    while(!__sync_bool_compare_and_swap(num, oldV, newV));

}

void pageRankPrint(float *pageRankArray, uint32_t num_vertices)
{
    uint32_t v;
    for(v = 0; v < num_vertices; v++)
    {
        printf("Rank[%d]=%f \n", v, pageRankArray[v]);
    }
}

// ********************************************************************************************
// ***************          GRID DataStructure               **************
// ********************************************************************************************


// function STREAMVERTICES(Fv,F)
//  Sum = 0
//    for each vertex do
//      if F(vertex) then
//        Sum += Fv(edge)
//      end if
//    end for
//  return Sum
// end function

// function STREAMEDGES(Fe,F)
//  Sum = 0
//    for each active block do >> block with active edges
//      for each edge ∈ block do
//        if F(edge.source) then
//          Sum += Fe(edge)
//        end if
//      end for
//    end for
//  return Sum
// end function
//we assume that the edges are not sorted in each partition

struct PageRankStats  *pageRankGraphGrid(struct Arguments *arguments, struct GraphGrid *graph)
{

    struct PageRankStats *stats = NULL;

    switch (arguments->pushpull)
    {
    case 0: // push
        stats = pageRankPullRowGraphGrid(arguments, graph);
        break;
    case 1: // pull
        stats = pageRankPushColumnGraphGrid(arguments, graph);
        break;
    case 2: // pull
        stats = pageRankPullRowFixedPointGraphGrid(arguments, graph);
        break;
    case 3: // push
        stats = pageRankPushColumnFixedPointGraphGrid(arguments, graph);
        break;
    default:// pull
        stats = pageRankPullRowGraphGrid(arguments, graph);
        break;
    }

    return stats;

}


struct PageRankStats *pageRankPullRowGraphGrid(struct Arguments *arguments, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;

    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Row (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->grid->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->grid->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        // pageRankStreamEdgesGraphGridRowWise(graph, riDividedOnDiClause, pageRanksNext);

        uint32_t i;
        // #pragma omp parallel for private(i)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t j;
            #pragma omp parallel for private(j)
            for (j = 0; j < totalPartitions; ++j)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

                    // #pragma omp atomic update
                    // __sync_fetch_and_add(&pageRanksNext[dest],riDividedOnDiClause[src]);
                    // addAtomicFloat(float *num, float value)

                    // #pragma omp atomic update
                    pageRanksNext[dest] +=  riDividedOnDiClause[src];
                }
            }
        }


        #pragma omp parallel for private(v) shared(arguments, pageRanksNext, stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}



struct PageRankStats *pageRankPullRowFixedPointGraphGrid(struct Arguments *arguments, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Row FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->grid->out_degree[v])
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->grid->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        // pageRankStreamEdgesGraphGridRowWise(graph, riDividedOnDiClause, pageRanksNext);

        uint32_t i;
        // #pragma omp parallel for private(i)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t j;
            #pragma omp parallel for private(j)
            for (j = 0; j < totalPartitions; ++j)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

                    // #pragma omp atomic update
                    pageRanksNext[dest] +=  riDividedOnDiClause[src];
                }
            }
        }


        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}



/******************************************************************/



struct PageRankStats *pageRankPushColumnGraphGrid(struct Arguments *arguments, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Col (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->grid->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->grid->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        // pageRankStreamEdgesGraphGridRowWise(graph, riDividedOnDiClause, pageRanksNext);

        uint32_t j;
        #pragma omp parallel for private(j)
        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t i;

            // #pragma omp parallel for private(i) // iterate over partitions columnwise
            for (i = 0; i < totalPartitions; ++i)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

                    // #pragma omp atomic update
                    pageRanksNext[dest] +=  riDividedOnDiClause[src];

                    // addAtomicFloat(&pageRanksNext[dest] , riDividedOnDiClause[src]);
                }
            }
        }


        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}



struct PageRankStats *pageRankPushColumnFixedPointGraphGrid(struct Arguments *arguments, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Col FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->grid->out_degree[v])
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->grid->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        // pageRankStreamEdgesGraphGridRowWise(graph, riDividedOnDiClause, pageRanksNext);

        uint32_t j;

        #pragma omp parallel for private(j)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions columnwise
        {
            uint32_t i;
            for (i = 0; i < totalPartitions; ++i)
            {
                uint32_t k;
                uint32_t src;
                uint32_t dest;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    src  = partition->edgeList->edges_array_src[k];
                    dest = partition->edgeList->edges_array_dest[k];

                    // #pragma omp atomic update
                    pageRanksNext[dest] +=  riDividedOnDiClause[src];

                    // addAtomicFloat(&pageRanksNext[dest] , riDividedOnDiClause[src]);
                }
            }
        }


        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}




// ********************************************************************************************
// ***************          CSR DataStructure                                    **************
// ********************************************************************************************


struct PageRankStats *pageRankGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    struct PageRankStats *stats = NULL;

    switch (arguments->pushpull)
    {

    case 0: // pull
        stats = pageRankPullGraphCSR(arguments, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphCSR(arguments, graph);
        break;
    case 2: // pull 64bit FP
        stats = pageRankPullFixedPoint64BitGraphCSR(arguments, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphCSR(arguments, graph);
        break;
    case 4: // pull 32bit Quant
        stats = pageRankPullQuant32BitGraphCSR(arguments, graph);
        break;
    case 5: // push
        stats = pageRankPushQuantGraphCSR(arguments, graph);
        break;
    case 6: // pull
        stats = pageRankDataDrivenPullGraphCSR(arguments, graph);
        break;
    case 7: // push
        stats = pageRankDataDrivenPushGraphCSR(arguments, graph);
        break;
    case 8: // pullpush
        stats = pageRankDataDrivenPullPushGraphCSR(arguments, graph);
        break;
    case 9: // pull 32bit FP
        stats = pageRankPullFixedPoint32BitGraphCSR(arguments, graph);
        break;
    case 10: // pull 16bit FP
        stats = pageRankPullFixedPoint16BitGraphCSR(arguments, graph);
        break;
    case 11: // pull 8bit FP
        stats = pageRankPullFixedPoint8BitGraphCSR(arguments, graph);
        break;
    case 12: // pull 16bit Quant
        stats = pageRankPullQuant16BitGraphCSR(arguments, graph);
        break;
    case 13: // pull 8bit Quant
        stats = pageRankPullQuant8BitGraphCSR(arguments, graph);
        break;

    // case 9: // push
    //     pageRankDataDrivenPullFixedPointGraphCSR(arguments, graph);
    // break;
    // case 10: // pull
    //     pageRankDataDrivenPushFixedPointGraphCSR(arguments, graph);
    // break;

    default:// pull
        stats = pageRankPullGraphCSR(arguments, graph);
        break;
    }

    return stats;

}

// topoligy driven approach
struct PageRankStats *pageRankPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }
            pageRanksNext[v] = nodeIncomingPR;
        }


        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{


    double error_total = 0.0;
    // uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) private(v) shared(stats,graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {

            uint32_t degree = graph->vertices->out_degree[v];
            uint32_t edge_idx = graph->vertices->edges_idx[v];
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                uint32_t u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {

            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}


// topoligy driven approach
struct PageRankStats *pageRankPullFixedPoint64BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // uint64_t stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // uint64_t epsilon_fp = DoubleToFixed64(arguments->epsilon);
    // uint64_t num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    // uint64_t* outDegreesFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));
    // uint64_t* pageRanksFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP_64 (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                pageRanksNext[v] += riDividedOnDiClause[u];
            }
        }



        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);
    free(pageRanksNext);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPullFixedPoint32BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // uint64_t stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // uint64_t epsilon_fp = DoubleToFixed64(arguments->epsilon);
    // uint64_t num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint32_t *riDividedOnDiClause = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    // uint64_t* outDegreesFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));
    // uint64_t* pageRanksFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP_32 (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }


    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = FloatToFixed32(stats->pageRanks[v] / graph->vertices->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                pageRanksNext[v] += riDividedOnDiClause[u];
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed32ToFloat(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);
    free(pageRanksNext);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPullFixedPoint16BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // uint64_t stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // uint64_t epsilon_fp = DoubleToFixed64(arguments->epsilon);
    // uint64_t num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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

    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint16_t *riDividedOnDiClause = (uint16_t *) my_malloc(graph->num_vertices * sizeof(uint16_t));
    // uint64_t* outDegreesFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));
    // uint64_t* pageRanksFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP_16 (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = FloatToFixed16(stats->pageRanks[v] / graph->vertices->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                pageRanksNext[v] += riDividedOnDiClause[u];
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed16ToFloat(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);
    free(pageRanksNext);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPullFixedPoint8BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // uint64_t stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // uint64_t epsilon_fp = DoubleToFixed64(arguments->epsilon);
    // uint64_t num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint8_t *riDividedOnDiClause = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    // uint64_t* outDegreesFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));
    // uint64_t* pageRanksFP = (uint64_t*) my_malloc(graph->num_vertices*sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP_8 (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = FloatToFixed8(stats->pageRanks[v] / graph->vertices->out_degree[v]);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                pageRanksNext[v] += riDividedOnDiClause[u];
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed8ToFloat(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);
    free(pageRanksNext);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPushFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    // uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    // uint64_t stats->base_prFP = DoubleToFixed(stats->base_pr);
    // uint64_t stats->dampFP = DoubleToFixed(stats->damp);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    // uint32_t* pageRanksFP = (uint32_t*) my_malloc(graph->num_vertices*sizeof(uint32_t));
    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        // pageRanksFP[v]=stats->base_prFP;
        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
            {
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices->out_degree[v]);
                // riDividedOnDiClause[v] = DIVFixed64V1(pageRanksFP[v],UInt64ToFixed(graph->vertices[v].out_degree));
            }
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) schedule(dynamic, 1024) private(v) shared(stats,graph,pageRanksNext,riDividedOnDiClause) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t degree = graph->vertices->out_degree[v];
            uint32_t edge_idx = graph->vertices->edges_idx[v];
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                uint32_t u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);
                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

//done by mohannad Ibranim
//v_0: No need for next iteration's quantization parameters. (eqn 1)
struct PageRankStats *pageRankPullQuant32BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{
    //QUANT_SCALE = 32;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double error_total = 0.0;

    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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

    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    uint32_t *riDividedOnDiClause_quant = (uint32_t *)my_malloc(graph->num_vertices * sizeof(uint32_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-30s %-19s| \n", "Starting Page Rank Pull Quant_32", "(tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0;
    }


    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        //1. Extract the quantization parameters from riDividedOnDiClause[]
        struct quant_params_32 rDivD_params;
        getMinMax_32(&rDivD_params, riDividedOnDiClause, graph->num_vertices);
        rDivD_params.scale = GetScale_32(rDivD_params.min, rDivD_params.max);
        rDivD_params.zero = 0;
        // printf("Iter %d quant parameters:\nMin = %.16f,\tMax = %.16f\nScale = %.24f,\tZero = %u\n",
        // stats->iterations,rDivD_params.min,rDivD_params.max,rDivD_params.scale,rDivD_params.zero);
        // printf(".........................................................\n");

        //2. Quantize riDividedOnDiClause[]
        #pragma omp parallel for private(v) shared(riDividedOnDiClause_quant,riDividedOnDiClause,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            riDividedOnDiClause_quant[v] = quantize_32(riDividedOnDiClause[v], rDivD_params.scale, rDivD_params.zero);
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint64_t nodeIncomingPR = 0;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                nodeIncomingPR += riDividedOnDiClause_quant[u];
            }
            pageRanksNext[v] = rDivD_params.scale * nodeIncomingPR;
        }

        //uint64_t temp_degree = 0;
        #pragma omp parallel for private(v) shared(arguments,pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + stats->damp * pageRanksNext[v];
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs(nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
                //temp_degree += vertices[v].in_degree;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);
    free(riDividedOnDiClause_quant);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPullQuant16BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{
    //QUANT_SCALE = 32;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double error_total = 0.0;

    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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

    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    uint16_t *riDividedOnDiClause_quant = (uint16_t *)my_malloc(graph->num_vertices * sizeof(uint16_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-30s %-19s| \n", "Starting Page Rank Pull Quant_16", "(tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        //1. Extract the quantization parameters from riDividedOnDiClause[]
        struct quant_params_16 rDivD_params;
        getMinMax_16(&rDivD_params, riDividedOnDiClause, graph->num_vertices);
        rDivD_params.scale = GetScale_16(rDivD_params.min, rDivD_params.max);
        rDivD_params.zero = 0;
        // printf("Iter %d quant parameters:\nMin = %.16f,\tMax = %.16f\nScale = %.24f,\tZero = %u\n",
        // stats->iterations,rDivD_params.min,rDivD_params.max,rDivD_params.scale,rDivD_params.zero);
        // printf(".........................................................\n");

        //2. Quantize riDividedOnDiClause[]
        #pragma omp parallel for private(v) shared(riDividedOnDiClause_quant,riDividedOnDiClause,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            riDividedOnDiClause_quant[v] = quantize_16(riDividedOnDiClause[v], rDivD_params.scale, rDivD_params.zero);
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint64_t nodeIncomingPR = 0;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                nodeIncomingPR += riDividedOnDiClause_quant[u];
            }
            pageRanksNext[v] = rDivD_params.scale * nodeIncomingPR;
        }

        //uint64_t temp_degree = 0;
        #pragma omp parallel for private(v) shared(arguments,pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + stats->damp * pageRanksNext[v];
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs(nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
                //temp_degree += vertices[v].in_degree;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);
    free(riDividedOnDiClause_quant);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPullQuant8BitGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{
    //QUANT_SCALE = 32;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double error_total = 0.0;

    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
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

    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    uint8_t *riDividedOnDiClause_quant = (uint8_t *)my_malloc(graph->num_vertices * sizeof(uint8_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-30s %-19s| \n", "Starting Page Rank Pull Quant_8", "(tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        //1. Extract the quantization parameters from riDividedOnDiClause[]
        struct quant_params_8 rDivD_params;
        getMinMax_8(&rDivD_params, riDividedOnDiClause, graph->num_vertices);
        rDivD_params.scale = GetScale_8(rDivD_params.min, rDivD_params.max);
        rDivD_params.zero = 0;
        // printf("Iter %d quant parameters:\nMin = %.16f,\tMax = %.16f\nScale = %.24f,\tZero = %u\n",
        // stats->iterations,rDivD_params.min,rDivD_params.max,rDivD_params.scale,rDivD_params.zero);
        // printf(".........................................................\n");

        //2. Quantize riDividedOnDiClause[]
        #pragma omp parallel for private(v) shared(riDividedOnDiClause_quant,riDividedOnDiClause,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            riDividedOnDiClause_quant[v] = quantize_8(riDividedOnDiClause[v], rDivD_params.scale, rDivD_params.zero);
        }

        #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint64_t nodeIncomingPR = 0;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = EXTRACT_VALUE(sorted_edges_array[j]);
                nodeIncomingPR += riDividedOnDiClause_quant[u];
            }
            //nodeIncomingPR -= (degree * rDivD_params.zero);
            pageRanksNext[v] = rDivD_params.scale * nodeIncomingPR;
        }

        //uint64_t temp_degree = 0;
        #pragma omp parallel for private(v) shared(arguments,pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + stats->damp * pageRanksNext[v];
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs(nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
                //temp_degree += vertices[v].in_degree;
            }
        }

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);
    free(riDividedOnDiClause_quant);

    stats->error_total = error_total;
    return stats;
}

//done by mohannad Ibranim
struct PageRankStats *pageRankPushQuantGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{
    //QUANT_SCALE = 16;
    // uint32_t i;
    uint32_t v;
    uint32_t activeVertices = 0;
    double error_total = 0.0;

    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    uint32_t *riDividedOnDiClause_quant = (uint32_t *)my_malloc(graph->num_vertices * sizeof(uint32_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-30s %-19s| \n", "Starting Page Rank Push Quant_32", "(tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0;
    }


    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        //1. Extract the quantization parameters from riDividedOnDiClause[]
        struct quant_params rDivD_params;
        getMinMax(&rDivD_params, riDividedOnDiClause, graph->num_vertices);
        rDivD_params.scale = GetScale(rDivD_params.min, rDivD_params.max);
        rDivD_params.zero = GetZeroPoint(rDivD_params.max, rDivD_params.scale);
        // printf("Itaration %d's quant parameters:\n\tMin = %f,\tMax = %f\nScale = %f,\tZero = %d\n......................",
        // stats->iterations,rDivD_params.min,rDivD_params.max,rDivD_params.scale,rDivD_params.zero);

        //2. Quantize riDividedOnDiClause[]
        #pragma omp parallel for private(v) shared(riDividedOnDiClause_quant,riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {
            riDividedOnDiClause_quant[v] = quantize(riDividedOnDiClause[v], rDivD_params.scale, rDivD_params.zero);
        }

        #pragma omp parallel for default(none) private(v) shared(stats,rDivD_params,riDividedOnDiClause_quant,graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {

            uint32_t degree = graph->vertices->out_degree[v];
            uint32_t edge_idx = graph->vertices->edges_idx[v];
            uint32_t j;

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                uint32_t u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);

                #pragma omp atomic update
                pageRanksNext[u] += rDivD_params.scale * (riDividedOnDiClause_quant[v] - rDivD_params.zero);
            }
        }

        #pragma omp parallel for private(v) shared(arguments, stats,pageRanksNext) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {

            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr +  stats->damp * pageRanksNext[v];
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0;
            double error = fabs(nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;
    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, Seconds(timer));
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}



struct PageRankStats *pageRankDataDrivenPullGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{



    double error_total = 0.0;
    uint32_t i;
    uint32_t v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif



    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);
    #pragma omp parallel for reduction(+:activeVertices)
    for(i = 0; i < graph->num_vertices; i++)
    {
        workListNext[i] = 1;
        activeVertices++;
    }

    swapWorkLists(&workListNext, &workListCurr);
    resetWorkList(workListNext, graph->num_vertices);
    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for default(none) shared(arguments,riDividedOnDiClause,sorted_edges_array,vertices,workListCurr,workListNext,stats,graph) private(v) reduction(+:activeVertices,error_total) schedule(dynamic, 1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                uint32_t edge_idx;
                uint32_t degree;
                uint32_t j;
                uint32_t u;
                double error = 0;
                float nodeIncomingPR = 0;
                degree = vertices->out_degree[v]; // when directed we use inverse graph out degree means in degree
                edge_idx = vertices->edges_idx[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = EXTRACT_VALUE(sorted_edges_array[j]);
                    nodeIncomingPR += riDividedOnDiClause[u]; // sum (PRi/outDegree(i))
                }
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  stats->base_pr + (stats->damp * nodeIncomingPR);
                error = fabs(newPageRank - oldPageRank);
                error_total += error / graph->num_vertices;
                if(error >= arguments->epsilon)
                {
                    stats->pageRanks[v] = newPageRank;
                    degree = graph->vertices->out_degree[v];
                    edge_idx = graph->vertices->edges_idx[v];
                    for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                    {
                        u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);

                        #pragma omp atomic write
                        workListNext[u] = 1;

                        // uint8_t old_val = workListNext[u];
                        // if(!old_val){
                        //    __sync_bool_compare_and_swap(&workListNext[u], 0, 1);
                        // }
                    }
                    activeVertices++;
                }
            }
        }
        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankDataDrivenPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t edge_idx;
    uint32_t degree;
    uint32_t j;
    uint32_t u;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;

    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));

    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif

    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);


    #pragma omp parallel for private(edge_idx,degree,v,j,u) shared(workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;
        degree = vertices->out_degree[v]; // when directed we use inverse graph out degree means in degree
        edge_idx = vertices->edges_idx[v];
        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = EXTRACT_VALUE(sorted_edges_array[j]);
            if(graph->vertices->out_degree[u])
                aResiduals[v] += 1.0f / graph->vertices->out_degree[u]; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(edge_idx,degree,v,j,u) shared(stats,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  aResiduals[v] + stats->pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                // #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                degree = graph->vertices->out_degree[v];
                float delta = stats->damp * (aResiduals[v] / degree);
                edge_idx = graph->vertices->edges_idx[v];

                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];
                    #pragma omp atomic update
                    aResiduals[u] += delta;
                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        if(!workListNext[u])
                        {
                            // #pragma omp atomic write
                            workListNext[u] = 1;
                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}


struct PageRankStats *pageRankDataDrivenPullPushGraphCSR(struct Arguments *arguments, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t edge_idx;
    uint32_t degree;
    uint32_t j;
    uint32_t u;



    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif



    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);


    #pragma omp parallel for private(edge_idx,degree,v,j,u) shared(workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0f;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;
        degree = vertices->out_degree[v]; // when directed we use inverse graph out degree means in degree
        edge_idx = vertices->edges_idx[v];
        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = EXTRACT_VALUE(sorted_edges_array[j]);
            if(graph->vertices->out_degree[u])
                aResiduals[v] += 1.0f / graph->vertices->out_degree[u]; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(edge_idx,degree,v,j,u) shared(stats,vertices,sorted_edges_array,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024) num_threads(arguments->ker_numThreads)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                float nodeIncomingPR = 0.0f;
                degree = vertices->out_degree[v];
                edge_idx = vertices->edges_idx[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = EXTRACT_VALUE(sorted_edges_array[j]);
                    nodeIncomingPR += stats->pageRanks[u] / graph->vertices->out_degree[u];
                }

                float newPageRank = stats->base_pr + (stats->damp * nodeIncomingPR);
                float oldPageRank =  stats->pageRanks[v];
                // float newPageRank =  aResiduals[v]+pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                degree = graph->vertices->out_degree[v];
                float delta = stats->damp * (aResiduals[v] / degree);
                edge_idx = graph->vertices->edges_idx[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = EXTRACT_VALUE(graph->sorted_edges_array->edges_array_dest[j]);
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        aResiduals[u] += delta;
                        if(!workListNext[u])
                        {
                            workListNext[u] = 1;
                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}


// float* pageRankDataDrivenPullFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPushFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPullPushFixedPointGraphCSR(struct Arguments *arguments, struct GraphCSR* graph){


// }



// ********************************************************************************************
// ***************          ArrayList DataStructure              **************
// ********************************************************************************************


struct PageRankStats *pageRankGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    struct PageRankStats *stats = NULL;

    switch (arguments->pushpull)
    {

    case 0: // pull
        stats = pageRankPullGraphAdjArrayList(arguments, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphAdjArrayList(arguments, graph);
        break;
    case 2: // pull
        stats = pageRankPullFixedPointGraphAdjArrayList(arguments, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphAdjArrayList(arguments, graph);
        break;
    case 4: // pull
        stats = pageRankDataDrivenPullGraphAdjArrayList(arguments, graph);
        break;
    case 5: // push
        stats = pageRankDataDrivenPushGraphAdjArrayList(arguments, graph);
        break;
    case 6: // pullpush
        stats = pageRankDataDrivenPullPushGraphAdjArrayList(arguments, graph);
        break;
    default:// push
        stats = pageRankPullGraphAdjArrayList(arguments, graph);
        break;
    }


    return stats;

}

struct PageRankStats *pageRankPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t activeVertices = 0;

    struct EdgeList *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;

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
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }

            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    struct EdgeList *Nodes;


    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graph->num_vertices * sizeof(omp_lock_t));




    #pragma omp parallel for default(none) private(i) shared(graph,vertex_lock)
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }



    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {


            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) private(v,Nodes) shared(graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {

            Nodes = graph->vertices[v].outNodes;
            uint32_t degree = graph->vertices[v].out_degree;
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = 0 ; j < (degree) ; j++)
            {
                uint32_t u = Nodes->edges_array_dest[j];

                // omp_set_lock(&(vertex_lock[u]));
                //   pageRanksNext[u] += riDividedOnDiClause[v];
                // omp_unset_lock((&vertex_lock[u]));

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

                // __atomic_fetch_add(&pageRanksNext[u], riDividedOnDiClause[v], __ATOMIC_RELAXED);
                // printf("tid %u degree %u edge_idx %u v %u u %u \n",tid,degree,edge_idx,v,u );

                // addAtomicFloat(&pageRanksNext[u] , riDividedOnDiClause[v]);
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    #pragma omp parallel for
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_destroy_lock(&(vertex_lock[i]));
    }

    free(timer);
    free(timer_inner);
    free(vertex_lock);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPullFixedPointGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t activeVertices = 0;

    struct EdgeList *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices[v].out_degree);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;

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
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }

            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPushFixedPointGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    struct EdgeList *Nodes;


    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graph->num_vertices * sizeof(omp_lock_t));



    #pragma omp parallel for default(none) private(i) shared(graph,vertex_lock)
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {


            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices[v].out_degree);
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) private(v,Nodes) shared(graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {

            Nodes = graph->vertices[v].outNodes;
            uint32_t degree = graph->vertices[v].out_degree;
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = 0 ; j < (degree) ; j++)
            {
                uint32_t u = Nodes->edges_array_dest[j];

                // omp_set_lock(&(vertex_lock[u]));
                //   pageRanksNext[u] += riDividedOnDiClause[v];
                // omp_unset_lock((&vertex_lock[u]));

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

                // __atomic_fetch_add(&pageRanksNext[u], riDividedOnDiClause[v], __ATOMIC_RELAXED);
                // printf("tid %u degree %u edge_idx %u v %u u %u \n",tid,degree,edge_idx,v,u );

                // addAtomicFloat(&pageRanksNext[u] , riDividedOnDiClause[v]);
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    #pragma omp parallel for
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_destroy_lock(&(vertex_lock[i]));
    }

    free(timer);
    free(timer_inner);
    free(vertex_lock);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankDataDrivenPullGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;
    struct EdgeList *Nodes;


    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);
    #pragma omp parallel for reduction(+:activeVertices)
    for(i = 0; i < graph->num_vertices; i++)
    {
        workListNext[i] = 1;
        activeVertices++;
    }

    swapWorkLists(&workListNext, &workListCurr);
    resetWorkList(workListNext, graph->num_vertices);
    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for default(none) shared(arguments,riDividedOnDiClause,workListCurr,workListNext,stats,graph) private(v,Nodes) reduction(+:activeVertices,error_total) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                uint32_t degree;
                uint32_t j;
                uint32_t u;
                double error = 0;
                float nodeIncomingPR = 0;

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
                    nodeIncomingPR += riDividedOnDiClause[u]; // sum (PRi/outDegree(i))
                }
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  stats->base_pr + (stats->damp * nodeIncomingPR);
                error = fabs(newPageRank - oldPageRank);
                error_total += error / graph->num_vertices;
                if(error >= arguments->epsilon)
                {
                    stats->pageRanks[v] = newPageRank;
                    Nodes = graph->vertices[v].outNodes;
                    degree = graph->vertices[v].out_degree;
                    for(j = 0 ; j < (degree) ; j++)
                    {
                        u = Nodes->edges_array_dest[j];

                        #pragma omp atomic write
                        workListNext[u] = 1;
                        // uint8_t old_val = workListNext[u];
                        // if(!old_val){
                        //    __sync_bool_compare_and_swap(&workListNext[u], 0, 1);
                        // }
                    }
                    activeVertices++;
                }
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankDataDrivenPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t degree;
    uint32_t j;
    uint32_t u;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;
    struct EdgeList *Nodes;

    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);

    #pragma omp parallel for private(Nodes,degree,v,j,u) shared(workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;


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
            if(graph->vertices[u].out_degree)
                aResiduals[v] += 1.0f / graph->vertices[u].out_degree; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  aResiduals[v] + stats->pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                // #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                Nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;
                float delta = stats->damp * (aResiduals[v] / degree);

                for(j = 0 ; j < (degree) ; j++)
                {
                    u = Nodes->edges_array_dest[j];
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        if(!workListNext[u])
                        {

                            // #pragma omp atomic write
                            workListNext[u] = 1;

                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjArrayList(struct Arguments *arguments, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t degree;
    uint32_t j;
    uint32_t u;
    struct EdgeList *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;

    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);


    #pragma omp parallel for private(Nodes,degree,v,j,u) shared(workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0f;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;


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
            if(graph->vertices[u].out_degree)
                aResiduals[v] += 1.0f / graph->vertices[u].out_degree; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                float nodeIncomingPR = 0.0f;

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
                    nodeIncomingPR += stats->pageRanks[u] / graph->vertices[u].out_degree;
                }

                float newPageRank = stats->base_pr + (stats->damp * nodeIncomingPR);
                float oldPageRank =  stats->pageRanks[v];
                // float newPageRank =  aResiduals[v]+pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                Nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;
                float delta = stats->damp * (aResiduals[v] / degree);



                for(j = 0 ; j < (degree) ; j++)
                {
                    uint32_t u = Nodes->edges_array_dest[j];
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        if(!workListNext[u])
                        {
                            workListNext[u] = 1;
                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}


// ********************************************************************************************
// ***************          LinkedList DataStructure           **************
// ********************************************************************************************


struct PageRankStats *pageRankGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    struct PageRankStats *stats = NULL;

    switch (arguments->pushpull)
    {
    case 0: // pull
        stats = pageRankPullGraphAdjLinkedList(arguments, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphAdjLinkedList(arguments, graph);
        break;
    case 2: // pull
        stats = pageRankPullFixedPointGraphAdjLinkedList(arguments, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphAdjLinkedList(arguments, graph);
        break;
    case 4: // pull
        stats = pageRankDataDrivenPullGraphAdjLinkedList(arguments, graph);
        break;
    case 5: // push
        stats = pageRankDataDrivenPushGraphAdjLinkedList(arguments, graph);
        break;
    case 6: // pullpush
        stats = pageRankDataDrivenPullPushGraphAdjLinkedList(arguments, graph);
        break;
    default:// push
        stats = pageRankPullGraphAdjLinkedList(arguments, graph);
        break;
    }


    return stats;

}

struct PageRankStats *pageRankPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t activeVertices = 0;

    struct AdjLinkedListNode *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;

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
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }

            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    struct AdjLinkedListNode *Nodes;


    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graph->num_vertices * sizeof(omp_lock_t));




    #pragma omp parallel for default(none) private(i) shared(graph,vertex_lock)
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }



    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {


            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) private(v,Nodes) shared(graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {

            Nodes = graph->vertices[v].outNodes;
            uint32_t degree = graph->vertices[v].out_degree;
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = 0 ; j < (degree) ; j++)
            {
                uint32_t u = Nodes->dest;
                Nodes = Nodes->next;

                // omp_set_lock(&(vertex_lock[u]));
                //   pageRanksNext[u] += riDividedOnDiClause[v];
                // omp_unset_lock((&vertex_lock[u]));

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

                // __atomic_fetch_add(&pageRanksNext[u], riDividedOnDiClause[v], __ATOMIC_RELAXED);
                // printf("tid %u degree %u edge_idx %u v %u u %u \n",tid,degree,edge_idx,v,u );

                // addAtomicFloat(&pageRanksNext[u] , riDividedOnDiClause[v]);
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    #pragma omp parallel for
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_destroy_lock(&(vertex_lock[i]));
    }

    free(timer);
    free(timer_inner);
    free(vertex_lock);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPullFixedPointGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t activeVertices = 0;

    struct AdjLinkedListNode *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        error_total = 0;
        activeVertices = 0;
        Start(timer_inner);
        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices[v].out_degree);
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;

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
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }

            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }


        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // printf(" -----------------------------------------------------\n");
    // printf("| %-10s | %-8lf | %-15s | %-9s | \n","PR Sum ",sum, stats->iterations, stats->time_total);
    // printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPushFixedPointGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    struct AdjLinkedListNode *Nodes;


    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graph->num_vertices * sizeof(omp_lock_t));




    #pragma omp parallel for default(none) private(i) shared(graph,vertex_lock)
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;
        #pragma omp parallel for private(v) shared(riDividedOnDiClause,stats,graph)
        for(v = 0; v < graph->num_vertices; v++)
        {


            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = DoubleToFixed64(stats->pageRanks[v] / graph->vertices[v].out_degree);
            else
                riDividedOnDiClause[v] = 0.0f;

        }

        #pragma omp parallel for default(none) private(v,Nodes) shared(graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {

            Nodes = graph->vertices[v].outNodes;
            uint32_t degree = graph->vertices[v].out_degree;
            // uint32_t tid = omp_get_thread_num();
            uint32_t j;

            for(j = 0 ; j < (degree) ; j++)
            {
                uint32_t  u = Nodes->dest;
                Nodes = Nodes->next;
                // omp_set_lock(&(vertex_lock[u]));
                //   pageRanksNext[u] += riDividedOnDiClause[v];
                // omp_unset_lock((&vertex_lock[u]));

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

                // __atomic_fetch_add(&pageRanksNext[u], riDividedOnDiClause[v], __ATOMIC_RELAXED);
                // printf("tid %u degree %u edge_idx %u v %u u %u \n",tid,degree,edge_idx,v,u );

                // addAtomicFloat(&pageRanksNext[u] , riDividedOnDiClause[v]);
            }
        }

        #pragma omp parallel for private(v) shared(arguments, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= arguments->epsilon)
            {
                activeVertices++;
            }
        }

        Stop(timer_inner);

        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");
    // pageRankPrint(pageRanks, graph->num_vertices);

    #pragma omp parallel for
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_destroy_lock(&(vertex_lock[i]));
    }

    free(timer);
    free(timer_inner);
    free(vertex_lock);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankDataDrivenPullGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t i;
    uint32_t v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;
    struct AdjLinkedListNode *Nodes;


    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);





    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);
    #pragma omp parallel for reduction(+:activeVertices)
    for(i = 0; i < graph->num_vertices; i++)
    {
        workListNext[i] = 1;
        activeVertices++;
    }

    swapWorkLists(&workListNext, &workListCurr);
    resetWorkList(workListNext, graph->num_vertices);
    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(graph->vertices[v].out_degree)
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices[v].out_degree;
            else
                riDividedOnDiClause[v] = 0.0f;
        }

        #pragma omp parallel for default(none) shared(arguments, riDividedOnDiClause, workListCurr, workListNext, stats, graph) private(v,Nodes) reduction(+:activeVertices,error_total) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                uint32_t degree;
                uint32_t j;
                uint32_t u;
                double error = 0;
                float nodeIncomingPR = 0;

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
                    nodeIncomingPR += riDividedOnDiClause[u]; // sum (PRi/outDegree(i))
                }
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  stats->base_pr + (stats->damp * nodeIncomingPR);
                error = fabs(newPageRank - oldPageRank);
                error_total += error / graph->num_vertices;
                if(error >= arguments->epsilon)
                {
                    stats->pageRanks[v] = newPageRank;
                    Nodes = graph->vertices[v].outNodes;
                    degree = graph->vertices[v].out_degree;
                    for(j = 0 ; j < (degree) ; j++)
                    {
                        u = Nodes->dest;
                        Nodes = Nodes->next;
                        #pragma omp atomic write
                        workListNext[u] = 1;
                        // uint8_t old_val = workListNext[u];
                        // if(!old_val){
                        //    __sync_bool_compare_and_swap(&workListNext[u], 0, 1);
                        // }
                    }
                    activeVertices++;
                }
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop

    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankDataDrivenPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t degree;
    uint32_t j;
    uint32_t u;



    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;
    struct AdjLinkedListNode *Nodes;



    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));



    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);


    #pragma omp parallel for private(Nodes,degree,v,j,u) shared(workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;


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
            if(graph->vertices[u].out_degree)
                aResiduals[v] += 1.0f / graph->vertices[u].out_degree; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  aResiduals[v] + stats->pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                // #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                Nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;
                float delta = stats->damp * (aResiduals[v] / degree);


                for(j = 0 ; j < (degree) ; j++)
                {
                    u = Nodes->dest;
                    Nodes = Nodes->next;
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        if(!workListNext[u])
                        {

                            // #pragma omp atomic write
                            workListNext[u] = 1;

                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");


    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjLinkedList(struct Arguments *arguments, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t degree;
    uint32_t j;
    uint32_t u;
    struct AdjLinkedListNode *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    uint8_t *workListCurr = NULL;
    uint8_t *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));
    workListNext  = (uint8_t *) my_malloc(graph->num_vertices * sizeof(uint8_t));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", arguments->epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    Start(timer_inner);


    #pragma omp parallel for private(Nodes,degree,v,j,u) shared(stats,workListCurr,workListNext,aResiduals) reduction(+:activeVertices)
    for(v = 0; v < graph->num_vertices; v++)
    {

        aResiduals[v] = 0.0f;
        workListCurr[v] = 1;
        workListNext[v] = 0;
        activeVertices++;


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
            if(graph->vertices[u].out_degree)
                aResiduals[v] += 1.0f / graph->vertices[u].out_degree; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < arguments->iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,arguments,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                float nodeIncomingPR = 0.0f;

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
                    nodeIncomingPR += stats->pageRanks[u] / graph->vertices[u].out_degree;
                }

                float newPageRank = stats->base_pr + (stats->damp * nodeIncomingPR);
                float oldPageRank =  stats->pageRanks[v];
                // float newPageRank =  aResiduals[v]+pageRanks[v];
                error_total += fabs(newPageRank / graph->num_vertices - oldPageRank / graph->num_vertices);

                #pragma omp atomic write
                stats->pageRanks[v] = newPageRank;

                Nodes = graph->vertices[v].outNodes;
                degree = graph->vertices[v].out_degree;

                float delta = stats->damp * (aResiduals[v] / degree);



                for(j = 0 ; j < (degree) ; j++)
                {
                    u = Nodes->dest;
                    Nodes = Nodes->next;
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= arguments->epsilon) && (prevResidual <= arguments->epsilon))
                    {
                        activeVertices++;
                        if(!workListNext[u])
                        {
                            workListNext[u] = 1;
                        }
                    }
                }
                aResiduals[v] = 0.0f;
            }
        }

        // activeVertices = getNumOfSetBits(workListNext);
        swapWorkLists(&workListNext, &workListCurr);
        resetWorkList(workListNext, graph->num_vertices);

        Stop(timer_inner);
        printf("| %-10u | %-8u | %-15.13lf | %-9f | \n", stats->iterations, activeVertices, error_total, Seconds(timer_inner));
        if(activeVertices == 0)
            break;

    }// end iteration loop


    double sum = 0.0f;
    #pragma omp parallel for reduction(+:sum)
    for(v = 0; v < graph->num_vertices; v++)
    {
        stats->pageRanks[v] = stats->pageRanks[v] / graph->num_vertices;
        sum += stats->pageRanks[v];
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iterations", "PR Sum", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-10u | %-8lf | %-15.13lf | %-9f | \n", stats->iterations, sum, error_total, stats->time_total);
    printf(" -----------------------------------------------------\n");

    // pageRankPrint(pageRanks, graph->num_vertices);
    free(workListCurr);
    free(workListNext);
    free(timer);
    free(timer_inner);
    free(aResiduals);
    free(riDividedOnDiClause);


    stats->error_total = error_total;
    return stats;

}
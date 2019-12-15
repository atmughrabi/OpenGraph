// -----------------------------------------------------------------------------
//
//      "OpenGraph"
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

#include "graphConfig.h"

#include "fixedPoint.h"

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

    __u32 v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));;
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

    __u32 v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));;
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

    __u32 v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));;
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

    __u32 v;

    struct PageRankStats *stats = (struct PageRankStats *) my_malloc(sizeof(struct PageRankStats));

    stats->damp = Damp;
    stats->base_pr = (1.0f - stats->damp);
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0;
    stats->error_total = 0.0;

    stats->realRanks = (__u32 *) my_malloc(graph->num_vertices * sizeof(__u32));;
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
    __u32 *lnewV;
    __u32 *loldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
        loldV = (__u32 *)&oldV;
        lnewV = (__u32 *)&newV;
    }
    while(!__sync_bool_compare_and_swap((__u32 *)num, *(loldV), *(lnewV)));

}


void addAtomicDouble(double *num, double value)
{

    double newV, oldV;
    __u64 *lnewV;
    __u64 *loldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
        loldV = (__u64 *)&oldV;
        lnewV = (__u64 *)&newV;
    }
    while(!__sync_bool_compare_and_swap((__u64 *)num, *(loldV), *(lnewV)));

}


void swapWorkLists (__u8 **workList1, __u8 **workList2)
{


    __u8 *workList_temp = *workList1;
    *workList1 = *workList2;
    *workList2 = workList_temp;

}

void resetWorkList(__u8 *workList, __u32 size)
{

    __u32 i;

    #pragma omp parallel for
    for(i = 0; i < size ; i++)
    {
        workList[i] = 0;

    }


}

void setWorkList(__u8 *workList,  __u32 size)
{

    __u32 i;

    #pragma omp parallel for
    for(i = 0; i < size ; i++)
    {
        workList[i] = 1;

    }


}

void setAtomic(__u64 *num, __u64 value)
{


    __u64 newV, oldV;

    do
    {
        oldV = *num;
        newV = value;
    }
    while(!__sync_bool_compare_and_swap(num, oldV, newV));

}

void addAtomicFixedPoint(__u64 *num, __u64 value)
{

    __u64 newV, oldV;

    do
    {
        oldV = *num;
        newV = oldV + value;
    }
    while(!__sync_bool_compare_and_swap(num, oldV, newV));

}

void pageRankPrint(float *pageRankArray, __u32 num_vertices)
{
    __u32 v;
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

struct PageRankStats  *pageRankGraphGrid(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphGrid *graph)
{

    struct PageRankStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // push
        stats = pageRankPullRowGraphGrid(epsilon, iterations, graph);
        break;
    case 1: // pull
        stats = pageRankPushColumnGraphGrid(epsilon, iterations, graph);
        break;
    case 2: // pull
        stats = pageRankPullRowFixedPointGraphGrid(epsilon, iterations, graph);
        break;
    case 3: // push
        stats = pageRankPushColumnFixedPointGraphGrid(epsilon, iterations, graph);
        break;
    default:// pull
        stats = pageRankPullRowGraphGrid(epsilon, iterations, graph);
        break;
    }

    return stats;

}


struct PageRankStats *pageRankPullRowGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;

    __u32 totalPartitions  = graph->grid->num_partitions;

    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Row (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {
        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        __u32 i;
        // #pragma omp parallel for private(i)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            __u32 j;
            #pragma omp parallel for private(j)
            for (j = 0; j < totalPartitions; ++j)
            {
                __u32 k;
                __u32 src;
                __u32 dest;
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


        #pragma omp parallel for private(v) shared(epsilon,pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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



struct PageRankStats *pageRankPullRowFixedPointGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    __u32 totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Row FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        __u32 i;
        // #pragma omp parallel for private(i)
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            __u32 j;
            #pragma omp parallel for private(j)
            for (j = 0; j < totalPartitions; ++j)
            {
                __u32 k;
                __u32 src;
                __u32 dest;
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


        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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



struct PageRankStats *pageRankPushColumnGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    __u32 totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Col (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        __u32 j;
        #pragma omp parallel for private(j)
        for (j = 0; j < totalPartitions; ++j)
        {
            __u32 i;

            // #pragma omp parallel for private(i) // iterate over partitions columnwise
            for (i = 0; i < totalPartitions; ++i)
            {
                __u32 k;
                __u32 src;
                __u32 dest;
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


        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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



struct PageRankStats *pageRankPushColumnFixedPointGraphGrid(double epsilon,  __u32 iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    __u32 totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));



    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Col FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        __u32 j;

        #pragma omp parallel for private(j)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions columnwise
        {
            __u32 i;
            for (i = 0; i < totalPartitions; ++i)
            {
                __u32 k;
                __u32 src;
                __u32 dest;
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


        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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


struct PageRankStats *pageRankGraphCSR(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphCSR *graph)
{

    struct PageRankStats *stats = NULL;

    switch (pushpull)
    {

    case 0: // pull
        stats = pageRankPullGraphCSR(epsilon, iterations, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphCSR(epsilon, iterations, graph);
        break;

    case 2: // pull
        stats = pageRankPullFixedPointGraphCSR(epsilon, iterations, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphCSR(epsilon, iterations, graph);
        break;
    case 4: // pull
        stats = pageRankDataDrivenPullGraphCSR(epsilon, iterations, graph);
        break;
    case 5: // push
        stats = pageRankDataDrivenPushGraphCSR(epsilon, iterations, graph);
        break;
    case 6: // pullpush
        stats = pageRankDataDrivenPullPushGraphCSR(epsilon, iterations, graph);
        break;

    // case 7: // push
    //     pageRankDataDrivenPullFixedPointGraphCSR(epsilon, iterations, graph);
    // break;
    // case 8: // pull
    //     pageRankDataDrivenPushFixedPointGraphCSR(epsilon, iterations, graph);
    // break;

    default:// pull
        stats = pageRankPullGraphCSR(epsilon, iterations, graph);
        break;
    }

    return stats;

}

// topoligy driven approach
struct PageRankStats *pageRankPullGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 edge_idx;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    __u32 *sorted_edges_array = NULL;
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


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float nodeIncomingPR = 0.0f;
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = sorted_edges_array[j];
                nodeIncomingPR += riDividedOnDiClause[u]; // stats->pageRanks[v]/graph->vertices[v].out_degree;
            }
            pageRanksNext[v] = nodeIncomingPR;
        }

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{


    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


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
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext,riDividedOnDiClause) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {

            __u32 degree = graph->vertices->out_degree[v];
            __u32 edge_idx = graph->vertices->edges_idx[v];
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                __u32 u = graph->sorted_edges_array->edges_array_dest[j];

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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {

            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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


// topoligy driven approach
struct PageRankStats *pageRankPullFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 edge_idx;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // __u64 stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // __u64 epsilon_fp = DoubleToFixed64(epsilon);
    // __u64 num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    __u32 *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

#if DIRECTED
    vertices = graph->inverse_vertices;
    sorted_edges_array = graph->inverse_sorted_edges_array->edges_array_dest;
#else
    vertices = graph->vertices;
    sorted_edges_array = graph->sorted_edges_array->edges_array_dest;
#endif



    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    // __u64* outDegreesFP = (__u64*) my_malloc(graph->num_vertices*sizeof(__u64));
    // __u64* pageRanksFP = (__u64*) my_malloc(graph->num_vertices*sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        // Stop(timer_inner);
        // printf("|A %-9u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));

        //  Start(timer_inner);
        #pragma omp parallel for reduction(+ : error_total,activeVertices) private(v,j,u,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            degree = vertices->out_degree[v];
            edge_idx = vertices->edges_idx[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = sorted_edges_array[j];
                pageRanksNext[v] += riDividedOnDiClause[u];
            }

        }
        // Stop(timer_inner);
        // printf("|B %-9u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));

        // Start(timer_inner);
        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
            {
                activeVertices++;
            }
        }

        // Stop(timer_inner);
        // printf("|C %-9u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));


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

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    // __u64 stats->base_prFP = DoubleToFixed(stats->base_pr);
    // __u64 stats->dampFP = DoubleToFixed(stats->damp);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graph->num_vertices * sizeof(omp_lock_t));



    #pragma omp parallel for default(none) private(i) shared(graph,vertex_lock)
    for (i = 0; i < graph->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }



    // __u32* pageRanksFP = (__u32*) my_malloc(graph->num_vertices*sizeof(__u32));
    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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
        // Stop(timer_inner);
        // printf("|A%-10u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));
        // Start(timer_inner);
        #pragma omp parallel for default(none) schedule(dynamic, 1024) private(v) shared(graph,pageRanksNext,riDividedOnDiClause)
        for(v = 0; v < graph->num_vertices; v++)
        {



            __u32 degree = graph->vertices->out_degree[v];
            __u32 edge_idx = graph->vertices->edges_idx[v];
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                __u32 u = graph->sorted_edges_array->edges_array_dest[j];

                // omp_set_lock(&(vertex_lock[u]));
                //   pageRanksNext[u] += riDividedOnDiClause[v];
                // omp_unset_lock((&vertex_lock[u]));

                #pragma omp atomic update
                pageRanksNext[u] += riDividedOnDiClause[v];

                // addAtomicFixedPoint(&pageRanksNext[u] , riDividedOnDiClause[v]);
            }
        }
        // Stop(timer_inner);
        // printf("|B%-10u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));
        // Start(timer_inner);
        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            // pageRanksFP[v] = FloatToFixed(nextPageRank);
            pageRanksNext[v] = 0;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankDataDrivenPullGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{



    double error_total = 0.0;
    __u32 i;
    __u32 v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    __u32 *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


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
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for default(none) shared(epsilon,riDividedOnDiClause,sorted_edges_array,vertices,workListCurr,workListNext,stats,graph) private(v) reduction(+:activeVertices,error_total) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {
                __u32 edge_idx;
                __u32 degree;
                __u32 j;
                __u32 u;
                double error = 0;
                float nodeIncomingPR = 0;
                degree = vertices->out_degree[v]; // when directed we use inverse graph out degree means in degree
                edge_idx = vertices->edges_idx[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = sorted_edges_array[j];
                    nodeIncomingPR += riDividedOnDiClause[u]; // sum (PRi/outDegree(i))
                }
                float oldPageRank =  stats->pageRanks[v];
                float newPageRank =  stats->base_pr + (stats->damp * nodeIncomingPR);
                error = fabs(newPageRank - oldPageRank);
                error_total += error / graph->num_vertices;
                if(error >= epsilon)
                {
                    stats->pageRanks[v] = newPageRank;
                    degree = graph->vertices->out_degree[v];
                    edge_idx = graph->vertices->edges_idx[v];
                    for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                    {
                        u = graph->sorted_edges_array->edges_array_dest[j];

                        #pragma omp atomic write
                        workListNext[u] = 1;
                        // __u8 old_val = workListNext[u];
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

struct PageRankStats *pageRankDataDrivenPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 edge_idx;
    __u32 degree;
    __u32 j;
    __u32 u;



    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    __u32 *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


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
    printf("| %-51.13lf | \n", epsilon);
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
            u = sorted_edges_array[j];
            if(graph->vertices->out_degree[u])
                aResiduals[v] += 1.0f / graph->vertices->out_degree[u]; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(edge_idx,degree,v,j,u) shared(stats,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
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
                    u = graph->sorted_edges_array->edges_array_dest[j];
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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


struct PageRankStats *pageRankDataDrivenPullPushGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 edge_idx;
    __u32 degree;
    __u32 j;
    __u32 u;



    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    __u32 *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


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
    printf("| %-51.13lf | \n", epsilon);
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
            u = sorted_edges_array[j];
            if(graph->vertices->out_degree[u])
                aResiduals[v] += 1.0f / graph->vertices->out_degree[u]; // sum (PRi/outDegree(i))
        }
        aResiduals[v] = (1.0f - stats->damp) * stats->damp * aResiduals[v];
    }

    Stop(timer_inner);
    printf("| %-10s | %-8u | %-15.13lf | %-9f | \n", "Init", activeVertices, error_total, Seconds(timer_inner));

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(edge_idx,degree,v,j,u) shared(stats,vertices,sorted_edges_array,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                float nodeIncomingPR = 0.0f;
                degree = vertices->out_degree[v];
                edge_idx = vertices->edges_idx[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = sorted_edges_array[j];
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
                    u = graph->sorted_edges_array->edges_array_dest[j];
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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


// float* pageRankDataDrivenPullFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPullPushFixedPointGraphCSR(double epsilon,  __u32 iterations, struct GraphCSR* graph){


// }



// ********************************************************************************************
// ***************          ArrayList DataStructure              **************
// ********************************************************************************************


struct PageRankStats *pageRankGraphAdjArrayList(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphAdjArrayList *graph)
{

    struct PageRankStats *stats = NULL;

    switch (pushpull)
    {

    case 0: // pull
        stats = pageRankPullGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 2: // pull
        stats = pageRankPullFixedPointGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 4: // pull
        stats = pageRankDataDrivenPullGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 5: // push
        stats = pageRankDataDrivenPushGraphAdjArrayList(epsilon, iterations, graph);
        break;
    case 6: // pullpush
        stats = pageRankDataDrivenPullPushGraphAdjArrayList(epsilon, iterations, graph);
        break;
    default:// push
        stats = pageRankPullGraphAdjArrayList(epsilon, iterations, graph);
        break;
    }


    return stats;

}

struct PageRankStats *pageRankPullGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 activeVertices = 0;

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
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

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
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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
            __u32 degree = graph->vertices[v].out_degree;
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = 0 ; j < (degree) ; j++)
            {
                __u32 u = Nodes->edges_array_dest[j];

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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPullFixedPointGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 activeVertices = 0;

    struct EdgeList *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPushFixedPointGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

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



    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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
            __u32 degree = graph->vertices[v].out_degree;
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = 0 ; j < (degree) ; j++)
            {
                __u32 u = Nodes->edges_array_dest[j];

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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankDataDrivenPullGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;
    struct EdgeList *Nodes;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for default(none) shared(epsilon,riDividedOnDiClause,workListCurr,workListNext,stats,graph) private(v,Nodes) reduction(+:activeVertices,error_total) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                __u32 degree;
                __u32 j;
                __u32 u;
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
                if(error >= epsilon)
                {
                    stats->pageRanks[v] = newPageRank;
                    Nodes = graph->vertices[v].outNodes;
                    degree = graph->vertices[v].out_degree;
                    for(j = 0 ; j < (degree) ; j++)
                    {
                        u = Nodes->edges_array_dest[j];

                        #pragma omp atomic write
                        workListNext[u] = 1;
                        // __u8 old_val = workListNext[u];
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

struct PageRankStats *pageRankDataDrivenPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 degree;
    __u32 j;
    __u32 u;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;
    struct EdgeList *Nodes;

    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);

    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
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

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjArrayList(double epsilon,  __u32 iterations, struct GraphAdjArrayList *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 degree;
    __u32 j;
    __u32 u;
    struct EdgeList *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;

    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
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
                    __u32 u = Nodes->edges_array_dest[j];
                    float prevResidual = 0.0f;

                    prevResidual = aResiduals[u];

                    #pragma omp atomic update
                    aResiduals[u] += delta;

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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


struct PageRankStats *pageRankGraphAdjLinkedList(double epsilon,  __u32 iterations, __u32 pushpull, struct GraphAdjLinkedList *graph)
{

    struct PageRankStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // pull
        stats = pageRankPullGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 1: // push
        stats = pageRankPushGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 2: // pull
        stats = pageRankPullFixedPointGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 3: // push
        stats = pageRankPushFixedPointGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 4: // pull
        stats = pageRankDataDrivenPullGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 5: // push
        stats = pageRankDataDrivenPushGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    case 6: // pullpush
        stats = pageRankDataDrivenPullPushGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    default:// push
        stats = pageRankPullGraphAdjLinkedList(epsilon, iterations, graph);
        break;
    }


    return stats;

}

struct PageRankStats *pageRankPullGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 activeVertices = 0;

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
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

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
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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
            __u32 degree = graph->vertices[v].out_degree;
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = 0 ; j < (degree) ; j++)
            {
                __u32 u = Nodes->dest;
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * pageRanksNext[v]);
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPullFixedPointGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 j;
    __u32 v;
    __u32 u;
    __u32 degree;
    __u32 activeVertices = 0;

    struct AdjLinkedListNode *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));




    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for default(none) private(v) shared(graph,pageRanksNext)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {
            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankPushFixedPointGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;

    // double error = 0;
    __u32 activeVertices = 0;

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



    __u64 *pageRanksNext = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));
    __u64 *riDividedOnDiClause = (__u64 *) my_malloc(graph->num_vertices * sizeof(__u64));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push FP (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
    printf(" -----------------------------------------------------\n");
    printf("| %-10s | %-8s | %-15s | %-9s | \n", "Iteration", "Active", "Error", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    #pragma omp parallel for default(none) private(v) shared(pageRanksNext,graph)
    for(v = 0; v < graph->num_vertices; v++)
    {

        pageRanksNext[v] = 0.0f;
    }

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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
            __u32 degree = graph->vertices[v].out_degree;
            // __u32 tid = omp_get_thread_num();
            __u32 j;

            for(j = 0 ; j < (degree) ; j++)
            {
                __u32  u = Nodes->dest;
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

        #pragma omp parallel for private(v) shared(epsilon, pageRanksNext,stats) reduction(+ : error_total, activeVertices)
        for(v = 0; v < graph->num_vertices; v++)
        {



            float prevPageRank =  stats->pageRanks[v];
            float nextPageRank =  stats->base_pr + (stats->damp * Fixed64ToDouble(pageRanksNext[v]));
            stats->pageRanks[v] = nextPageRank;
            pageRanksNext[v] = 0.0f;
            double error = fabs( nextPageRank - prevPageRank);
            error_total += (error / graph->num_vertices);

            if(error >= epsilon)
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

struct PageRankStats *pageRankDataDrivenPullGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 i;
    __u32 v;




    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;
    struct AdjLinkedListNode *Nodes;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);





    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
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

        #pragma omp parallel for default(none) shared(epsilon,riDividedOnDiClause,workListCurr,workListNext,stats,graph) private(v,Nodes) reduction(+:activeVertices,error_total) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            if(workListCurr[v])
            {

                __u32 degree;
                __u32 j;
                __u32 u;
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
                if(error >= epsilon)
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
                        // __u8 old_val = workListNext[u];
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

struct PageRankStats *pageRankDataDrivenPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 degree;
    __u32 j;
    __u32 u;



    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;
    struct AdjLinkedListNode *Nodes;



    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));



    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
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

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjLinkedList(double epsilon,  __u32 iterations, struct GraphAdjLinkedList *graph)
{

    double error_total = 0.0;
    __u32 v;
    __u32 degree;
    __u32 j;
    __u32 u;
    struct AdjLinkedListNode *Nodes;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    __u8 *workListCurr = NULL;
    __u8 *workListNext = NULL;
    int activeVertices = 0;


    workListCurr  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));
    workListNext  = (__u8 *) my_malloc(graph->num_vertices * sizeof(__u8));


    resetWorkList(workListNext, graph->num_vertices);
    resetWorkList(workListCurr, graph->num_vertices);




    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *aResiduals = (float *) my_malloc(graph->num_vertices * sizeof(float));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Push DD (tolerance/epsilon)");
    printf(" -----------------------------------------------------\n");
    printf("| %-51.13lf | \n", epsilon);
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

    for(stats->iterations = 0; stats->iterations < iterations; stats->iterations++)
    {
        Start(timer_inner);
        error_total = 0;
        activeVertices = 0;

        #pragma omp parallel for default(none) private(Nodes,degree,v,j,u) shared(stats,epsilon,graph,workListCurr,workListNext,aResiduals) reduction(+:error_total,activeVertices) schedule(dynamic,1024)
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

                    if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
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
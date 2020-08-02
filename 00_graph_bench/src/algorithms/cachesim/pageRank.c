#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>

#include "timer.h"
#include "myMalloc.h"
#include "boolean.h"
#include "arrayQueue.h"
#include "bitmap.h"
#include "worklist.h"

#include "graphConfig.h"

#include "fixedPoint.h"
#include "quantization.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjArrayList.h"
#include "graphAdjLinkedList.h"

#include "cache.h"
#include "pageRank_Kernels.h"
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

struct PageRankStats  *pageRankGraphGrid(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph)
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


struct PageRankStats *pageRankPullRowGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;

    uint32_t totalPartitions  = graph->grid->num_partitions;

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif


    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Row (tolerance/epsilon)");
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

#ifdef CACHE_HARNESS
        pageRankPullRowGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

#ifdef CPU_HARNESS
        pageRankPullRowGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif


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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}



struct PageRankStats *pageRankPullRowFixedPointGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Pull-Row FP (tolerance/epsilon)");
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
                riDividedOnDiClause[v] = 0;
        }


#ifdef CACHE_HARNESS
        pageRankPullRowFixedPointGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

#ifdef CPU_HARNESS
        pageRankPullRowFixedPointGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;

}



/******************************************************************/



struct PageRankStats *pageRankPushColumnGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));


    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push-Col (tolerance/epsilon)");
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


#ifdef CACHE_HARNESS
        pageRankPushColumnGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

#ifdef CPU_HARNESS
        pageRankPushColumnGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

    stats->error_total = error_total;
    return stats;
}



struct PageRankStats *pageRankPushColumnFixedPointGraphGrid(double epsilon,  uint32_t iterations, struct GraphGrid *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;


#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphGrid(graph);
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));



    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Page Rank Push-Col FP (tolerance/epsilon)");
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
                riDividedOnDiClause[v] = 0;
        }


#ifdef CACHE_HARNESS
        pageRankPushColumnFixedPointGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif

#ifdef CPU_HARNESS
        pageRankPushColumnFixedPointGraphGridKernel(riDividedOnDiClause, pageRanksNext,  graph->grid->partitions, totalPartitions);
#endif


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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

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


struct PageRankStats *pageRankGraphCSR(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph)
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
        stats = pageRankPullFixedPointGraphCSR(epsilon, iterations, graph);
        break;
    case 5: // push
        stats = pageRankPushFixedPointGraphCSR(epsilon, iterations, graph);
        break;

    case 6: // pull
        stats = pageRankDataDrivenPullGraphCSR(epsilon, iterations, graph);
        break;
    case 7: // push
        stats = pageRankDataDrivenPushGraphCSR(epsilon, iterations, graph);
        break;
    case 8: // pullpush
        stats = pageRankDataDrivenPullPushGraphCSR(epsilon, iterations, graph);
        break;

    // case 9: // push
    //     pageRankDataDrivenPullFixedPointGraphCSR(epsilon, iterations, graph);
    // break;
    // case 10: // pull
    //     pageRankDataDrivenPushFixedPointGraphCSR(epsilon, iterations, graph);
    // break;

    default:// pull
        stats = pageRankPullGraphCSR(epsilon, iterations, graph);
        break;
    }

    return stats;

}


// topoligy driven approach
struct PageRankStats *pageRankPullGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    // uint32_t j;
    uint32_t v;
    // uint32_t u;
    // uint32_t degree;
    // uint32_t edge_idx;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));
    float *pageRanksNext = (float *) my_malloc(graph->num_vertices * sizeof(float));
    float *riDividedOnDiClause = (float *) my_malloc(graph->num_vertices * sizeof(float));

#ifdef CACHE_HARNESS
    uint32_t numPropertyRegions = 2;
    struct PropertyMetaData *propertyMetaData = (struct PropertyMetaData *) my_malloc(graph->num_vertices * sizeof(struct PropertyMetaData));
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, numPropertyRegions);

    propertyMetaData[0].base_address = (uint64_t)&riDividedOnDiClause[0];
    propertyMetaData[0].size = graph->num_vertices;
    propertyMetaData[0].data_type_size = sizeof(float);

    propertyMetaData[1].base_address = (uint64_t)&pageRanksNext[0];
    propertyMetaData[1].size = graph->num_vertices;
    propertyMetaData[1].data_type_size = sizeof(float);

    initDoubleTaggedCacheRegion(cache, propertyMetaData);
#endif

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
            {
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            }
            else
            {
                riDividedOnDiClause[v] = 0.0f;
            }
        }


#ifdef CACHE_HARNESS
        pageRankPullGraphCSRKernelCache(cache, riDividedOnDiClause, pageRanksNext, vertices->out_degree, vertices->edges_idx, sorted_edges_array, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        pageRankPullGraphCSRKernel(riDividedOnDiClause, pageRanksNext, vertices->out_degree, vertices->edges_idx, sorted_edges_array, graph->num_vertices);
#endif


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



#ifdef CACHE_HARNESS
    printStatsDoubleTaggedCache(cache, graph->vertices->in_degree, graph->vertices->out_degree);
    freeDoubleTaggedCache(cache);
    free(propertyMetaData);
#endif


    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);
    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
{


    double error_total = 0.0;
    uint32_t v;

    // double error = 0;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif


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
            if(graph->vertices->out_degree[v])
                riDividedOnDiClause[v] = stats->pageRanks[v] / graph->vertices->out_degree[v];
            else
                riDividedOnDiClause[v] = 0;
        }


#ifdef CACHE_HARNESS
        pageRankPushGraphCSRKernelCache(cache, riDividedOnDiClause, pageRanksNext, graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        pageRankPushGraphCSRKernel(riDividedOnDiClause, pageRanksNext, graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest, graph->num_vertices);
#endif


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


    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;
}


// topoligy driven approach
struct PageRankStats *pageRankPullFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t v;
    uint32_t activeVertices = 0;

    // float init_pr = 1.0f / (float)graph->num_vertices;


    // uint64_t stats->base_pr_fp = FloatToFixed64(stats->base_pr);
    // uint64_t epsilon_fp = DoubleToFixed64(epsilon);
    // uint64_t num_vertices_fp = UInt32ToFixed64();
    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Vertex *vertices = NULL;
    uint32_t *sorted_edges_array = NULL;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif


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
                riDividedOnDiClause[v] = 0;
        }

        // Stop(timer_inner);
        // printf("|A %-9u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));

        //  Start(timer_inner);

#ifdef CACHE_HARNESS
        pageRankPullFixedPointGraphCSRKernelCache(cache, riDividedOnDiClause, pageRanksNext, vertices->out_degree, vertices->edges_idx, sorted_edges_array, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        pageRankPullFixedPointGraphCSRKernel(riDividedOnDiClause, pageRanksNext, vertices->out_degree, vertices->edges_idx, sorted_edges_array, graph->num_vertices);
#endif

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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;

}

struct PageRankStats *pageRankPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
{

    double error_total = 0.0;
    uint32_t v;

    uint32_t activeVertices = 0;

    struct PageRankStats *stats = newPageRankStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) malloc(sizeof(struct Timer));

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

    uint64_t *pageRanksNext = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));
    uint64_t *riDividedOnDiClause = (uint64_t *) my_malloc(graph->num_vertices * sizeof(uint64_t));


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
                riDividedOnDiClause[v] = 0;

        }
        // Stop(timer_inner);
        // printf("|A%-10u | %-8u | %-15.13lf | %-9f | \n",stats->iterations, activeVertices,error_total, Seconds(timer_inner));
        // Start(timer_inner);


#ifdef CACHE_HARNESS
        pageRankPushFixedPointGraphCSRKernelCache(cache, riDividedOnDiClause, pageRanksNext, graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        pageRankPushFixedPointGraphCSRKernel(riDividedOnDiClause, pageRanksNext, graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest, graph->num_vertices);
#endif
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


    free(timer);
    free(timer_inner);
    free(pageRanksNext);
    free(riDividedOnDiClause);

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;
}


struct PageRankStats *pageRankDataDrivenPullGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
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

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

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


#ifdef CACHE_HARNESS
        activeVertices += pageRankDataDrivenPullGraphCSRKernelCache(cache, riDividedOnDiClause, stats->pageRanks,
                          vertices->out_degree, vertices->edges_idx, sorted_edges_array,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        activeVertices += pageRankDataDrivenPullGraphCSRKernel(riDividedOnDiClause, stats->pageRanks,
                          vertices->out_degree, vertices->edges_idx, sorted_edges_array,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif


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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;
}

struct PageRankStats *pageRankDataDrivenPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
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

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

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

#ifdef CACHE_HARNESS
        activeVertices += pageRankDataDrivenPushGraphCSRKernelCache(cache, aResiduals, stats->pageRanks,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        activeVertices += pageRankDataDrivenPushGraphCSRKernel(aResiduals, stats->pageRanks,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif

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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;
}


struct PageRankStats *pageRankDataDrivenPullPushGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR *graph)
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

#ifdef CACHE_HARNESS
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, graph->num_vertices, POLICY, 0);
#endif

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


#ifdef CACHE_HARNESS
        activeVertices += pageRankDataDrivenPullPushGraphCSRKernelCache(cache, aResiduals, stats->pageRanks,
                          vertices->out_degree, vertices->edges_idx, sorted_edges_array,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif

#ifdef CPU_HARNESS
        activeVertices += pageRankDataDrivenPullPushGraphCSRKernel(aResiduals, stats->pageRanks,
                          vertices->out_degree, vertices->edges_idx, sorted_edges_array,
                          graph->vertices->out_degree, graph->vertices->edges_idx, graph->sorted_edges_array->edges_array_dest,
                          workListCurr, workListNext, &error_total, epsilon, graph->num_vertices);
#endif


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

#ifdef CACHE_HARNESS
    // printStats(cache->cold_cache);
    freeDoubleTaggedCache(cache);
#endif

    stats->error_total = error_total;
    return stats;

}


// float* pageRankDataDrivenPullFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph){


// }

// float* pageRankDataDrivenPullPushFixedPointGraphCSR(double epsilon,  uint32_t iterations, struct GraphCSR* graph){


// }



// ********************************************************************************************
// ***************          ArrayList DataStructure              **************
// ********************************************************************************************


struct PageRankStats *pageRankGraphAdjArrayList(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankPullGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankPullFixedPointGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankPushFixedPointGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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
                riDividedOnDiClause[v] = 0;

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

struct PageRankStats *pageRankDataDrivenPullGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankDataDrivenPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjArrayList(double epsilon,  uint32_t iterations, struct GraphAdjArrayList *graph)
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
                    uint32_t u = Nodes->edges_array_dest[j];
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


struct PageRankStats *pageRankGraphAdjLinkedList(double epsilon,  uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankPullGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankPullFixedPointGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankPushFixedPointGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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
                riDividedOnDiClause[v] = 0;

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

struct PageRankStats *pageRankDataDrivenPullGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankDataDrivenPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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

struct PageRankStats *pageRankDataDrivenPullPushGraphAdjLinkedList(double epsilon,  uint32_t iterations, struct GraphAdjLinkedList *graph)
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
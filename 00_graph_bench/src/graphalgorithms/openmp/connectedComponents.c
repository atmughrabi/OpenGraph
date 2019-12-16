// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : connectedComponents.c
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
#include <Judy.h>

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
#include "reorder.h"

#include "connectedComponents.h"


Pvoid_t JArray = (PWord_t) NULL; // Declare static hash table

// ********************************************************************************************
// ***************                  Stats DataStructure                          **************
// ********************************************************************************************



struct CCStats *newCCStatsGraphCSR(struct GraphCSR *graph)
{

    uint32_t v;

    struct CCStats *stats = (struct CCStats *) my_malloc(sizeof(struct CCStats));

    stats->iterations = 0;
    stats->neighbor_rounds = 2;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->components = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->counts = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->components[v] = v;
        stats->labels[v] = v;
        stats->counts[v] = 0;
    }

    return stats;

}
struct CCStats *newCCStatsGraphGrid(struct GraphGrid *graph)
{

    uint32_t v;

    struct CCStats *stats = (struct CCStats *) my_malloc(sizeof(struct CCStats));

    stats->neighbor_rounds = 2;
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->components = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->counts = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->components[v] = v;
        stats->labels[v] = v;
        stats->counts[v] = 0;
    }

    return stats;

}
struct CCStats *newCCStatsGraphAdjArrayList(struct GraphAdjArrayList *graph)
{
    uint32_t v;

    struct CCStats *stats = (struct CCStats *) my_malloc(sizeof(struct CCStats));

    stats->neighbor_rounds = 2;
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->components = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->counts = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->components[v] = v;
        stats->labels[v] = v;
        stats->counts[v] = 0;
    }
    return stats;

}
struct CCStats *newCCStatsGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{
    uint32_t v;

    struct CCStats *stats = (struct CCStats *) my_malloc(sizeof(struct CCStats));

    stats->neighbor_rounds = 2;
    stats->iterations = 0;
    stats->num_vertices = graph->num_vertices;
    stats->time_total = 0.0f;
    stats->components = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->counts = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    stats->labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));

    #pragma omp parallel for default(none) private(v) shared(stats)
    for(v = 0; v < stats->num_vertices; v++)
    {
        stats->components[v] = v;
        stats->labels[v] = v;
        stats->counts[v] = 0;
    }

    return stats;
}

void freeCCStats(struct CCStats *stats)
{

    if(stats)
    {
        if(stats->components)
            free(stats->components);
        if(stats->counts)
            free(stats->counts);
        if(stats->labels)
            free(stats->labels);
        free(stats);
    }
}

void printCCStats(struct CCStats *stats)
{

    Word_t *PValue;
    Word_t   Index;
    uint32_t k = 5;
    uint32_t numComp = 0;
    uint32_t i;

    for(i = 0; i < stats->num_vertices; i++)
    {
        addSample(stats->components[i]);
    }


    Index = 0;
    JLF(PValue, JArray, Index);
    while (PValue != NULL)
    {
        // printf("%lu %lu\n", Index, *PValue);
        stats->counts[Index] = *PValue;
        * PValue = 0;
        JLN(PValue, JArray, Index);

    }

    for(i = 0; i < stats->num_vertices; i++)
    {
        if(stats->counts[i])
            numComp++;
    }

    stats->labels = radixSortEdgesByDegree(stats->counts, stats->labels, stats->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Top Clusters", "Count");
    printf(" -----------------------------------------------------\n");

    for(i = (stats->num_vertices - 1); i > (stats->num_vertices - 1 - k); i--)
    {

        printf("| %-21u | %-27u | \n", stats->labels[i], stats->counts[i] );

    }
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27u | \n", "Num Components", numComp);
    printf(" -----------------------------------------------------\n");
}

void printComponents(struct CCStats *stats)
{

    uint32_t i;
    for(i = 0 ; i < stats->num_vertices; i++)
    {
        printf("v : %u comp : %u \n", i, stats->components[i]);
    }

}

// ********************************************************************************************
// ***************                       Helper Functions                        **************
// ********************************************************************************************

uint32_t atomicMin(uint32_t *oldValue, uint32_t newValue)
{

    uint32_t oldTemp;
    uint32_t flag = 0;

    do
    {
        oldTemp = *oldValue;
    }
    while(oldTemp > newValue && !(flag = __sync_bool_compare_and_swap(oldValue, oldTemp, newValue)));

    return flag;

}

void linkNodes(uint32_t u, uint32_t v, uint32_t *components)
{
    uint32_t p1 = components[u];
    uint32_t p2 = components[v];

    while(p1 != p2)
    {
        uint32_t high = p1 > p2 ? p1 : p2;
        uint32_t low = p1 + (p2 - high);
        uint32_t phigh = components[high];

        if ((phigh == low) ||
                (phigh == high && __sync_bool_compare_and_swap(&(components[high]), high, low)))
            break;
        p1 = components[components[high]];
        p2 = components[low];

    }

}


void compressNodes(uint32_t num_vertices, uint32_t *components)
{
    uint32_t n;
    #pragma omp parallel for schedule(dynamic, 2048)
    for (n = 0; n < num_vertices; n++)
    {
        while (components[n] != components[components[n]])
        {
            components[n] = components[components[n]];
        }
    }
}


void addSample(uint32_t id)
{
    Word_t *PValue;

    JLI(PValue, JArray, id);
    *PValue += 1;

}

uint32_t sampleFrequentNode(uint32_t num_vertices, uint32_t num_samples, uint32_t *components)
{

    Word_t *PValue;
    Word_t   Index;
    uint32_t i;
    initializeMersenneState (mt19937var, 27491095); 
    for (i = 0; i < num_samples; i++)
    {
        uint32_t n = generateRandInt(mt19937var) % num_vertices;
        addSample(components[n]);
    }

    uint32_t maxKey = 0;
    uint32_t maxCount = 0;

    Index = 0;
    JLF(PValue, JArray, Index);
    while (PValue != NULL)
    {
        // printf("%lu %lu\n", Index, *PValue);
        if(*PValue > maxCount)
        {
            maxCount = *PValue;
            maxKey = Index;

        }
        *PValue = 0;
        JLN(PValue, JArray, Index);

    }

    float fractiongraph = ((float)maxCount / num_samples);

    printf("| %-21s | %-27u | \n", "Skipping(%)", (int)fractiongraph * 100);



    return maxKey;
}

// ********************************************************************************************
// ***************                  CSR DataStructure                            **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphCSR(uint32_t iterations, uint32_t pushpull, struct GraphCSR *graph)
{

    struct CCStats *stats = NULL;

    switch (pushpull)
    {

    case 0: // Shiloach Vishkin
        stats = connectedComponentsShiloachVishkinGraphCSR( iterations, graph);
        break;
    case 1: // Afforest
        stats = connectedComponentsAfforestGraphCSR( iterations, graph);
        break;
    case 2: // WCC
        stats = connectedComponentsWeaklyGraphCSR( iterations, graph);
        break;
    default:// Afforest
        stats = connectedComponentsAfforestGraphCSR( iterations, graph);
        break;
    }

    return stats;

}

struct CCStats *connectedComponentsShiloachVishkinGraphCSR( uint32_t iterations, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t componentsCount = 0;
    uint32_t change = 0;

    struct CCStats *stats = newCCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Shiloach-Vishkin Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;

    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            degree = graph->vertices->out_degree[src];
            edge_idx = graph->vertices->edges_idx[src];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                dest = graph->sorted_edges_array->edges_array_dest[j];
                uint32_t comp_src = stats->components[src];
                uint32_t comp_dest = stats->components[dest];

                if(comp_src == comp_dest)
                    continue;

                uint32_t comp_high = comp_src > comp_dest ? comp_src : comp_dest;
                uint32_t comp_low = comp_src + (comp_dest - comp_high);

                if(comp_high == stats->components[comp_high])
                {
                    change = 1;
                    stats->components[comp_high] = comp_low;
                }
            }
        }


        compressNodes( stats->num_vertices, stats->components);

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);
    return stats;

}

struct CCStats *connectedComponentsAfforestGraphCSR( uint32_t iterations, struct GraphCSR *graph)
{

    uint32_t u;
    uint32_t componentsCount = 0;
    Word_t    Bytes;
    uint32_t num_samples = 1024;

    if(num_samples > graph->num_vertices)
        num_samples = graph->num_vertices / 2;

    struct CCStats *stats = newCCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));


    stats->neighbor_rounds = 2;


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Afforest Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Neighbor Round", "Time (S)");
    printf(" -----------------------------------------------------\n");

    uint32_t r = 0;

    Start(timer);
    for(r = 0; r < stats->neighbor_rounds; r++)
    {
        Start(timer_inner);
        #pragma omp parallel for schedule(dynamic, 2048)
        for(u = 0; u < graph->num_vertices; u++)
        {
            uint32_t j;
            uint32_t v;
            uint32_t degree_out =  graph->vertices->out_degree[u];
            uint32_t edge_idx_out =  graph->vertices->edges_idx[u];

            for(j = (edge_idx_out + r) ; j < (edge_idx_out + degree_out) ; j++)
            {
                v =  graph->sorted_edges_array->edges_array_dest[j];
                linkNodes(u, v, stats->components);
                break;
            }
        }
        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));

        Start(timer_inner);
        compressNodes(graph->num_vertices, stats->components);
        Stop(timer_inner);
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27f | \n", "", Seconds(timer_inner));
        printf(" -----------------------------------------------------\n");

    }// end neighbor_rounds loop


    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Sampling Components", "");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    uint32_t sampleComp = sampleFrequentNode(graph->num_vertices, num_samples,  stats->components);
    Stop(timer_inner);
    printf("| Most freq ID: %-7u | %-27f | \n", sampleComp, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Final Link Phase", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
#if DIRECTED
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;
        uint32_t degree_out;
        uint32_t degree_in;
        uint32_t edge_idx_out;
        uint32_t edge_idx_in;

        if(stats->components[u] == sampleComp)
            continue;

        degree_out =  graph->vertices->out_degree[u];
        edge_idx_out =  graph->vertices->edges_idx[u];

        for(j = (edge_idx_out + stats->neighbor_rounds) ; j < (edge_idx_out + degree_out) ; j++)
        {
            v =  graph->sorted_edges_array->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }

        degree_in =  graph->inverse_vertices->out_degree[u];
        edge_idx_in =  graph->inverse_vertices->edges_idx[u];

        for(j = (edge_idx_in) ; j < (edge_idx_in + degree_in) ; j++)
        {
            v =  graph->inverse_sorted_edges_array->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }

    }
#else
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;
        uint32_t degree_out;
        uint32_t edge_idx_out;

        if(stats->components[u] == sampleComp)
            continue;

        degree_out =  graph->vertices->out_degree[u];
        edge_idx_out =  graph->vertices->edges_idx[u];

        for(j = (edge_idx_out + stats->neighbor_rounds) ; j < (edge_idx_out + degree_out) ; j++)
        {
            v =  graph->sorted_edges_array->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }
    }
#endif
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", componentsCount, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    compressNodes(graph->num_vertices, stats->components);
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));
    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->neighbor_rounds, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);

    JSLFA(Bytes, JArray);
    return stats;

}

struct CCStats *connectedComponentsWeaklyGraphCSR( uint32_t iterations, struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t componentsCount = 0;
    uint32_t change = 0;

    struct CCStats *stats = newCCStatsGraphCSR(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Weakly Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;


    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v,degree,edge_idx) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            degree = graph->vertices->out_degree[src];
            edge_idx = graph->vertices->edges_idx[src];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                dest = graph->sorted_edges_array->edges_array_dest[j];

                if(atomicMin(&(stats->components[dest]), stats->components[src]))
                {
                    setBitAtomic(bitmapNext, dest);
                }

                if(atomicMin(&(stats->components[src]), stats->components[dest]))
                {
                    setBitAtomic(bitmapNext, src);
                }
            }
        }


        // compressNodes( stats->num_vertices, stats->components);

        #pragma omp parallel for reduction (+:change)
        for(v = 0 ; v < ((bitmapNext->size + kBitsPerWord - 1) / kBitsPerWord); v++)
        {
            change += bitmapNext->bitarray[v];
            bitmapNext->bitarray[v] = 0;
        }



        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    printCCStats(stats);
    // connectedComponentsVerifyGraphCSR(stats, graph);
    return stats;

}

uint32_t connectedComponentsVerifyGraphCSR(struct CCStats *stats, struct GraphCSR *graph)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Connected Components Verification");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    uint32_t pass = 1;
    struct ArrayQueue *frontier = newArrayQueue(graph->num_vertices);
    uint32_t *inverselabels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    uint32_t iter;
    uint32_t j;
    uint32_t i;

    for(iter = 0; iter < stats->num_vertices; iter++)
    {
        inverselabels[stats->components[iter]] = iter;
    }

    uint32_t n;
    uint32_t comp;

    for(iter = 0 ; iter < graph->num_vertices; iter++)
    {

        comp = stats->components[iter];
        n = inverselabels[comp];

        softResetArrayQueue(frontier);
        enArrayQueueWithBitmap(frontier, n);

        for(i = frontier->head ; i < frontier->tail; i++)
        {
            uint32_t u = frontier->queue[i];
            uint32_t edge_idx = graph->vertices->edges_idx[u];
            uint32_t out_degree = graph->vertices->out_degree[u];

            for(j = edge_idx ; j < (edge_idx + out_degree) ; j++)
            {
                uint32_t v = graph->sorted_edges_array->edges_array_dest[j];

                if(stats->components[v] != comp)
                {
                    pass = 0;
                }
                if(!isEnArrayQueued(frontier, v))
                    enArrayQueueWithBitmap(frontier, v);
            }

#if DIRECTED

            uint32_t in_edge_idx = graph->inverse_vertices->edges_idx[u];
            uint32_t in_degree = graph->inverse_vertices->out_degree[u];

            for(j = in_edge_idx ; j < (in_edge_idx + in_degree) ; j++)
            {
                uint32_t v = graph->inverse_sorted_edges_array->edges_array_dest[j];

                if(stats->components[v] != comp)
                {
                    pass = 0;
                }
                if(!isEnArrayQueued(frontier, v))
                    enArrayQueueWithBitmap(frontier, v);
            }


#endif


        }
    }


    for(iter = 0 ; iter < (frontier->q_bitmap->size); iter++)
    {
        if(!getBit(frontier->q_bitmap, iter))
            pass++;
    }

    if(!pass)
    {

        printf("PASS\n");
        pass = 1;
    }
    else
    {

        printf("FAIL %u\n", pass);
        pass = 0;


    }




    free(inverselabels);
    freeArrayQueue(frontier);

    return pass;
}


// ********************************************************************************************
// ***************                  GRID DataStructure                           **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphGrid(uint32_t iterations, uint32_t pushpull, struct GraphGrid *graph)
{

    struct CCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // Shiloach Vishkin
        stats = connectedComponentsShiloachVishkinGraphGrid( iterations, graph);
        break;
    case 1: // Afforest
        stats = connectedComponentsAfforestGraphGrid( iterations, graph);
        break;
    case 2: // Weakly Connected
        stats = connectedComponentsWeaklyGraphGrid( iterations, graph);
        break;
    default:// Afforest
        stats = connectedComponentsWeaklyGraphGrid( iterations, graph);
        break;
    }

    return stats;

}

struct CCStats *connectedComponentsShiloachVishkinGraphGrid( uint32_t iterations, struct GraphGrid *graph)
{

    uint32_t i;
    uint32_t componentsCount = 0;
    uint32_t change = 0;
    uint32_t totalPartitions  = graph->grid->num_partitions;
    struct CCStats *stats = newCCStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Shiloach-Vishkin Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;

    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(i) schedule (dynamic,numThreads)
        for (i = 0; i < totalPartitions; ++i)
        {
            uint32_t j;
            // #pragma omp parallel for private(j) schedule (dynamic,numThreads)
            for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
            {
                uint32_t k;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    uint32_t src = partition->edgeList->edges_array_src[k];
                    uint32_t dest = partition->edgeList->edges_array_dest[k];


                    uint32_t comp_src = stats->components[src];
                    uint32_t comp_dest = stats->components[dest];

                    if(comp_src != comp_dest)
                    {
                        uint32_t comp_high = comp_src > comp_dest ? comp_src : comp_dest;
                        uint32_t comp_low = comp_src + (comp_dest - comp_high);

                        if(comp_high == stats->components[comp_high])
                        {
                            change = 1;
                            stats->components[comp_high] = comp_low;
                        }
                    }
                }
            }
        }

        compressNodes( stats->num_vertices, stats->components);

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);
    return stats;

}


struct CCStats *connectedComponentsAfforestGraphGrid( uint32_t iterations, struct GraphGrid *graph)
{

    uint32_t i;
    uint32_t v;
    uint32_t componentsCount = 0;
    Word_t    Bytes;
    uint32_t num_samples = 1024;
    uint32_t totalPartitions  = graph->grid->num_partitions;

    if(num_samples > graph->num_vertices)
        num_samples = graph->num_vertices / 2;

    struct CCStats *stats = newCCStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));
    uint32_t *neighbor = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    struct Bitmap *linked = newBitmap(graph->num_vertices);

    stats->neighbor_rounds = 2;
    #pragma omp parallel for default(none) private(v) shared(graph,neighbor)
    for(v = 0; v < graph->num_vertices; v++)
    {
        neighbor[v] = 0;
    }


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Afforest Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Neighbor Round", "Time (S)");
    printf(" -----------------------------------------------------\n");

    uint32_t r = 0;

    Start(timer);
    for(r = 0; r < stats->neighbor_rounds; r++)
    {
        Start(timer_inner);
        #pragma omp parallel for private(i) schedule (dynamic,numThreads)
        for (i = 0; i < totalPartitions; ++i)
        {
            uint32_t j;
            // #pragma omp parallel for private(j) schedule (dynamic,numThreads)
            for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
            {
                uint32_t k;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {
                    uint32_t src = partition->edgeList->edges_array_src[k];
                    uint32_t dest = partition->edgeList->edges_array_dest[k];

                    if(neighbor[src] >= r && !getBit(linked, src))
                    {
                        linkNodes(src, dest, stats->components);
                        setBit(linked, src);
                    }
                    else
                    {
                        neighbor[src]++;
                    }
                }
            }
        }
        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));

        #pragma omp parallel for default(none) private(v) shared(stats,neighbor)
        for(v = 0; v < stats->num_vertices; v++)
        {
            neighbor[v] = 0;
        }
        clearBitmap(linked);

        Start(timer_inner);
        compressNodes(graph->num_vertices, stats->components);
        Stop(timer_inner);
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27f | \n", "", Seconds(timer_inner));
        printf(" -----------------------------------------------------\n");

    }// end neighbor_rounds loop


    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Sampling Components", "");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    uint32_t sampleComp = sampleFrequentNode(graph->num_vertices, num_samples,  stats->components);
    Stop(timer_inner);
    printf("| Most freq ID: %-7u | %-27f | \n", sampleComp, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Final Link Phase", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
#if DIRECTED
    #pragma omp parallel for private(i) schedule (dynamic,numThreads)
    for (i = 0; i < totalPartitions; ++i)
    {
        uint32_t j;
        // #pragma omp parallel for private(j) schedule (dynamic,numThreads)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
        {
            uint32_t k;
            struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
            for (k = 0; k < partition->num_edges; ++k)
            {
                uint32_t src = partition->edgeList->edges_array_src[k];
                uint32_t dest = partition->edgeList->edges_array_dest[k];

                if(stats->components[src] != sampleComp)
                {

                    if(neighbor[src] >= stats->neighbor_rounds)
                    {
                        linkNodes(src, dest, stats->components);
                    }
                    else
                    {
                        neighbor[src]++;
                    }

                }

                if(stats->components[dest] != sampleComp)
                {
                    linkNodes(dest, src, stats->components);
                }

            }
        }
    }
#else
    #pragma omp parallel for private(i) schedule (dynamic,numThreads)
    for (i = 0; i < totalPartitions; ++i)
    {
        uint32_t j;
        // #pragma omp parallel for private(j) schedule (dynamic,numThreads)
        for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
        {
            uint32_t k;
            struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
            for (k = 0; k < partition->num_edges; ++k)
            {
                uint32_t src = partition->edgeList->edges_array_src[k];
                uint32_t dest = partition->edgeList->edges_array_dest[k];

                if(stats->components[src] != sampleComp)
                {

                    if(neighbor[src] >= stats->neighbor_rounds)
                    {
                        linkNodes(src, dest, stats->components);
                    }
                    else
                    {
                        neighbor[src]++;
                    }
                }
            }
        }
    }
#endif
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", componentsCount, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    compressNodes(graph->num_vertices, stats->components);
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));
    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->neighbor_rounds, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    free(neighbor);
    printCCStats(stats);
    freeBitmap(linked);
    JSLFA(Bytes, JArray);
    return stats;

}


struct CCStats *connectedComponentsWeaklyGraphGrid(uint32_t iterations, struct GraphGrid *graph)
{

    uint32_t v;
    uint32_t componentsCount = 0;
    uint32_t change = 0;
    uint32_t totalPartitions  = graph->grid->num_partitions;

    struct CCStats *stats = newCCStatsGraphGrid(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Weakly Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    stats->iterations = 0;
    change = 1;


    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        uint32_t i;
        #pragma omp parallel for private(i) schedule (dynamic,numThreads)
        for (i = 0; i < totalPartitions; ++i)
        {
            uint32_t j;
            // #pragma omp parallel for private(j) schedule (dynamic,numThreads)
            for (j = 0; j < totalPartitions; ++j)  // iterate over partitions colwise
            {
                uint32_t k;
                struct Partition *partition = &graph->grid->partitions[(i * totalPartitions) + j];
                for (k = 0; k < partition->num_edges; ++k)
                {

                    uint32_t src = partition->edgeList->edges_array_src[k];
                    uint32_t dest = partition->edgeList->edges_array_dest[k];

                    if(atomicMin(&(stats->components[dest]), stats->components[src]))
                    {
                        setBitAtomic(bitmapNext, dest);
                    }

                    if(atomicMin(&(stats->components[src]), stats->components[dest]))
                    {
                        setBitAtomic(bitmapNext, src);
                    }

                }
            }
        }


        // compressNodes( stats->num_vertices, stats->components);
        #pragma omp parallel for reduction (+:change)
        for(v = 0 ; v < ((bitmapNext->size + kBitsPerWord - 1) / kBitsPerWord); v++)
        {
            change += bitmapNext->bitarray[v];
            bitmapNext->bitarray[v] = 0;
        }

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    printCCStats(stats);
    // connectedComponentsVerifyGraphCSR(stats, graph);
    return stats;



}


// ********************************************************************************************
// ***************                  ArrayList DataStructure                      **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjArrayList(uint32_t iterations, uint32_t pushpull, struct GraphAdjArrayList *graph)
{

    struct CCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // Shiloach Vishkin
        stats = connectedComponentsShiloachVishkinGraphAdjArrayList( iterations, graph);
        break;
    case 1: // Afforest
        stats = connectedComponentsAfforestGraphAdjArrayList( iterations, graph);
        break;
    case 2: // Weakly Connected
        stats = connectedComponentsWeaklyGraphAdjArrayList( iterations, graph);
        break;
    default:// Afforest
        stats = connectedComponentsAfforestGraphAdjArrayList( iterations, graph);
        break;
    }

    return stats;

}

struct CCStats *connectedComponentsShiloachVishkinGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph)
{
    uint32_t v;
    uint32_t degree;
    uint32_t componentsCount = 0;
    uint32_t change = 0;
    struct EdgeList *Nodes;

    struct CCStats *stats = newCCStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Shiloach-Vishkin Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;

    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            Nodes = graph->vertices[v].outNodes;
            degree = graph->vertices[v].out_degree;

            for(j = 0 ; j < (degree) ; j++)
            {
                dest = Nodes->edges_array_dest[j];
                uint32_t comp_src = stats->components[src];
                uint32_t comp_dest = stats->components[dest];

                if(comp_src == comp_dest)
                    continue;

                uint32_t comp_high = comp_src > comp_dest ? comp_src : comp_dest;
                uint32_t comp_low = comp_src + (comp_dest - comp_high);

                if(comp_high == stats->components[comp_high])
                {
                    change = 1;
                    stats->components[comp_high] = comp_low;
                }
            }
        }


        compressNodes( stats->num_vertices, stats->components);

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);
    return stats;


}
struct CCStats *connectedComponentsAfforestGraphAdjArrayList(uint32_t iterations, struct GraphAdjArrayList *graph)
{

    uint32_t u;
    uint32_t componentsCount = 0;
    Word_t    Bytes;
    uint32_t num_samples = 1024;

    if(num_samples > graph->num_vertices)
        num_samples = graph->num_vertices / 2;

    struct CCStats *stats = newCCStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));


    stats->neighbor_rounds = 2;


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Afforest Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Neighbor Round", "Time (S)");
    printf(" -----------------------------------------------------\n");

    uint32_t r = 0;

    Start(timer);
    for(r = 0; r < stats->neighbor_rounds; r++)
    {
        Start(timer_inner);
        #pragma omp parallel for schedule(dynamic, 2048)
        for(u = 0; u < graph->num_vertices; u++)
        {
            uint32_t j;
            uint32_t v;

            struct EdgeList *Nodes = graph->vertices[u].outNodes;
            uint32_t degree_out = graph->vertices[u].out_degree;

            for(j = (0 + r) ; j < (degree_out) ; j++)
            {
                v =  Nodes->edges_array_dest[j];
                linkNodes(u, v, stats->components);
                break;
            }
        }
        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));

        Start(timer_inner);
        compressNodes(graph->num_vertices, stats->components);
        Stop(timer_inner);
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27f | \n", "", Seconds(timer_inner));
        printf(" -----------------------------------------------------\n");

    }// end neighbor_rounds loop


    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Sampling Components", "");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    uint32_t sampleComp = sampleFrequentNode(graph->num_vertices, num_samples,  stats->components);
    Stop(timer_inner);
    printf("| Most freq ID: %-7u | %-27f | \n", sampleComp, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Final Link Phase", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
#if DIRECTED
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;

        if(stats->components[u] == sampleComp)
            continue;

        struct EdgeList *Nodes_out = graph->vertices[u].outNodes;
        uint32_t degree_out = graph->vertices[u].out_degree;

        for(j = ( 0 + stats->neighbor_rounds) ; j < (degree_out) ; j++)
        {
            v =  Nodes_out->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }

        struct EdgeList *Nodes_in = graph->vertices[u].inNodes;
        uint32_t degree_in = graph->vertices[u].in_degree;

        for(j = (0) ; j < (degree_in) ; j++)
        {
            v = Nodes_in->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }

    }
#else
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;

        if(stats->components[u] == sampleComp)
            continue;

        struct EdgeList *Nodes_out = graph->vertices[u].outNodes;
        uint32_t degree_out = graph->vertices[u].out_degree;

        for(j = ( 0 + stats->neighbor_rounds) ; j < (degree_out) ; j++)
        {
            v =  Nodes_out->edges_array_dest[j];
            linkNodes(u, v, stats->components);
        }
    }
#endif
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", componentsCount, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    compressNodes(graph->num_vertices, stats->components);
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));
    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->neighbor_rounds, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);

    JSLFA(Bytes, JArray);
    return stats;


}
struct CCStats *connectedComponentsWeaklyGraphAdjArrayList( uint32_t iterations, struct GraphAdjArrayList *graph)
{

    uint32_t v;
    uint32_t componentsCount = 0;
    uint32_t change = 0;

    struct CCStats *stats = newCCStatsGraphAdjArrayList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Weakly Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;


    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            struct EdgeList *Nodes_out = graph->vertices[v].outNodes;
            uint32_t degree_out = graph->vertices[v].out_degree;

            for(j = 0 ; j < (degree_out) ; j++)
            {
                dest = Nodes_out->edges_array_dest[j];

                if(atomicMin(&(stats->components[dest]), stats->components[src]))
                {
                    setBitAtomic(bitmapNext, dest);
                }

                if(atomicMin(&(stats->components[src]), stats->components[dest]))
                {
                    setBitAtomic(bitmapNext, src);
                }
            }
        }


        // compressNodes( stats->num_vertices, stats->components);

        #pragma omp parallel for reduction (+:change)
        for(v = 0 ; v < ((bitmapNext->size + kBitsPerWord - 1) / kBitsPerWord); v++)
        {
            change += bitmapNext->bitarray[v];
            bitmapNext->bitarray[v] = 0;
        }



        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    printCCStats(stats);
    // connectedComponentsVerifyGraphCSR(stats, graph);
    return stats;

}

// ********************************************************************************************
// ***************                  LinkedList DataStructure                     **************
// ********************************************************************************************

struct CCStats *connectedComponentsGraphAdjLinkedList(uint32_t iterations, uint32_t pushpull, struct GraphAdjLinkedList *graph)
{

    struct CCStats *stats = NULL;

    switch (pushpull)
    {
    case 0: // Shiloach Vishkin
        stats = connectedComponentsShiloachVishkinGraphAdjLinkedList( iterations, graph);
        break;
    case 1: // Afforest
        stats = connectedComponentsAfforestGraphAdjLinkedList( iterations, graph);
        break;
    case 2: // Weakly Connected
        stats = connectedComponentsWeaklyGraphAdjLinkedList( iterations, graph);
        break;
    default:// Afforest
        stats = connectedComponentsAfforestGraphAdjLinkedList( iterations, graph);
        break;
    }

    return stats;

}

struct CCStats *connectedComponentsShiloachVishkinGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph)
{
    uint32_t v;
    uint32_t degree;
    uint32_t componentsCount = 0;
    uint32_t change = 0;
    struct AdjLinkedListNode *Nodes;

    struct CCStats *stats = newCCStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Shiloach-Vishkin Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;

    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v,degree,Nodes) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            Nodes = graph->vertices[v].outNodes;
            degree = graph->vertices[v].out_degree;

            for(j = 0 ; j < (degree) ; j++)
            {
                
                dest = Nodes->dest;
                Nodes = Nodes->next;

                uint32_t comp_src = stats->components[src];
                uint32_t comp_dest = stats->components[dest];

                if(comp_src == comp_dest)
                    continue;

                uint32_t comp_high = comp_src > comp_dest ? comp_src : comp_dest;
                uint32_t comp_low = comp_src + (comp_dest - comp_high);

                if(comp_high == stats->components[comp_high])
                {
                    change = 1;
                    stats->components[comp_high] = comp_low;
                }
            }
        }


        compressNodes( stats->num_vertices, stats->components);

        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);
    return stats;
}
struct CCStats *connectedComponentsAfforestGraphAdjLinkedList(uint32_t iterations, struct GraphAdjLinkedList *graph)
{

    uint32_t u;
    uint32_t componentsCount = 0;
    Word_t    Bytes;
    uint32_t num_samples = 1024;

    if(num_samples > graph->num_vertices)
        num_samples = graph->num_vertices / 2;

    struct CCStats *stats = newCCStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));


    stats->neighbor_rounds = 2;


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Afforest Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Neighbor Round", "Time (S)");
    printf(" -----------------------------------------------------\n");

    uint32_t r = 0;

    Start(timer);
    for(r = 0; r < stats->neighbor_rounds; r++)
    {
        Start(timer_inner);
        #pragma omp parallel for schedule(dynamic, 2048)
        for(u = 0; u < graph->num_vertices; u++)
        {
            uint32_t j;
            uint32_t v;

            struct AdjLinkedListNode *Nodes = graph->vertices[u].outNodes;
            uint32_t degree_out = graph->vertices[u].out_degree;

            for(j = (0 + r) ; j < (degree_out) ; j++)
            {
                v = Nodes->dest;
                Nodes = Nodes->next;

                linkNodes(u, v, stats->components);
                break;
            }
        }
        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));

        Start(timer_inner);
        compressNodes(graph->num_vertices, stats->components);
        Stop(timer_inner);
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
        printf(" -----------------------------------------------------\n");
        printf("| %-21s | %-27f | \n", "", Seconds(timer_inner));
        printf(" -----------------------------------------------------\n");

    }// end neighbor_rounds loop


    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Sampling Components", "");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    uint32_t sampleComp = sampleFrequentNode(graph->num_vertices, num_samples,  stats->components);
    Stop(timer_inner);
    printf("| Most freq ID: %-7u | %-27f | \n", sampleComp, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Final Link Phase", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
#if DIRECTED
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;

        if(stats->components[u] == sampleComp)
            continue;

        struct AdjLinkedListNode *Nodes_out = graph->vertices[u].outNodes;
        uint32_t degree_out = graph->vertices[u].out_degree;

        for(j = ( 0 + stats->neighbor_rounds) ; j < (degree_out) ; j++)
        {
            v = Nodes_out->dest;
            Nodes_out = Nodes_out->next;

            linkNodes(u, v, stats->components);
        }

        struct AdjLinkedListNode *Nodes_in = graph->vertices[u].inNodes;
        uint32_t degree_in = graph->vertices[u].in_degree;

        for(j = (0) ; j < (degree_in) ; j++)
        {
            v = Nodes_in->dest;
            Nodes_in = Nodes_in->next;

            linkNodes(u, v, stats->components);
        }

    }
#else
    #pragma omp parallel for schedule(dynamic, 2048)
    for(u = 0; u < graph->num_vertices; u++)
    {
        uint32_t j;
        uint32_t v;

        if(stats->components[u] == sampleComp)
            continue;

        struct AdjLinkedListNode *Nodes_out = graph->vertices[u].outNodes;
        uint32_t degree_out = graph->vertices[u].out_degree;

        for(j = ( 0 + stats->neighbor_rounds) ; j < (degree_out) ; j++)
        {
            v = Nodes_out->dest;
            Nodes_out = Nodes_out->next;


            linkNodes(u, v, stats->components);
        }
    }
#endif
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", componentsCount, Seconds(timer_inner));

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Compress", "Time (S)");
    printf(" -----------------------------------------------------\n");
    Start(timer_inner);
    compressNodes(graph->num_vertices, stats->components);
    Stop(timer_inner);
    printf("| %-21u | %-27f | \n", r, Seconds(timer_inner));
    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->neighbor_rounds, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);

    printCCStats(stats);

    JSLFA(Bytes, JArray);
    return stats;


}
struct CCStats *connectedComponentsWeaklyGraphAdjLinkedList( uint32_t iterations, struct GraphAdjLinkedList *graph)
{

    uint32_t v;
    uint32_t componentsCount = 0;
    uint32_t change = 0;

    struct CCStats *stats = newCCStatsGraphAdjLinkedList(graph);
    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));
    struct Timer *timer_inner = (struct Timer *) my_malloc(sizeof(struct Timer));

    struct Bitmap *bitmapNext = newBitmap(graph->num_vertices);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Weakly Connected Components");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27s | \n", "Iteration", "Time (S)");
    printf(" -----------------------------------------------------\n");


    Start(timer);
    stats->iterations = 0;
    change = 1;


    while(change)
    {
        Start(timer_inner);
        change = 0;
        stats->iterations++;

        #pragma omp parallel for private(v) schedule(dynamic, 1024)
        for(v = 0; v < graph->num_vertices; v++)
        {
            uint32_t j;
            uint32_t src = v;
            uint32_t dest;

            struct AdjLinkedListNode *Nodes_out = graph->vertices[v].outNodes;
            uint32_t degree_out = graph->vertices[v].out_degree;

            for(j = 0 ; j < (degree_out) ; j++)
            {
                dest = Nodes_out->dest;
                Nodes_out = Nodes_out->next;

                if(atomicMin(&(stats->components[dest]), stats->components[src]))
                {
                    setBitAtomic(bitmapNext, dest);
                }

                if(atomicMin(&(stats->components[src]), stats->components[dest]))
                {
                    setBitAtomic(bitmapNext, src);
                }
            }
        }


        // compressNodes( stats->num_vertices, stats->components);

        #pragma omp parallel for reduction (+:change)
        for(v = 0 ; v < ((bitmapNext->size + kBitsPerWord - 1) / kBitsPerWord); v++)
        {
            change += bitmapNext->bitarray[v];
            bitmapNext->bitarray[v] = 0;
        }



        Stop(timer_inner);
        printf("| %-21u | %-27f | \n", stats->iterations, Seconds(timer_inner));
    }

    Stop(timer);
    stats->time_total = Seconds(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Iterations", "Components", "Time (S)");
    printf(" -----------------------------------------------------\n");
    printf("| %-15u | %-15u | %-15f | \n", stats->iterations, componentsCount, stats->time_total);
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(timer_inner);
    freeBitmap(bitmapNext);
    printCCStats(stats);
    // connectedComponentsVerifyGraphCSR(stats, graph);
    return stats;

}
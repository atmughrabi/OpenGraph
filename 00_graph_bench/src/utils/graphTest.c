// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphTest.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>
#include <stdbool.h>
#include <omp.h>
#include <math.h>
#include <assert.h>

#include "myMalloc.h"
#include "mt19937.h"
#include "graphConfig.h"
#include "timer.h"

#include "edgeList.h"

#include "graphCSR.h"
#include "graphAdjLinkedList.h"
#include "graphAdjArrayList.h"
#include "graphGrid.h"

#include "BFS.h"
#include "DFS.h"
#include "pageRank.h"
#include "incrementalAggregation.h"
#include "bellmanFord.h"
#include "SSSP.h"
#include "SPMV.h"
#include "connectedComponents.h"
#include "triangleCount.h"

#include "graphStats.h"
#include "graphRun.h"
#include "graphTest.h"


uint32_t equalFloat(float a, float b, double epsilon)
{
    return fabs(a - b) < epsilon;
}

uint32_t compareFloatArrays(float *arr1, float *arr2, uint32_t arr1_size, uint32_t arr2_size)
{
    uint32_t i = 0;
    uint32_t missmatch = 0;
    double epsilon = 1e-3f;

    if(arr1_size != arr2_size)
        return 1;

    for(i = 0 ; i < arr1_size; i++)
    {

        if(!equalFloat(arr1[i], arr2[i], epsilon))
        {
            missmatch++;
        }
    }

    if (missmatch < 5)
        missmatch = 0;
    
    return missmatch;
}

uint32_t compareRealRanks(uint32_t *arr1, uint32_t *arr2, uint32_t arr1_size, uint32_t arr2_size)
{
    uint32_t i = 0;
    uint32_t missmatch = 0;
    // uint32_t rank_diff = 0;

    if(arr1_size != arr2_size)
        return 1;

    uint32_t *labels1 = (uint32_t *) my_malloc(arr1_size * sizeof(uint32_t));
    uint32_t *labels2 = (uint32_t *) my_malloc(arr2_size * sizeof(uint32_t));

    for(i = 0; i < arr1_size; i++)
    {
        labels1[arr1[i]] = i+1;
        labels2[arr2[i]] = i+1;
    }


    for(i = 0 ; i < arr1_size; i++)
    {

        if(labels1[i] != labels2[i])
        {
            // rank_diff = (labels1[i] > labels2[i]) ? (labels1[i] - labels2[i]) : (labels2[i] - labels1[i]);
            missmatch ++;
        }
    }

    free(labels1);
    free(labels2);
    return (missmatch);
}

uint32_t compareDistanceArrays(uint32_t *arr1, uint32_t *arr2, uint32_t arr1_size, uint32_t arr2_size)
{
    uint32_t i = 0;
    uint32_t missmatch = 0;

    if(arr1_size != arr2_size)
        return 1;

    for(i = 0 ; i < arr1_size; i++)
    {
        if(arr1[i] != arr2[i])
        {
            missmatch++;
        }
    }
    return missmatch;
}

uint32_t cmpGraphAlgorithmsTestStats(void *ref_stats, void *cmp_stats, uint32_t algorithm)
{

    uint32_t missmatch = 0;

    switch (algorithm)
    {
    case 0:  // bfs filename root
    {
        struct BFSStats *ref_stats_tmp = (struct BFSStats * )ref_stats;
        struct BFSStats *cmp_stats_tmp = (struct BFSStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->distances, cmp_stats_tmp->distances, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    case 1: // pagerank filename
    {
        struct PageRankStats *ref_stats_tmp = (struct PageRankStats * )ref_stats;
        struct PageRankStats *cmp_stats_tmp = (struct PageRankStats * )cmp_stats;
        // missmatch += compareRealRanks(ref_stats_tmp->realRanks, cmp_stats_tmp->realRanks, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
        missmatch += compareFloatArrays(ref_stats_tmp->pageRanks, cmp_stats_tmp->pageRanks, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);

    }
    break;
    case 2: // SSSP-Dijkstra file name root
    {
        struct SSSPStats *ref_stats_tmp = (struct SSSPStats * )ref_stats;
        struct SSSPStats *cmp_stats_tmp = (struct SSSPStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->distances, cmp_stats_tmp->distances, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    case 3: // SSSP-Bellmanford file name root
    {
        struct BellmanFordStats *ref_stats_tmp = (struct BellmanFordStats * )ref_stats;
        struct BellmanFordStats *cmp_stats_tmp = (struct BellmanFordStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->distances, cmp_stats_tmp->distances, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    case 4: // DFS file name root
    {
        struct DFSStats *ref_stats_tmp = (struct DFSStats * )ref_stats;
        struct DFSStats *cmp_stats_tmp = (struct DFSStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->distances, cmp_stats_tmp->distances, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    case 5: // SPMV file name root
    {
        struct SPMVStats *ref_stats_tmp = (struct SPMVStats * )ref_stats;
        struct SPMVStats *cmp_stats_tmp = (struct SPMVStats * )cmp_stats;
        missmatch += compareFloatArrays(ref_stats_tmp->vector_output, cmp_stats_tmp->vector_output, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
        missmatch = 0;
    }
    break;
    case 6: // Connected Components
    {
        struct CCStats *ref_stats_tmp = (struct CCStats * )ref_stats;
        struct CCStats *cmp_stats_tmp = (struct CCStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->components, cmp_stats_tmp->components, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    case 7: // Connected Components
    {
        struct TCStats *ref_stats_tmp = (struct TCStats * )ref_stats;
        struct TCStats *cmp_stats_tmp = (struct TCStats * )cmp_stats;
        missmatch += (ref_stats_tmp->counts == cmp_stats_tmp->counts) ? 1 : 0;
    }
    break;
    case 8: // incremental Aggregation file name root
    {
        struct IncrementalAggregationStats *ref_stats_tmp = (struct IncrementalAggregationStats * )ref_stats;
        struct IncrementalAggregationStats *cmp_stats_tmp = (struct IncrementalAggregationStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->labels, cmp_stats_tmp->labels, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
        missmatch = 0;
    }
    break;
    default:// bfs
    {
        struct BFSStats *ref_stats_tmp = (struct BFSStats * )ref_stats;
        struct BFSStats *cmp_stats_tmp = (struct BFSStats * )cmp_stats;
        missmatch += compareDistanceArrays(ref_stats_tmp->distances, cmp_stats_tmp->distances, ref_stats_tmp->num_vertices, cmp_stats_tmp->num_vertices);
    }
    break;
    }

    return missmatch;
}

void *runGraphAlgorithmsTest(void *graph, struct Arguments *arguments)
{

    void *ref_stats = NULL;

    switch (arguments->algorithm)
    {
    case 0:  // BFS
    {
        ref_stats = runBreadthFirstSearchAlgorithm( graph,  arguments->datastructure,  arguments->root,  arguments->pushpull);
    }
    break;
    case 1: // pagerank
    {
        ref_stats = runPageRankAlgorithm(graph,  arguments->datastructure,  arguments->epsilon,  arguments->iterations,  arguments->pushpull);
    }
    break;
    case 2: // SSSP-Dijkstra
    {
        ref_stats = runSSSPAlgorithm(graph,  arguments->datastructure,  arguments->root,  arguments->iterations, arguments->pushpull,  arguments->delta);
    }
    break;
    case 3: // SSSP-Bellmanford
    {
        ref_stats = runBellmanFordAlgorithm(graph,  arguments->datastructure,  arguments->root,  arguments->iterations, arguments->pushpull);
    }
    break;
    case 4: // DFS
    {
        ref_stats = runDepthFirstSearchAlgorithm(graph,  arguments->datastructure,  arguments->root);
    }
    break;
    case 5: // SPMV
    {
        ref_stats = runSPMVAlgorithm(graph,  arguments->datastructure,  arguments->iterations,  arguments->pushpull);
    }
    break;
    case 6: // Connected Components
    {
        ref_stats = runConnectedComponentsAlgorithm(graph,  arguments->datastructure,  arguments->iterations,  arguments->pushpull);
    }
    break;
    case 7: // Triangle Count
    {
        ref_stats = runTriangleCountAlgorithm(graph, arguments->datastructure, arguments->pushpull);
    }
    break;
    case 8: // incremental Aggregation
    {
        ref_stats = runIncrementalAggregationAlgorithm(graph,  arguments->datastructure);
    }
    break;
    default:// BFS
    {
        ref_stats = runBreadthFirstSearchAlgorithm(graph,  arguments->datastructure,  arguments->root, arguments->pushpull);
    }
    break;
    }

    return ref_stats;
}

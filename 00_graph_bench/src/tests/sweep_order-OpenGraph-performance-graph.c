// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : test_accel-graph.c
// Create : 2019-07-29 16:52:00
// Revise : 2019-09-28 15:36:29
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>
#include <stdbool.h>
#include <omp.h>
#include <assert.h>

#include "graphStats.h"
#include "edgeList.h"
#include "myMalloc.h"

#include "graphCSR.h"
#include "graphAdjLinkedList.h"
#include "graphAdjArrayList.h"
#include "graphGrid.h"

#include "mt19937.h"
#include "graphConfig.h"
#include "timer.h"
#include "graphRun.h"

#include "BFS.h"
#include "DFS.h"
#include "pageRank.h"
#include "incrementalAggregation.h"
#include "bellmanFord.h"
#include "SSSP.h"
#include "connectedComponents.h"
#include "triangleCount.h"

#include "graphTest.h"
// #define GRAPH_NUM 4
// #define GRAPH_NUM 14
#define GRAPH_NUM 9

#define THREAD_POINTS 1
#define THREAD_POINTS_START 5
#define THREAD_POINTS_TOTAL 9

#define ORDER_CONFIG 12
#define ORDER_CONFIG_TOTAL 12

void sweepSPMV(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepCC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepTC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepSSSP(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepSSSPDelta(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepBFS(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepPR(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);
void sweepBC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL]);

char graph_dir[1024];
char label_dir[1024];
char unified_perf_file[1024];

float PLRU_stats_BFS[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]            = {0};
float PLRU_stats_PR[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]             = {0};
float PLRU_stats_PR_DELTA[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]       = {0};
float PLRU_stats_BC[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]             = {0};
float PLRU_stats_SPMV[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]           = {0};
// float PLRU_stats_TC[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS]    = {0};
float PLRU_stats_CC[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]             = {0};
float PLRU_stats_CC_Afforest[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]    = {0};
float PLRU_stats_SSSP[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]           = {0};
float PLRU_stats_SSSP_DELTA[GRAPH_NUM][ORDER_CONFIG_TOTAL][THREAD_POINTS_TOTAL]     = {0};

uint32_t lmode_l2[ORDER_CONFIG_TOTAL] = {0, 4, 5, 11, 11, 11, 11, 11, 11, 11, 11, 11};
uint32_t lmode_l3[ORDER_CONFIG_TOTAL] = {0, 0, 0, 0, 4, 5, 0, 4, 5, 0, 4, 5};
uint32_t mmode[ORDER_CONFIG_TOTAL]    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
uint32_t threadPoints[THREAD_POINTS_TOTAL]    = {1, 2, 4, 8, 16, 32, 64, 160, 48};

int
main (int argc, char **argv)
{
    char *config_labels[ORDER_CONFIG_TOTAL] =
    {
        "Rand-Order",
        "DBG-out",
        "DBG-in",
        "Rabbit",
        "Rabbit+DBG-out",
        "Rabbit+DBG-in",
        "Gorder",
        "Gorder+DBG-out",
        "Gorder+DBG-in",
        "Corder",
        "Corder+DBG-out",
        "Corder+DBG-in"
    };

    char *reorder_labels[ORDER_CONFIG_TOTAL] =
    {
        "NO.labels",
        "NO.labels",
        "NO.labels",
        "graph_Rabbit.rand.labels",
        "graph_Rabbit.rand.labels",
        "graph_Rabbit.rand.labels",
        "graph_Gorder.rand.labels",
        "graph_Gorder.rand.labels",
        "graph_Gorder.rand.labels",
        "graph.rand.bin.txt.csr_Corder.labels",
        "graph.rand.bin.txt.csr_Corder.labels",
        "graph.rand.bin.txt.csr_Corder.labels"
    };

// -----------------------------------------------------------------------------


    // char *benchmarks_graphs[GRAPH_NUM] =
    // {
    //     "LAW-amazon-2008",
    //     "LAW-cnr-2000",
    //     "LAW-dblp-2010",
    //     "LAW-enron"
    // };

    // char *benchmarks_dir[GRAPH_NUM] =
    // {
    //     "../01_test_graphs/LAW/LAW-amazon-2008",
    //     "../01_test_graphs/LAW/LAW-cnr-2000",
    //     "../01_test_graphs/LAW/LAW-dblp-2010",
    //     "../01_test_graphs/LAW/LAW-enron"
    // };

// -----------------------------------------------------------------------------

    // char *benchmarks_graphs[GRAPH_NUM] =
    // {
    //     "GAP-road",
    //     "GAP-twitter",
    //     "GONG-gplus",
    //     "KONECT-wikipedia_link_en",
    //     "LAW-in-2004",
    //     "LAW-it-2004",
    //     "LAW-uk-2002",
    //     "LAW-uk-2005",
    //     "LAW-webbase-2001",
    //     "SNAP-cit-Patents",
    //     "SNAP-com-Orkut",
    //     "SNAP-soc-LiveJournal1",
    //     "SNAP-soc-Pokec",
    //     "SNAP-web-Google"
    // };

    // char *benchmarks_dir[GRAPH_NUM] =
    // {
    //     "../../01_GraphDatasets/GAP/GAP-road",
    //     "../../01_GraphDatasets/GAP/GAP-twitter",
    //     "../../01_GraphDatasets/GONG/GONG-gplus",
    //     "../../01_GraphDatasets/KONECT/KONECT-wikipedia_link_en",
    //     "../../01_GraphDatasets/LAW/LAW-in-2004",
    //     "../../01_GraphDatasets/LAW/LAW-it-2004",
    //     "../../01_GraphDatasets/LAW/LAW-uk-2002",
    //     "../../01_GraphDatasets/LAW/LAW-uk-2005",
    //     "../../01_GraphDatasets/LAW/LAW-webbase-2001",
    //     "../../01_GraphDatasets/SNAP/SNAP-cit-Patents",
    //     "../../01_GraphDatasets/SNAP/SNAP-com-Orkut",
    //     "../../01_GraphDatasets/SNAP/SNAP-soc-LiveJournal1",
    //     "../../01_GraphDatasets/SNAP/SNAP-soc-Pokec",
    //     "../../01_GraphDatasets/SNAP/SNAP-web-Google"
    // };

// -----------------------------------------------------------------------------

    char *benchmarks_graphs[GRAPH_NUM] =
    {
        "GAP-road",
        "GAP-twitter",
        "GONG-gplus",
        "KONECT-wikipedia_link_en",
        "LAW-uk-2005",
        "LAW-webbase-2001",
        "SNAP-cit-Patents",
        "SNAP-com-Orkut",
        "SNAP-soc-LiveJournal1"
    };

    char *benchmarks_dir[GRAPH_NUM] =
    {
        "../../01_GraphDatasets/GAP/GAP-road",
        "../../01_GraphDatasets/GAP/GAP-twitter",
        "../../01_GraphDatasets/GONG/GONG-gplus",
        "../../01_GraphDatasets/KONECT/KONECT-wikipedia_link_en",
        "../../01_GraphDatasets/LAW/LAW-uk-2005",
        "../../01_GraphDatasets/LAW/LAW-webbase-2001",
        "../../01_GraphDatasets/SNAP/SNAP-cit-Patents",
        "../../01_GraphDatasets/SNAP/SNAP-com-Orkut",
        "../../01_GraphDatasets/SNAP/SNAP-soc-LiveJournal1"
    };

// -----------------------------------------------------------------------------

    struct Arguments arguments;
    /* Default values. */

    arguments.wflag = 0;
    arguments.xflag = 0;
    arguments.sflag = 0;
    arguments.dflag = 1;

    arguments.iterations = 100;
    arguments.epsilon = 0.00000001;
    arguments.source = 5319;
    arguments.algorithm = 0; // BFS
    arguments.datastructure = 0;
    arguments.pushpull = 0;
    arguments.sort = 1;

    arguments.symmetric = 0;
    arguments.weighted = 0;
    arguments.delta = 1;

    arguments.pre_numThreads  = omp_get_max_threads();
    arguments.algo_numThreads = omp_get_max_threads();
    arguments.ker_numThreads = arguments.algo_numThreads ;

    arguments.fnameb = "../01_test_graphs/LAW/LAW-enron/graph.rand.bin";
    arguments.fnamel = "../01_test_graphs/LAW/LAW-enron/graph_Gorder.labels";
    arguments.fnameb_format = 1;
    arguments.convert_format = 1;
    arguments.trials = 4; // random number of trials
    initializeMersenneState (&(arguments.mt19937var), 27491095);
    omp_set_nested(1);

    arguments.lmode = 0;
    arguments.lmode_l2 = 0;
    arguments.lmode_l3 = 0;
    arguments.mmode = 0;

    void *graph = NULL;

    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t k = 0;
    // uint32_t kk = 0;
    // void *ref_data;

    sprintf(unified_perf_file, "%s/results_algo_all_time.unified.%s", "./openmp-results", "perf");

    for ( i = 0; i < GRAPH_NUM; ++i)
    {
        printf("graph %s\n", benchmarks_dir[i]);

        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            sprintf (graph_dir, "%s/%s", benchmarks_dir[i], "graph.rand.bin");
            sprintf (label_dir, "%s/%s", benchmarks_dir[i], reorder_labels[j]);
            arguments.lmode = 0; // base is random order
            arguments.lmode_l2 =  lmode_l2[j];
            arguments.lmode_l3 = lmode_l3[j];
            arguments.mmode = mmode[j];
            arguments.fnameb = graph_dir;
            arguments.fnamel = label_dir;
            printf("graph config %5u - %u %u %u - %s\n", j, arguments.lmode_l2, arguments.lmode_l3, arguments.mmode, config_labels[j]);

            graph = generateGraphDataStructure(&arguments);

            arguments.pushpull = 2; //
            sweepBFS(arguments, graph, &(PLRU_stats_BFS[i][j][0]));

            arguments.pushpull = 0; //
            sweepPR(arguments, graph, &(PLRU_stats_PR[i][j][0]));

            arguments.pushpull = 7; //
            sweepPR(arguments, graph, &(PLRU_stats_PR_DELTA[i][j][0]));

            arguments.pushpull = 0; //
            sweepBC(arguments, graph, &(PLRU_stats_BC[i][j][0]));

            arguments.pushpull = 0; //
            sweepSPMV(arguments, graph, &(PLRU_stats_SPMV[i][j][0]));

            // arguments.pushpull = 1; //
            // sweepTC(arguments, graph, &(PLRU_stats_TC[i][j][0]));

            arguments.pushpull = 0; //
            sweepCC(arguments, graph, &(PLRU_stats_CC[i][j][0]));

            arguments.pushpull = 1; //
            sweepCC(arguments, graph, &(PLRU_stats_CC_Afforest[i][j][0]));

            arguments.pushpull = 0; //
            sweepSSSP(arguments, graph, &(PLRU_stats_SSSP[i][j][0]));

            arguments.pushpull = 0; //
            sweepSSSPDelta(arguments, graph, &(PLRU_stats_SSSP_DELTA[i][j][0]));

            freeGraphDataStructure(graph, arguments.datastructure);
        }
    }
    // print out stats to file each graph processed
    FILE *fptr1;
    fptr1 = fopen(unified_perf_file, "a+");

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) BFS, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_BFS[i][j][k] / arguments.trials );
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) PR, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ", threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_PR[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) PR_DELTA, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_PR_DELTA[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) BC, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_BC[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) CC ShiloachVishkin, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_CC[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) CC Afforest, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_CC_Afforest[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }


    // for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    // {
    //     fprintf(fptr1, " -----------------------------------------------------\n");
    //     fprintf(fptr1, " Performance (Seconds) TC, Num Threads %u \n",   threadPoints[k]);
    //     fprintf(fptr1, " -----------------------------------------------------\n");

    //     fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
    //     for (j = 0; j < ORDER_CONFIG; ++j)
    //     {
    //         fprintf(fptr1, "%-14s, ",  config_labels[j]);
    //     }
    //     fprintf(fptr1, " \n");

    //     for ( i = 0; i < GRAPH_NUM; ++i)
    //     {
    //         fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
    //         for (j = 0; j < ORDER_CONFIG; ++j)
    //         {
    //             fprintf(fptr1, "%-14f, ",  PLRU_stats_TC[i][j][k]);
    //         }
    //         fprintf(fptr1, " \n");
    //     }
    //     fprintf(fptr1, " -----------------------------------------------------\n");
    // }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) SPMV, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_SPMV[i][j][k]);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) SSSP bellmanFord, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_SSSP[i][j][k] / arguments.trials);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }

    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        fprintf(fptr1, " -----------------------------------------------------\n");
        fprintf(fptr1, " Performance (Seconds) SSSP Delta, Num Threads %u \n",   threadPoints[k]);
        fprintf(fptr1, " -----------------------------------------------------\n");

        fprintf(fptr1, "NumThreads%-15u, ",  threadPoints[k]);
        for (j = 0; j < ORDER_CONFIG; ++j)
        {
            fprintf(fptr1, "%-14s, ",  config_labels[j]);
        }
        fprintf(fptr1, " \n");

        for ( i = 0; i < GRAPH_NUM; ++i)
        {
            fprintf(fptr1, "%-25s, ",  benchmarks_graphs[i]);
            for (j = 0; j < ORDER_CONFIG; ++j)
            {
                fprintf(fptr1, "%-14f, ",  PLRU_stats_SSSP_DELTA[i][j][k] / arguments.trials);
            }
            fprintf(fptr1, " \n");
        }
        fprintf(fptr1, " -----------------------------------------------------\n");
    }


    fclose(fptr1);

    exit (0);
}


void sweepBFS(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    uint32_t kk = 0;
    void *ref_data;
    arguments.algorithm = 0; // BFS
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        initializeMersenneState (&(arguments.mt19937var), 27491095);
        for(kk = 0 ; kk < arguments.trials; kk++)
        {
            arguments.source = generateRandomRootGeneral(&arguments, graph); // random root each trial
            ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
            PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
            // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
            freeGraphStatsGeneral(ref_data, arguments.algorithm);
        }
    }
}

void sweepSSSP(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    uint32_t kk = 0;
    void *ref_data;
    arguments.algorithm = 3; // SSSP-bellmanford
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        initializeMersenneState (&(arguments.mt19937var), 27491095);
        for(kk = 0 ; kk < arguments.trials; kk++)
        {
            arguments.source = generateRandomRootGeneral(&arguments, graph); // random root each trial
            ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
            PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
            // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
            freeGraphStatsGeneral(ref_data, arguments.algorithm);
        }
    }
}

void sweepSSSPDelta(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    uint32_t kk = 0;
    void *ref_data;
    arguments.algorithm = 2; // SSSP-delta
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        initializeMersenneState (&(arguments.mt19937var), 27491095);
        for(kk = 0 ; kk < arguments.trials; kk++)
        {
            arguments.source = generateRandomRootGeneral(&arguments, graph); // random root each trial
            ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
            PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
            // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
            freeGraphStatsGeneral(ref_data, arguments.algorithm);
        }
    }
}

void sweepPR(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    void *ref_data;
    arguments.algorithm = 1; // PR
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
        PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
        // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
        freeGraphStatsGeneral(ref_data, arguments.algorithm);
    }

}

void sweepBC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    void *ref_data;
    arguments.algorithm = 7; // BC
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
        PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
        // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
        freeGraphStatsGeneral(ref_data, arguments.algorithm);
    }

}

void sweepTC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    // case 3: // With binary intersection
    //    stats = triangleCountBinaryIntersectionGraphCSR(graph);
    //    break;}
    uint32_t k = 0;
    void *ref_data;
    arguments.algorithm = 8; // TC
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
        PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
        // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
        freeGraphStatsGeneral(ref_data, arguments.algorithm);
    }

}


void sweepSPMV(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    void *ref_data;
    arguments.algorithm = 5; // SPMV
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
        PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
        // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
        freeGraphStatsGeneral(ref_data, arguments.algorithm);
    }

}

void sweepCC(struct Arguments arguments, void *graph, float PLRU_stats[THREAD_POINTS_TOTAL])
{
    uint32_t k = 0;
    void *ref_data;
    arguments.algorithm = 6; // CC
    for(k = THREAD_POINTS_START; k < THREAD_POINTS+THREAD_POINTS_START; k ++)
    {
        arguments.algo_numThreads = threadPoints[k];
        arguments.ker_numThreads = arguments.algo_numThreads ;
        ref_data = runGraphAlgorithmsTest(&arguments, graph); // ref stats should mach oother algo
        PLRU_stats[k] += getGraphAlgorithmsTestTime(ref_data, arguments.algorithm);
        // printStatsDoubleTaggedCacheToFile(ref_stats_tmp->cache, unified_perf_file);
        freeGraphStatsGeneral(ref_data, arguments.algorithm);
    }

}
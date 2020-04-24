// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Mohannad Ibrahim
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
#include <string.h>
#include <math.h>
#include <stdint.h>


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

#include <assert.h>
#include "graphTest.h"

int numThreads;
mt19937state *mt19937var;
#define GRAPH_NUM 26

// "   mm                        ""#             mmm                       #     \n"
// "   ##    mmm    mmm    mmm     #           m"   "  m mm   mmm   mmmm   # mm  \n"
// "  #  #  #"  "  #"  "  #"  #    #           #   mm  #"  " "   #  #" "#  #"  # \n"
// "  #mm#  #      #      #""""    #     """   #    #  #     m"""#  #   #  #   # \n"
// " #    # "#mm"  "#mm"  "#mm"    "mm          "mmm"  #     "mm"#  ##m#"  #   # \n"
// "                                                                #            \n"


int main (int argc, char **argv)
{

    uint32_t benchmark_count = GRAPH_NUM;
    char *benchmarks_law[GRAPH_NUM] =
    {
        "../../01_GraphDatasets/amazon-2008/graph.wbin",
        "../../01_GraphDatasets/arabic-2005/graph.wbin",
        "../../01_GraphDatasets/cnr-2000/graph.wbin",
        "../../01_GraphDatasets/dblp-2010/graph.wbin",
        "../../01_GraphDatasets/enron/graph.wbin",
        "../../01_GraphDatasets/eu-2005/graph.wbin",
        "../../01_GraphDatasets/hollywood-2009/graph.wbin",
        "../../01_GraphDatasets/in-2004/graph.wbin",
        "../../01_GraphDatasets/indochina-2004/graph.wbin",
        "../../01_GraphDatasets/it-2004/graph.wbin",
        "../../01_GraphDatasets/ljournal-2008/graph.wbin",
        "../../01_GraphDatasets/sk-2005/graph.wbin",
        "../../01_GraphDatasets/uk-2002/graph.wbin",
        "../../01_GraphDatasets/uk-2005/graph.wbin",
        "../../01_GraphDatasets/webbase-2001/graph.wbin",

        "../../01_GraphDatasets/Gong-gplus/graph.wbin",
        "../../01_GraphDatasets/GAP-road/graph.wbin",
        "../../01_GraphDatasets/SNAP-soc-pokec/graph.wbin",
        "../../01_GraphDatasets/SNAP-cit-Patents/graph.wbin",
        "../../01_GraphDatasets/SNAP-com-orkut/graph.wbin",
        "../../01_GraphDatasets/SNAP-soc-LiveJournal1/graph.wbin",
        "../../01_GraphDatasets/KONECT-wikipedia_link_en/graph.wbin",

        "../../01_GraphDatasets/gplus/graph.wbin",
        "../../01_GraphDatasets/USA-road/graph.wbin",
        "../../01_GraphDatasets/enwiki-2013/graph.wbin",
        "../../01_GraphDatasets/twitter/graph.wbin"
        
    };

    // char *benchmarks_test[7] =
    // {
    //     "../04_test_graphs/test/graph.wbin",
    //     "../04_test_graphs/v51_e1021/graph.wbin",
    //     "../04_test_graphs/v300_e2730/graph.wbin",
    //     "../04_test_graphs/amazon/graph.wbin",
    //     "../04_test_graphs/dblp/graph.wbin",
    //     "../04_test_graphs/euall/graph.wbin",
    //     "../04_test_graphs/Gnutella/graph.wbin"
    // };

    uint32_t pushpull_size = 8;
    uint32_t pushpull[] = {4, 2, 4, 9, 10, 11, 12, 13} ;

    struct Arguments arguments;
    arguments.wflag = 0;
    arguments.xflag = 0;
    arguments.sflag = 0;
    arguments.Sflag = 0;

    arguments.iterations = 300;
    arguments.trials = 1;
    arguments.epsilon = 1e-8;
    arguments.root = 5319;
    arguments.algorithm = 1;
    arguments.datastructure = 0;
    arguments.pushpull = 0;
    arguments.sort = 0;
    arguments.lmode = 0;
    arguments.symmetric = 0;
    arguments.weighted = 0;
    arguments.delta = 1;
    arguments.numThreads = omp_get_max_threads();
    arguments.fnameb_format = 1;
    arguments.convert_format = 1;
    arguments.verbosity = 0;
    arguments.binSize = 100;

    void *graph = NULL;

    numThreads =  arguments.numThreads;

    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));

    mt19937var = (mt19937state *) my_malloc(sizeof(mt19937state));
    initializeMersenneState (mt19937var, 27491095);

    omp_set_nested(1);
    omp_set_num_threads(numThreads);

    printf("*-----------------------------------------------------*\n");
    printf("| %-20s %-30u | \n", "Number of Threads :", numThreads);
    printf(" -----------------------------------------------------\n");

    void *ref_data;
    void *cmp_data;

    //for every benchmark
    uint32_t i;
    uint32_t j;

    Start(timer);
    for(i = 0; i < benchmark_count; i++)
    {
        arguments.fnameb = benchmarks_law[i];
        arguments.pushpull = 0;

        // arguments.fnameb = benchmarks_test[i];

        //appropriate filename
        printf("Begin tests for %s\n", arguments.fnameb);

        graph = generateGraphDataStructure(&arguments);

        if(graph == NULL) continue;

        ref_data = runGraphAlgorithmsTest(graph, &arguments); // ref stats should mach oother algo

        if(ref_data == NULL) continue;

        for(j = 0 ;  j < pushpull_size; j++)
        {
            arguments.pushpull = pushpull[j];
            arguments.verbosity = 0;

            cmp_data = runGraphAlgorithmsTest(graph, &arguments);

            struct PageRankStats *ref_stats_tmp = (struct PageRankStats * )ref_data;
            struct PageRankStats *cmp_stats_tmp = (struct PageRankStats * )cmp_data;

            // arguments.verbosity = 3; // will calculate avg for whole graph a bit slower/ and dump the topK ranks in file ref vs calc
            arguments.verbosity = 2; // will calculate avg for whole graph a bit slower
            // arguments.verbosity = 1;
            collectStatsPageRank(&arguments, ref_stats_tmp, cmp_stats_tmp, 1);
            arguments.verbosity = 0;

            freeGraphStatsGeneral(cmp_data, arguments.algorithm);
        }

        freeGraphStatsGeneral(ref_data, arguments.algorithm);
        freeGraphDataStructure(graph, arguments.datastructure);

        printf("Finished tests for %s\n Time\n", arguments.fnameb);
    }

    Stop(timer);
    printf("Page Rank Error Test Done ....... Time (%-9f)\n", Seconds(timer));
    free(timer);

    return 0;
}






// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : accel-graph.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:28
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>
#include <stdbool.h>
#include <omp.h>
#include <stdint.h>

#include "myMalloc.h"
#include "timer.h"
#include "mt19937.h"
#include "graphConfig.h"
#include "graphRun.h"
#include "graphStats.h"
#include "edgeList.h"


// "   mm                        ""#             mmm                       #     \n"
// "   ##    mmm    mmm    mmm     #           m"   "  m mm   mmm   mmmm   # mm  \n"
// "  #  #  #"  "  #"  "  #"  #    #           #   mm  #"  " "   #  #" "#  #"  # \n"
// "  #mm#  #      #      #""""    #     """   #    #  #     m"""#  #   #  #   # \n"
// " #    # "#mm"  "#mm"  "#mm"    "mm          "mmm"  #     "mm"#  ##m#"  #   # \n"
// "                                                                #            \n"

uint64_t afu_config;
uint64_t cu_config;
uint64_t afu_config_2;
uint64_t cu_config_2;

int numThreads;
mt19937state *mt19937var;

const char *argp_program_version =
    "AccelGraph_CAPI v3.0";
const char *argp_program_bug_address =
    "<atmughra@ncsu.edu>|<atmughrabi@gmail.com>";
/* Program documentation. */
static char doc[] =
    "\nAccelGraph is an open source graph processing framework, it is designed to be a benchmarking suite for various graph processing algorithms on FPGAs.\n";

/* A description of the arguments we accept. */
static char args_doc[] = "-f <graph file> -d [data structure] -a [algorithm] -r [root] -n [num threads] [-h -c -s -w]";

/* The options we understand. */
static struct argp_option options[] =
{
    {
        "graph-file",         'f', "<FILE>\n",      0,
        "\nEdge list represents the graph binary format to run the algorithm textual format change graph-file-format.\n"
    },
    {
        "graph-file-format",  'z', "[TEXT|BIN|CSR:1]\n",      0,
        "\nSpecify file format to be read, is it textual edge list, or a binary file edge list. This is specifically useful if you have Graph CSR/Grid structure already saved in a binary file format to skip the preprocessing step. [0]-text edgeList [1]-binary edgeList [2]-graphCSR binary.\n"
    },
    {
        "algorithm",         'a', "[DEFAULT:0]\n",      0,
        "\n[0]-BFS, [1]-Page-rank, [2]-SSSP-DeltaStepping, [3]-SSSP-BellmanFord, [4]-DFS,[5]-SPMV, [6]-Connected-Components, [7]-Triangle Counting, [8]-IncrementalAggregation.\n"
    },
    {
        "data-structure",    'd', "[DEFAULT:0]\n",      0,
        "\n[0]-CSR, [1]-Grid, [2]-Adj LinkedList, [3]-Adj ArrayList [4-5] same order bitmap frontiers.\n"
    },
    {
        "root",              'r', "[DEFAULT:0]\n",      0,
        "\nBFS, DFS, SSSP root"
    },
    {
        "direction",         'p', "[DEFAULT:0]\n",      0,
        "\n[0]-PULL, [1]-PUSH,[2]-HYBRID. NOTE: Please consult the function switch table for each algorithm.\n"
    },
    {
        "sort",              'o', "[DEFAULT:0]\n",      0,
        "\n[0]-radix-src [1]-radix-src-dest [2]-count-src [3]-count-src-dst.\n"
    },
    {
        "num-threads",       'n', "[DEFAULT:MAX]\n",      0,
        "\nDefault:max number of threads the system has"
    },
    {
        "num-iterations",    'i', "[DEFAULT:20]\n",      0,
        "\nNumber of iterations for page rank to converge [default:20] SSSP-BellmanFord [default:V-1].\n"
    },
    {
        "num-trials",        't', "[DEFAULT:1]\n",      0,
        "\nNumber of trials for whole run (graph algorithm run) [default:0].\n"
    },
    {
        "tolerance",         'e', "[EPSILON:0.0001]\n",      0,
        "\nTolerance value of for page rank [default:0.0001].\n"
    },
    {
        "delta",             'b', "[DELTA:1]\n",      0,
        "\nSSSP Delta value [Default:1].\n"
    },
    {
        "light-reorder",     'l', "[ORDER:0]\n",      0,
        "\nRelabels the graph for better cache performance. [default:0]-no-reordering [1]-page-rank-order [2]-in-degree [3]-out-degree [4]-in/out degree [5]-Rabbit [6]-Epoch-pageRank [7]-Epoch-BFS [8]-LoadFromFile\n"
    },
    {
        "convert-format",    'c', "[TEXT|BIN|CSR:1]\n",      0,
        "\n[serialize flag must be on --serialize to write] Serialize graph text format (edge list format) to binary graph file on load example:-f <graph file> -c this is specifically useful if you have Graph CSR/Grid structure and want to save in a binary file format to skip the preprocessing step for future runs. [0]-text edgeList [1]-binary edgeList [2]-graphCSR binary.\n"
    },
    {
        "generate-weights",  'w', 0,      0,
        "\nGenerate random weights don't load from graph file. Check ->graphConfig.h #define WEIGHTED 1 beforehand then recompile using this option.\n"
    },
    {
        "symmetrize",        's', 0,      0,
        "\nSymmetric graph, create a set of incoming edges.\n"
    },
    {
        "serialize",         'x', 0,      0,
        "\nEnable file conversion/serialization use with --convert-format.\n"
    },
    {
        "stats",             'S', 0,      0,
        "\nWrite algorithm stats to file. same directory as the graph.\nPageRank: Dumps top-k ranks matching using QPR similarity metrics.\n"
    },
    {
        "bin-size",         'g', "[SIZE:512]\n",      0,
        "\nYou bin vertices's histogram according to this parameter, if you have a large graph you want to illustrate.\n"
    },
    {
        "verbosity",    'j', "[DEFAULT:0]\n",      0,
        "\nFor now it controls the output of .perf file and PageRank .stats (needs --stats enabled) files\nPageRank .stat [1:top-k results] [2:top-k results and top-k ranked vertices listed.\n"
    },
    {
        "remove-duplicate", 'k', 0,      0,
        "\nRemovers duplicate edges and self loops from the graph.\n"
    },
    {
        "afu-config",            'm', "[DEFAULT:0x1]\n",      0,
        "\nAFU-Control buffers(read/write/prefetcher) arbitration 0x01 round robin 0x10 fixed priority.\n"
    },
    {
        "cu-config",             'q', "[DEFAULT:0x01]\n",      0,
        "\nCU configurations for requests cached/non cached/prefetcher active or not check README for more explanation.\n"
    },
    { 0 }
};



/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    char *eptr;

    switch (key)
    {
    case 'f':
        arguments->fnameb = arg;
        break;
    case 'z':
        arguments->fnameb_format = atoi(arg);
        break;
    case 'd':
        arguments->datastructure = atoi(arg);
        break;
    case 'a':
        arguments->algorithm = atoi(arg);
        break;
    case 'r':
        arguments->root = atoi(arg);
        break;
    case 'n':
        arguments->numThreads = atoi(arg);
        break;
    case 'i':
        arguments->iterations = atoi(arg);
        break;
    case 't':
        arguments->trials = atoi(arg);
        break;
    case 'e':
        arguments->epsilon = atof(arg);
        break;
    case 'p':
        arguments->pushpull = atoi(arg);
        break;
    case 'o':
        arguments->sort = atoi(arg);
        break;
    case 'l':
        arguments->lmode = atoi(arg);
        break;
    case 'b':
        arguments->delta = atoi(arg);
        break;
    case 's':
        arguments->symmetric = 1;
        break;
    case 'S':
        arguments->Sflag = 1;
        break;
    case 'w':
        arguments->weighted = 1;
        break;
    case 'x':
        arguments->xflag = 1;
        break;
    case 'c':
        arguments->convert_format = atoi(arg);
        break;
    case 'g':
        arguments->binSize = atoi(arg);
        break;
    case 'j':
        arguments->verbosity = atoi(arg);
        break;
    case 'k':
        arguments->dflag = 1;
        break;
    case 'm':
        arguments->afu_config = strtoll(arg, &eptr, 0);
        break;
    case 'q':
        arguments->cu_config = strtoll(arg, &eptr, 0);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


static struct argp argp = { options, parse_opt, args_doc, doc };

int
main (int argc, char **argv)
{

    struct Arguments arguments;
    /* Default values. */

    arguments.wflag = 0;
    arguments.xflag = 0;
    arguments.sflag = 0;
    arguments.Sflag = 0;
    arguments.dflag = 0;
    arguments.binSize = 512;
    arguments.verbosity = 0;
    arguments.iterations = 20;
    arguments.trials = 1;
    arguments.epsilon = 0.0001;
    arguments.root = 0;
    arguments.algorithm = 0;
    arguments.datastructure = 0;
    arguments.pushpull = 0;
    arguments.sort = 0;
    arguments.lmode = 0;
    arguments.symmetric = 0;
    arguments.weighted = 0;
    arguments.delta = 1;
    arguments.numThreads = omp_get_max_threads();
    arguments.fnameb = NULL;
    arguments.fnameb_format = 1;
    arguments.convert_format = 1;
    arguments.afu_config = 0x01;
    arguments.cu_config  = 0x01;
    arguments.afu_config_2 = 0x01;
    arguments.cu_config_2  = 0x01;

    void *graph = NULL;

    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    numThreads =  omp_get_max_threads();
    afu_config =  arguments.afu_config;
    cu_config  =  arguments.cu_config;
    afu_config_2  =  arguments.afu_config_2;
    cu_config_2   =  arguments.cu_config_2;

    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));

    mt19937var = (mt19937state *) my_malloc(sizeof(mt19937state));
    initializeMersenneState (mt19937var, 27491095);

    omp_set_nested(1);
    omp_set_num_threads(numThreads);



    printf("*-----------------------------------------------------*\n");
    printf("| %-25s %-25d | \n", "Number of Threads Pre :", numThreads);
    printf(" -----------------------------------------------------\n");

    if(arguments.xflag) // if stats flag is on collect stats or serialize your graph
    {
        writeSerializedGraphDataStructure(&arguments);
    }
    else
    {

        graph = generateGraphDataStructure(&arguments);


        numThreads =  arguments.numThreads;
        omp_set_num_threads(numThreads);

        printf("*-----------------------------------------------------*\n");
        printf("| %-25s %-25d | \n", "Number of Threads Algo :", numThreads);
        printf(" -----------------------------------------------------\n");


        runGraphAlgorithms(graph, &arguments);
        freeGraphDataStructure(graph, arguments.datastructure);
    }




    free(timer);
    exit (0);
}






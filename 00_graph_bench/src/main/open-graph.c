// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : open-graph.c
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


int numThreads;
mt19937state *mt19937var;

const char *argp_program_version =
    "OpenGraph 1.0";
const char *argp_program_bug_address =
    "<atmughrabi@gmail.com> <atmughra@ncsu.edu>";
/* Program documentation. */
static char doc[] =
    "OpenGraph is an open source graph processing framework, it is designed to be a portable benchmarking suite for various graph processing algorithms.";

/* A description of the arguments we accept. */
static char args_doc[] = "-f <graph file> -d [data structure] -a [algorithm] -r [root] -n [num threads] [-h -c -s -w]";

/* The options we understand. */
static struct argp_option options[] =
{
    {
        "graph-file",         'f', "<FILE>",      0,
        "\nEdgelist representing the graph,in binary format to run the algorithm. Or if in textual format change graph-file-format"
    },
    {
        "graph-file-format",  'z', "[TEXT|BIN|CSR:1]",      0,
        "\nSpecify file format to be read, is it textual edge list, or a binary file edge list. This is specifically useful if you have Graph CSR/Grid structure already saved in a binary file format to skip the preprocessing step. [0]-text edgeList [1]-binary edgeList [2]-graphCSR binary"
    },
    {
        "algorithm",         'a', "[DEFAULT:0]",      0,
        "\n[0]-BFS, [1]-Page-rank, [2]-SSSP-DeltaStepping, [3]-SSSP-BellmanFord, [4]-DFS,[5]-SPMV, [6]-Connected-Components, [7]-Triangle Counting, [8]-IncrementalAggregation"
    },
    {
        "data-structure",    'd', "[DEFAULT:0]",      0,
        "\n[0]-CSR, [1]-Grid, [2]-Adj LinkedList, [3]-Adj ArrayList [4-5] same order bitmap frontiers"
    },
    {
        "root",              'r', "[DEFAULT:0]",      0,
        "\nBFS, DFS, SSSP root"
    },
    {
        "direction",         'p', "[DEFAULT:0]",      0,
        "\n[0]-PULL, [1]-PUSH,[2]-HYBRID. NOTE: Please consult the function switch table for each algorithm"
    },
    {
        "sort",              'o', "[DEFAULT:0]",      0,
        "\n[0]-radix-src [1]-radix-src-dest [2]-count-src [3]-count-src-dst"
    },
    {
        "num-threads",       'n', "[DEFAULT:MAX]",      0,
        "\nDefault:max number of threads the system has"
    },
    {
        "num-iterations",    'i', "[DEFAULT:20]",      0,
        "\nNumber of iterations for page rank to converge [default:20] SSSP-BellmanFord [default:V-1] "
    },
    {
        "num-trials",        't', "[DEFAULT:1]",      0,
        "\nNumber of random trials for each whole run (graph algorithm run) [default:0] "
    },
    {
        "tolerance",         'e', "[EPSILON:0.0001]",      0,
        "\nTolerance value of for page rank [default:0.0001] "
    },
    {
        "epsilon",           'e', "[DEFAULT:0.0001]",      OPTION_ALIAS
    },
    {
        "delta",             'b', "[DELTA:1]",      0,
        "\nSSSP Delta value [Default:1]"
    },
    {
        "light-reorder",     'l', "[ORDER:0]",      0,
        "\nRelabels the graph for better cache performance. [default:0]-no-reordering [1]-page-rank-order [2]-in-degree [3]-out-degree [4]-in/out degree [5]-Rabbit [6]-Epoch-pageRank [7]-Epoch-BFS [8]-LoadFromFile "
    },
    {
        "convert-format",    'c', "[TEXT|BIN|CSR:1]",      0,
        "\n[stats flag must be on --stats to write]Serialize graph text format (edge list format) to binary graph file on load example:-f <graph file> -c this is specifically useful if you have Graph CSR/Grid structure and want to save in a binary file format to skip the preprocessing step for future runs. [0]-text edgeList [1]-binary edgeList [2]-graphCSR binary"
    },
    {
        "generate-weights",  'w', 0,      0,
        "\nGenerate random weights don't load from graph file. Check ->graphConfig.h #define WEIGHTED 1 beforehand then recompile using this option"
    },
    {
        "symmetrize",        's', 0,      0,
        "\nSymmetrize graph, create a set of incoming edges"
    },
    {
        "stats",             'x', 0,      0,
        "\nDump a histogram to file based on in-out degree count bins / sorted according to in/out-degree or page-ranks\n--bin-size/--in-out-degree are related bin verticies in  terms of in or out degrees "
    },
    {
        "bin-size",         'g', "SIZE:512",      0,
        "\nYou bin vertices's histogram according to this parameter, if you have a large graph you want to illustrate "
    },
    {
        "in-out-degree",    'j', "[DEFAULT:2]",      0,
        "\n[1]-in-degree, [2]-out-degree, bin histogram with out/in-degree binned. [DEFAULT:2]  "
    },
    {
        "remove-duplicate", 'k', 0,      0,
        "\nRemovers duplicate edges and self loops from the graph"
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
        arguments->inout_degree = atoi(arg);
        break;
    case 'k':
        arguments->dflag = 1;
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
    arguments.dflag = 0;
    arguments.binSize = 512;
    arguments.inout_degree = 0;
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

    void *graph = NULL;

    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    numThreads =  arguments.numThreads;

    struct Timer *timer = (struct Timer *) my_malloc(sizeof(struct Timer));

    mt19937var = (mt19937state *) my_malloc(sizeof(mt19937state));
    initializeMersenneState (mt19937var, 27491095);

    omp_set_nested(1);
    omp_set_num_threads(numThreads);




    printf("*-----------------------------------------------------*\n");
    printf("| %-20s %-30u | \n", "Number of Threads :", numThreads);
    printf(" -----------------------------------------------------\n");

    if(arguments.xflag) // if stats flag is on collect stats or serialize your graph
    {
        // uint32_t binSize = arguments.iterations;
        // uint32_t inout_degree = arguments.pushpull;
        // uint32_t inout_lmode = arguments.lmode;
        // collectStats(binSize, arguments.fnameb, arguments.sort, inout_lmode, arguments.symmetric, arguments.weighted, inout_degree);
        writeSerializedGraphDataStructure(&arguments);
    }
    else
    {

        graph = generateGraphDataStructure(&arguments);
        runGraphAlgorithms(graph, &arguments);
        freeGraphDataStructure(graph, arguments.datastructure);
    }




    free(timer);
    exit (0);
}






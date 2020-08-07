// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphRun.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "mt19937.h"
#include "graphConfig.h"
#include "timer.h"

#include "graphCSR.h"
#include "graphGrid.h"
#include "graphAdjLinkedList.h"
#include "graphAdjArrayList.h"

#include "BFS.h"
#include "DFS.h"
#include "pageRank.h"
#include "incrementalAggregation.h"
#include "bellmanFord.h"
#include "SSSP.h"
#include "SPMV.h"
#include "connectedComponents.h"
#include "triangleCount.h"
#include "betweennessCentrality.h"


#include "reorder.h"
#include "graphStats.h"
#include "graphRun.h"


void generateGraphPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51lf | \n", time);
    printf(" -----------------------------------------------------\n");

}



void writeSerializedGraphDataStructure(struct Arguments *arguments)  // for now this only support graph CSR
{

    // check input type edgelist text/bin or graph csr
    // read input file create CSR graph then write to binaryfile
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    if(arguments->fnameb_format == 0 && arguments->convert_format == 1)  // for now it edge list is text only convert to binary
    {
        Start(timer);
        arguments->fnameb = readEdgeListstxt(arguments->fnameb, arguments->weighted);
        arguments->fnameb_format = 1; // now you have a bin file
#if WEIGHTED
        arguments->weighted = 1; // no need to generate weights again this affects readedgelistbin
#else
        arguments->weighted = 0;
#endif
        Stop(timer);
        generateGraphPrintMessageWithtime("Serialize EdgeList text to binary (Seconds)", Seconds(timer));
    }
    else if(arguments->fnameb_format == 1 && arguments->convert_format == 0)  // for now it edge list is text only convert to binary
    {
        Start(timer);
        struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted);  // read edglist from binary file

        if(arguments->lmode)
            edgeList = reorderGraphProcess(edgeList, arguments);

        writeEdgeListToTXTFile(edgeList, arguments->fnameb);
        arguments->fnameb_format = 1; // now you have a bin file
#if WEIGHTED
        arguments->weighted = 1; // no need to generate weights again this affects readedgelistbin
#else
        arguments->weighted = 0;
#endif
        Stop(timer);
        generateGraphPrintMessageWithtime("Serialize EdgeList binary to text (Seconds)", Seconds(timer));

        freeEdgeList(edgeList);
    }
    else if(arguments->fnameb_format == 0 && arguments->convert_format == 2)  // for now it edge list is text only convert to binary
    {
        void *graph = NULL;
        struct GraphCSR *graphCSR = NULL;

        Start(timer);
        arguments->fnameb = readEdgeListstxt(arguments->fnameb, arguments->weighted);
        arguments->fnameb_format = 1; // now you have a bin file
#if WEIGHTED
        arguments->weighted = 1; // no need to generate weights again this affects readedgelistbin
#else
        arguments->weighted = 0;
#endif
        Stop(timer);
        generateGraphPrintMessageWithtime("Serialize EdgeList text to binary (Seconds)", Seconds(timer));

        Start(timer);
        graph = (void *)graphCSRPreProcessingStep (arguments);
        Stop(timer);
        generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));

        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        writeToBinFileGraphCSR (arguments->fnameb, graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));

        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
    }
    else if(arguments->fnameb_format == 1 && arguments->convert_format == 2)   // for now it edge list is text only convert to binary
    {
        void *graph = NULL;
        struct GraphCSR *graphCSR = NULL;


        Start(timer);
        graph = (void *)graphCSRPreProcessingStep (arguments);
        Stop(timer);
        generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));


        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        writeToBinFileGraphCSR (arguments->fnameb, graph);
        Stop(timer);
        generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));

        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
    }
    else if(arguments->fnameb_format == 0 && arguments->convert_format == 0)
    {
        Start(timer);
        arguments->fnameb = readEdgeListstxt(arguments->fnameb, arguments->weighted);
        arguments->fnameb_format = 1; // now you have a bin file
#if WEIGHTED
        arguments->weighted = 1; // no need to generate weights again this affects readedgelistbin
#else
        arguments->weighted = 0;
#endif
        struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted);

        if(arguments->lmode)
            edgeList = reorderGraphProcess(edgeList, arguments);

        writeEdgeListToTXTFile(edgeList, arguments->fnameb);
        Stop(timer);

        generateGraphPrintMessageWithtime("INPUT and OUTPUT Same format", Seconds(timer));


    }
    else if(arguments->fnameb_format == arguments->convert_format)    // for now it edge list is text only convert to binary
    {

        Start(timer);



        Stop(timer);

        generateGraphPrintMessageWithtime("INPUT and OUTPUT Same format no need to serialize", Seconds(timer));
    }


    free(timer);

}

void readSerializeGraphDataStructure(struct Arguments *arguments)  // for now this only support graph CSR
{

    // check input type edgelist text/bin or graph csr
    // read input file create to the correct structure without preprocessing

}


void *generateGraphDataStructure(struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    void *graph = NULL;

    if(arguments->algorithm == 8)  // Triangle counting depends on order
    {

        arguments->sort = 1;
        // arguments->lmode = 2;
    }

    if(arguments->algorithm == 9)  // Incremental aggregation order
    {

        arguments->sort = 1;
        // arguments->lmode = 2;
    }

    if(arguments->fnameb_format == 0)  // for now it edge list is text only convert to binary
    {
        Start(timer);
        arguments->fnameb = readEdgeListstxt(arguments->fnameb, arguments->weighted);
        arguments->fnameb_format = 1; // now you have a bin file
#if WEIGHTED
        arguments->weighted = 1; // no need to generate weights again this affects readedgelistbin
#else
        arguments->weighted = 0;
#endif
        Stop(timer);
        generateGraphPrintMessageWithtime("Serialize EdgeList text to binary (Seconds)", Seconds(timer));
    }

    if(arguments->fnameb_format == 1 ) // if it is a graphCSR binary file
    {

        switch (arguments->datastructure)
        {
        case 0: // CSR
        case 4:
            Start(timer);
            graph = (void *)graphCSRPreProcessingStep (arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        case 1: // Grid
        case 5:
            Start(timer);
            graph = (void *)graphGridPreProcessingStep (arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphGrid Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        case 2: // Adj Linked List
            Start(timer);
            graph = (void *)graphAdjLinkedListPreProcessingStep (arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphAdjLinkedList Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        case 3: // Adj Array List
            Start(timer);
            graph = (void *)graphAdjArrayListPreProcessingStep ( arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphAdjArrayList Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        default:// CSR
            Start(timer);
            graph = (void *)graphCSRPreProcessingStep ( arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));

            break;
        }
    }
    else if(arguments->fnameb_format == 2)
    {
        Start(timer);
        graph = (void *)readFromBinFileGraphCSR (arguments->fnameb);
        Stop(timer);
        generateGraphPrintMessageWithtime("GraphCSR Preprocessing Step Time (Seconds)", Seconds(timer));
    }
    else
    {
        Start(timer);
        Stop(timer);
        generateGraphPrintMessageWithtime("UNKOWN Graph format Preprocessing Step Time (Seconds)", Seconds(timer));
    }

    free(timer);
    return graph;

}


void runGraphAlgorithms(void *graph, struct Arguments *arguments)
{

    // print total average stats to an external file fnameb.stats numthreads avg trial time

    double time_total = 0.0f;
    uint32_t  trials = arguments->trials;


    while(trials)
    {
        switch (arguments->algorithm)
        {
        case 0:  // BFS
        {
            struct BFSStats *stats = runBreadthFirstSearchAlgorithm( graph,  arguments->datastructure,  arguments->root,  arguments->pushpull);
            time_total += stats->time_total;
            freeBFSStats(stats);
        }
        break;
        case 1: // pagerank
        {
            struct PageRankStats *stats = runPageRankAlgorithm(graph,  arguments->datastructure,  arguments->epsilon,  arguments->iterations,  arguments->pushpull);
            time_total += stats->time_total;

            if(arguments->Sflag) // output page rank error statistics
            {
                struct PageRankStats *ref_stats = runPageRankAlgorithm(graph,  arguments->datastructure,  arguments->epsilon,  arguments->iterations,  0);
                collectStatsPageRank(arguments, ref_stats, stats, trials);
                freePageRankStats(ref_stats);
            }

            freePageRankStats(stats);
        }
        break;
        case 2: // SSSP-Delta
        {
            struct SSSPStats *stats = runSSSPAlgorithm(graph,  arguments->datastructure,  arguments->root,  arguments->iterations, arguments->pushpull,  arguments->delta);
            time_total += stats->time_total;
            freeSSSPStats(stats);
        }
        break;
        case 3: // SSSP-Bellmanford
        {
            struct BellmanFordStats *stats = runBellmanFordAlgorithm(graph,  arguments->datastructure,  arguments->root,  arguments->iterations, arguments->pushpull);
            time_total += stats->time_total;
            freeBellmanFordStats(stats);
        }
        break;
        case 4: // DFS
        {
            struct DFSStats *stats = runDepthFirstSearchAlgorithm(graph,  arguments->datastructure,  arguments->root);
            time_total += stats->time_total;
            freeDFSStats(stats);
        }
        break;
        case 5: // SPMV
        {
            struct SPMVStats *stats = runSPMVAlgorithm(graph,  arguments->datastructure,  arguments->iterations,  arguments->pushpull);
            time_total += stats->time_total;
            freeSPMVStats(stats);
        }
        break;
        case 6: // Connected Components
        {
            struct CCStats *stats = runConnectedComponentsAlgorithm(graph,  arguments->datastructure,  arguments->iterations,  arguments->pushpull);
            time_total += stats->time_total;
            freeCCStats(stats);
        }
        break;
        case 7: // Betweenness Centrality
        {
            struct BetweennessCentralityStats *stats = runBetweennessCentralityAlgorithm(graph,  arguments->datastructure, arguments->iterations, arguments->pushpull);
            time_total += stats->time_total;
            freeBetweennessCentralityStats(stats);
        }
        break;
        case 8: // Triangle Counting
        {
            struct TCStats *stats = runTriangleCountAlgorithm(graph, arguments->datastructure, arguments->pushpull);
            time_total += stats->time_total;
            freeTCStats(stats);
        }
        break;
        case 9: // incremental Aggregation
        {
            struct IncrementalAggregationStats *stats = runIncrementalAggregationAlgorithm(graph,  arguments->datastructure);
            time_total += stats->time_total;
            freeIncrementalAggregationStats(stats);
        }
        break;

        default: // BFS
        {
            struct BFSStats *stats = runBreadthFirstSearchAlgorithm(graph,  arguments->datastructure,  arguments->root, arguments->pushpull);
            time_total += stats->time_total;
            freeBFSStats(stats);
        }
        break;
        }

        arguments->root = generateRandomRootGeneral(graph, arguments);
        trials--;
    }

    generateGraphPrintMessageWithtime("*     -----> Trials Avg Time (Seconds) <-----", (time_total / (double)arguments->trials));

    if(arguments->verbosity > 0)
    {
        char *fname_txt = (char *) malloc((strlen(arguments->fnameb) + 50) * sizeof(char));
        sprintf(fname_txt, "%s_%d_%d_%d_%d.%s", arguments->fnameb, arguments->algorithm, arguments->datastructure, arguments->trials, arguments->pushpull, "perf");
        FILE *fptr;
        fptr = fopen(fname_txt, "a+");
        fprintf(fptr, "%u %lf \n", arguments->numThreads, (time_total / (double)arguments->trials));
        fclose(fptr);
        free(fname_txt);
    }
}

uint32_t generateRandomRootGraphCSR(struct GraphCSR *graph)
{

    uint32_t root = 0;

    while(1)
    {
        root = generateRandInt(mt19937var);
        if(root < graph->num_vertices)
        {
            if(graph->vertices->out_degree[root] > 0)
                break;
        }
    }

    return root;

}


uint32_t generateRandomRootGraphGrid(struct GraphGrid *graph)
{

    uint32_t root = 0;

    while(1)
    {
        root = generateRandInt(mt19937var);
        if(root < graph->num_vertices)
        {
            if(graph->grid->out_degree[root] > 0)
                break;
        }
    }

    return root;

}

uint32_t generateRandomRootGraphAdjLinkedList(struct GraphAdjLinkedList *graph)
{

    uint32_t root = 0;

    while(1)
    {
        root = generateRandInt(mt19937var);
        if(root < graph->num_vertices)
        {
            if(graph->vertices[root].out_degree > 0)
                break;
        }
    }

    return root;

}

uint32_t generateRandomRootGraphAdjArrayList(struct GraphAdjArrayList *graph)
{

    uint32_t root = 0;

    while(1)
    {
        root = generateRandInt(mt19937var);
        if(root < graph->num_vertices)
        {
            if(graph->vertices[root].out_degree > 0)
                break;
        }
    }

    return root;

}

uint32_t generateRandomRootGeneral(void *graph, struct Arguments *arguments)
{

    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        arguments->root = generateRandomRootGraphCSR(graphCSR);
        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        arguments->root = generateRandomRootGraphGrid(graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        arguments->root = generateRandomRootGraphAdjLinkedList(graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        arguments->root = generateRandomRootGraphAdjArrayList(graphAdjArrayList);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        arguments->root = generateRandomRootGraphCSR(graphCSR);
        break;
    }

    return arguments->root;

}

struct BFSStats *runBreadthFirstSearchAlgorithm(void *graph, uint32_t datastructure, int root, uint32_t pushpull)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BFSStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = breadthFirstSearchGraphCSR(root, pushpull, graphCSR);
        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;

        stats = breadthFirstSearchGraphGrid(root, pushpull, graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        stats = breadthFirstSearchGraphAdjLinkedList(root,  pushpull, graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        stats = breadthFirstSearchGraphAdjArrayList(root, pushpull, graphAdjArrayList);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = breadthFirstSearchGraphCSR(root, pushpull, graphCSR);
        break;
    }

    return stats;

}

struct DFSStats *runDepthFirstSearchAlgorithm(void *graph, uint32_t datastructure, int root)
{


    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct DFSStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = depthFirstSearchGraphCSR(root, graphCSR);
        break;

    case 1: // Grid
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);

        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);

        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);

        break;


    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = depthFirstSearchGraphCSR(root, graphCSR);
        break;
    }

    return stats;

}

struct CCStats *runConnectedComponentsAlgorithm(void *graph, uint32_t datastructure, uint32_t iterations, uint32_t pushpull)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct CCStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = connectedComponentsGraphCSR(iterations, pushpull, graphCSR);
        break;
    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        stats = connectedComponentsGraphGrid(iterations, pushpull, graphGrid);
        break;
    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = connectedComponentsGraphAdjLinkedList(iterations, pushpull, graphAdjLinkedList);
        break;
    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = connectedComponentsGraphAdjArrayList(iterations, pushpull, graphAdjArrayList);
        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = connectedComponentsGraphCSR(iterations, pushpull, graphCSR);
        break;
    }


    return stats;

}

struct TCStats *runTriangleCountAlgorithm(void *graph, uint32_t datastructure, uint32_t pushpull)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct TCStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = triangleCountGraphCSR(pushpull, graphCSR);
        break;
    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        stats = triangleCountGraphGrid(pushpull, graphGrid);
        break;
    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = triangleCountGraphAdjLinkedList(pushpull, graphAdjLinkedList);
        break;
    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = triangleCountGraphAdjArrayList(pushpull, graphAdjArrayList);
        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = triangleCountGraphCSR(pushpull, graphCSR);
        break;
    }


    return stats;

}



struct SPMVStats *runSPMVAlgorithm(void *graph, uint32_t datastructure, uint32_t iterations, uint32_t pushpull)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct SPMVStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = SPMVGraphCSR( iterations, pushpull, graphCSR);

        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        stats = SPMVGraphGrid(iterations, pushpull, graphGrid);

        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = SPMVGraphAdjLinkedList(iterations, pushpull, graphAdjLinkedList);

        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = SPMVGraphAdjArrayList(iterations, pushpull, graphAdjArrayList);

        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = SPMVGraphCSR(iterations, pushpull, graphCSR);

        break;
    }


    return stats;

}


struct IncrementalAggregationStats *runIncrementalAggregationAlgorithm(void *graph, uint32_t datastructure)
{


    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct IncrementalAggregationStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = incrementalAggregationGraphCSR(graphCSR);

        generateGraphPrintMessageWithtime("BUGGY IMPLEMENTATION UNCOMMENT IF YOU NEED IT", 0);
        break;

    case 1: // Grid
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;


    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = incrementalAggregationGraphCSR(graphCSR);

        generateGraphPrintMessageWithtime("BUGGY IMPLEMENTATION UNCOMMENT IF YOU NEED IT", 0);
        break;
    }

    return stats;

}


struct BetweennessCentralityStats *runBetweennessCentralityAlgorithm(void *graph, uint32_t datastructure, uint32_t iterations, uint32_t pushpull)
{
    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BetweennessCentralityStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = betweennessCentralityGraphCSR(iterations, pushpull, graphCSR);
        break;

    case 1: // Grid
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;


    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = betweennessCentralityGraphCSR(iterations, pushpull, graphCSR);
        break;
    }

    // if you want to output pageranks and rankins sorted use this
    stats->realRanks = radixSortEdgesByPageRank (stats->betweennessCentrality, stats->realRanks, stats->num_vertices);

    return stats;

}


struct PageRankStats *runPageRankAlgorithm(void *graph, uint32_t datastructure, double epsilon, uint32_t iterations, uint32_t pushpull)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct PageRankStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = pageRankGraphCSR(epsilon, iterations, pushpull, graphCSR);

        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        stats = pageRankGraphGrid(epsilon, iterations, pushpull, graphGrid);

        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = pageRankGraphAdjLinkedList(epsilon, iterations, pushpull, graphAdjLinkedList);

        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = pageRankGraphAdjArrayList(epsilon, iterations, pushpull, graphAdjArrayList);

        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = pageRankGraphCSR(epsilon, iterations, pushpull, graphCSR);

        break;
    }


    // if you want to output pageranks and rankins sorted use this
    stats->realRanks = radixSortEdgesByPageRank (stats->pageRanks, stats->realRanks, stats->num_vertices);
    return stats;


}

struct BellmanFordStats *runBellmanFordAlgorithm(void *graph, uint32_t datastructure, uint32_t root, uint32_t iterations, uint32_t pushpull)
{

    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BellmanFordStats *stats = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = bellmanFordGraphCSR(root, iterations, pushpull, graphCSR);
        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;

        stats = bellmanFordGraphGrid(root, iterations, pushpull, graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        stats = bellmanFordGraphAdjLinkedList(root, iterations, pushpull, graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        stats = bellmanFordGraphAdjArrayList(root, iterations, pushpull, graphAdjArrayList);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = bellmanFordGraphCSR(root, iterations, pushpull, graphCSR);
        break;
    }

    return stats;

}


struct SSSPStats *runSSSPAlgorithm(void *graph, uint32_t datastructure, uint32_t root, uint32_t iterations, uint32_t pushpull, uint32_t delta)
{

    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;

    struct SSSPStats *stats = NULL;
    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = SSSPGraphCSR(root, iterations, pushpull, graphCSR, delta);

        break;

    case 1: // Grid
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = SSSPGraphCSR(root, iterations, pushpull, graphCSR, delta);
        break;
    }

    return stats;

}



void freeGraphDataStructure(void *graph, uint32_t datastructure)
{
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;

    switch (datastructure)
    {
    case 0: // CSR
        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
        break;

    case 1: // Grid
        graphGrid = (struct GraphGrid *)graph;
        Start(timer);
        graphGridFree(graphGrid);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph Grid (Seconds)", Seconds(timer));
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        Start(timer);
        graphAdjLinkedListFree(graphAdjLinkedList);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph Adjacency Linked List (Seconds)", Seconds(timer));
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        Start(timer);
        graphAdjArrayListFree(graphAdjArrayList);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph Adjacency Array List (Seconds)", Seconds(timer));
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
        break;
    }

    free(timer);

}


void freeGraphStatsGeneral(void *stats, uint32_t algorithm)
{

    switch (algorithm)
    {
    case 0:  // bfs
    {
        struct BFSStats *freeStatsBFS = (struct BFSStats * )stats;
        freeBFSStats(freeStatsBFS);
    }
    break;
    case 1: // pagerank
    {
        struct PageRankStats *freeStatsPageRank = (struct PageRankStats * )stats;
        freePageRankStats(freeStatsPageRank);
    }
    break;
    case 2: // SSSP-Dijkstra
    {
        struct SSSPStats *freeStatsSSSP = (struct SSSPStats * )stats;
        freeSSSPStats(freeStatsSSSP);
    }
    break;
    case 3: // SSSP-Bellmanford
    {
        struct BellmanFordStats *freeStatsBellmanFord = (struct BellmanFordStats * )stats;
        freeBellmanFordStats(freeStatsBellmanFord);
    }
    break;
    case 4: // DFS
    {
        struct DFSStats *freeStatsDFS = (struct DFSStats * )stats;
        freeDFSStats(freeStatsDFS);
    }
    break;
    case 5: //SPMV
    {
        struct SPMVStats *freeStats = (struct SPMVStats *)stats;
        freeSPMVStats(freeStats);
    }
    break;
    case 6: // Connected Components
    {
        struct CCStats *freeStats = (struct CCStats *)stats;
        freeCCStats(freeStats);
    }
    break;
    case 7: // Triangle Counting
    {
        struct TCStats *freeStats = (struct TCStats *)stats;
        freeTCStats(freeStats);
    }
    break;
    case 8: // incremental Aggregation
    {
        struct IncrementalAggregationStats *freeStats = (struct IncrementalAggregationStats *)stats;
        freeIncrementalAggregationStats(freeStats);
    }
    break;
    default:// bfs file
    {
        struct BFSStats *freeStatsBFS = (struct BFSStats *)stats;
        freeBFSStats(freeStatsBFS);
    }
    break;
    }

}
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
#include <omp.h>

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



#include "sortRun.h"
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

        arguments->lmode = arguments->lmode_l2;
        if(arguments->lmode)
            edgeList = reorderGraphProcess(edgeList, arguments);

        arguments->lmode = arguments->lmode_l3;
        if(arguments->lmode)
            edgeList = reorderGraphProcess(edgeList, arguments);

        if(arguments->mmode)
            edgeList = maskGraphProcess(edgeList, arguments);

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
        Start(timer);
        struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted); // read edglist from binary file
        Stop(timer);
        // edgeListPrint(edgeList);
        graphCSRPrintMessageWithtime("Read Edge List From File (Seconds)", Seconds(timer));

        edgeList = sortRunAlgorithms(edgeList, arguments->sort);

        if(arguments->dflag)
        {
            Start(timer);
            edgeList = removeDulpicatesSelfLoopEdges(edgeList);
            Stop(timer);
            graphCSRPrintMessageWithtime("Removing duplicate edges (Seconds)", Seconds(timer));
        }

        if(arguments->lmode)
        {
            edgeList = reorderGraphProcess(edgeList, arguments);
            edgeList = sortRunAlgorithms(edgeList, arguments->sort);
        }

        // add another layer 2 of reordering to test how DBG affect Gorder, or Gorder affect Rabbit order ...etc
        arguments->lmode = arguments->lmode_l2;
        if(arguments->lmode)
        {
            edgeList = reorderGraphProcess(edgeList, arguments);
            edgeList = sortRunAlgorithms(edgeList, arguments->sort);
        }

        arguments->lmode = arguments->lmode_l3;
        if(arguments->lmode)
        {
            edgeList = reorderGraphProcess(edgeList, arguments);
            edgeList = sortRunAlgorithms(edgeList, arguments->sort);
        }

        if(arguments->mmode)
            edgeList = maskGraphProcess(edgeList, arguments);


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

    printf("*-----------------------------------------------------*\n");
    printf("| %-35s %-15d | \n", "Number of Threads Preprocessing:", arguments->pre_numThreads);
    printf(" -----------------------------------------------------\n");

    omp_set_num_threads(arguments->pre_numThreads);

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
            graph = (void *)graphAdjArrayListPreProcessingStep (arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphAdjArrayList Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        case 6: // CSR Dual Order
            Start(timer);
            graph = (void *)graphCSRPreProcessingStepDualOrder(arguments);
            Stop(timer);
            generateGraphPrintMessageWithtime("GraphCSR DO Preprocessing Step Time (Seconds)", Seconds(timer));
            break;
        default:// CSR
            Start(timer);
            graph = (void *)graphCSRPreProcessingStep (arguments);
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

void runGraphAlgorithms(struct Arguments *arguments, void *graph)
{
    printf("*-----------------------------------------------------*\n");
    printf("| %-35s %-15d | \n", "Number of Threads Algorithm :", arguments->algo_numThreads);
    printf(" -----------------------------------------------------\n");
    printf("*-----------------------------------------------------*\n");
    printf("| %-35s %-15d | \n", "Number of Threads Kernel    :", arguments->ker_numThreads);
    printf(" -----------------------------------------------------\n");
    omp_set_num_threads(arguments->algo_numThreads);

    double time_total = 0.0f;
    uint32_t  trials = arguments->trials;

    while(trials)
    {
        switch (arguments->algorithm)
        {
        case 0:  // BFS
        {
            struct BFSStats *stats = runBreadthFirstSearchAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeBFSStats(stats);
            }
        }
        break;
        case 1: // pagerank
        {
            struct PageRankStats *stats = runPageRankAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;

                if(arguments->Sflag) // output page rank error statistics
                {
                    arguments->pushpull = 0;
                    struct PageRankStats *ref_stats = runPageRankAlgorithm(arguments, graph);
                    collectStatsPageRank(arguments, ref_stats, stats, trials);
                    freePageRankStats(ref_stats);
                }

                freePageRankStats(stats);
            }
        }
        break;
        case 2: // SSSP-Delta
        {
            struct SSSPStats *stats = runSSSPAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeSSSPStats(stats);
            }
        }
        break;
        case 3: // SSSP-Bellmanford
        {
            struct BellmanFordStats *stats = runBellmanFordAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeBellmanFordStats(stats);
            }
        }
        break;
        case 4: // DFS
        {
            struct DFSStats *stats = runDepthFirstSearchAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeDFSStats(stats);
            }
        }
        break;
        case 5: // SPMV
        {
            struct SPMVStats *stats = runSPMVAlgorithm(arguments,  graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeSPMVStats(stats);
            }
        }
        break;
        case 6: // Connected Components
        {
            struct CCStats *stats = runConnectedComponentsAlgorithm(arguments,  graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeCCStats(stats);
            }
        }
        break;
        case 7: // Betweenness Centrality
        {
            struct BetweennessCentralityStats *stats = runBetweennessCentralityAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeBetweennessCentralityStats(stats);
            }
        }
        break;
        case 8: // Triangle Counting
        {
            struct TCStats *stats = runTriangleCountAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeTCStats(stats);
            }
        }
        break;
        case 9: // incremental Aggregation
        {
            struct IncrementalAggregationStats *stats = runIncrementalAggregationAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeIncrementalAggregationStats(stats);
            }
        }
        break;

        default: // BFS
        {
            struct BFSStats *stats = runBreadthFirstSearchAlgorithm(arguments, graph);
            if(stats)
            {
                time_total += stats->time_total;
                freeBFSStats(stats);
            }
        }
        break;
        }

        arguments->source = generateRandomRootGeneral(arguments, graph);
        trials--;
    }

    generateGraphPrintMessageWithtime("*     -----> Trials Avg Time (Seconds) <-----", (time_total / (double)arguments->trials));

    if(arguments->verbosity > 0)
    {
        char *fname_txt = (char *) malloc((strlen(arguments->fnameb) + 50) * sizeof(char));
        sprintf(fname_txt, "%s_%d_%d_%d_%d.%s", arguments->fnameb, arguments->algorithm, arguments->datastructure, arguments->trials, arguments->pushpull, "perf");
        FILE *fptr;
        fptr = fopen(fname_txt, "a+");
        fprintf(fptr, "%u %lf \n", arguments->algo_numThreads, (time_total / (double)arguments->trials));
        fclose(fptr);
        free(fname_txt);
    }
}

uint32_t generateRandomRootGraphCSR(mt19937state *mt19937var, struct GraphCSR *graph)
{

    uint32_t source = 0;
    uint32_t source_temp = 0;

    while(1)
    {
        source = generateRandInt(mt19937var);
        if(source < graph->num_vertices)
        {
            source_temp = graph->sorted_edges_array->label_array[source];
            if(graph->vertices->out_degree[source_temp] > graph->avg_degree)
                break;
        }
    }

    return source;

}


uint32_t generateRandomRootGraphGrid(mt19937state *mt19937var, struct GraphGrid *graph)
{

    uint32_t source = 0;

    while(1)
    {
        source = generateRandInt(mt19937var);
        if(source < graph->num_vertices)
        {
            if(graph->grid->out_degree[source] > graph->avg_degree)
                break;
        }
    }

    return source;

}

uint32_t generateRandomRootGraphAdjLinkedList(mt19937state *mt19937var, struct GraphAdjLinkedList *graph)
{

    uint32_t source = 0;

    while(1)
    {
        source = generateRandInt(mt19937var);
        if(source < graph->num_vertices)
        {
            if(graph->vertices[source].out_degree > graph->avg_degree)
                break;
        }
    }

    return source;

}

uint32_t generateRandomRootGraphAdjArrayList(mt19937state *mt19937var, struct GraphAdjArrayList *graph)
{

    uint32_t source = 0;

    while(1)
    {
        source = generateRandInt(mt19937var);
        if(source < graph->num_vertices)
        {
            if(graph->vertices[source].out_degree > graph->avg_degree)
                break;
        }
    }

    return source;

}

uint32_t generateRandomRootGeneral(struct Arguments *arguments, void *graph)
{

    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
    case 6:
        graphCSR = (struct GraphCSR *)graph;
        arguments->source = generateRandomRootGraphCSR(&(arguments->mt19937var), graphCSR);
        break;

    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;
        arguments->source = generateRandomRootGraphGrid(&(arguments->mt19937var), graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        arguments->source = generateRandomRootGraphAdjLinkedList(&(arguments->mt19937var), graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        arguments->source = generateRandomRootGraphAdjArrayList(&(arguments->mt19937var), graphAdjArrayList);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        arguments->source = generateRandomRootGraphCSR(&(arguments->mt19937var), graphCSR);
        break;
    }

    return arguments->source;

}

struct BFSStats *runBreadthFirstSearchAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BFSStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;

        stats = breadthFirstSearchGraphCSR(arguments, graphCSR);
        break;

    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;

        stats = breadthFirstSearchGraphGrid(arguments, graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        stats = breadthFirstSearchGraphAdjLinkedList(arguments, graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        stats = breadthFirstSearchGraphAdjArrayList(arguments, graphAdjArrayList);
        break;

    case 6: // CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = breadthFirstSearchGraphCSRDualOrder(arguments, graphCSR);
        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = breadthFirstSearchGraphCSR(arguments, graphCSR);
        break;
    }

    return stats;

}

struct DFSStats *runDepthFirstSearchAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct DFSStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = depthFirstSearchGraphCSR(arguments, graphCSR);
        break;

    case 1: // Grid
    case 5:
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

    case 6: // CSR
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);

        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = depthFirstSearchGraphCSR(arguments, graphCSR);
        break;
    }

    return stats;

}

struct CCStats *runConnectedComponentsAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct CCStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4: // CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = connectedComponentsGraphCSR(arguments, graphCSR);
        break;
    case 1: // Grid
    case 5: // Grid
        graphGrid = (struct GraphGrid *)graph;
        stats = connectedComponentsGraphGrid(arguments, graphGrid);
        break;
    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = connectedComponentsGraphAdjLinkedList(arguments, graphAdjLinkedList);
        break;
    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = connectedComponentsGraphAdjArrayList(arguments, graphAdjArrayList);
        break;
    case 6: // CSR
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = connectedComponentsGraphCSR(arguments, graphCSR);
        break;
    }


    return stats;

}

struct TCStats *runTriangleCountAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct TCStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = triangleCountGraphCSR(arguments, graphCSR);
        break;
    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;
        stats = triangleCountGraphGrid(arguments, graphGrid);
        break;
    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = triangleCountGraphAdjLinkedList(arguments, graphAdjLinkedList);
        break;
    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = triangleCountGraphAdjArrayList(arguments, graphAdjArrayList);
        break;
    case 6: // CSR
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = triangleCountGraphCSR(arguments, graphCSR);
        break;
    }


    return stats;

}



struct SPMVStats *runSPMVAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct SPMVStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = SPMVGraphCSR(arguments, graphCSR);

        break;

    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;
        stats = SPMVGraphGrid(arguments, graphGrid);

        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = SPMVGraphAdjLinkedList(arguments, graphAdjLinkedList);

        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = SPMVGraphAdjArrayList(arguments, graphAdjArrayList);

        break;

    case 6: // CSR

        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);

        break;
    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = SPMVGraphCSR(arguments, graphCSR);

        break;
    }


    return stats;

}


struct IncrementalAggregationStats *runIncrementalAggregationAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct IncrementalAggregationStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = incrementalAggregationGraphCSR(graphCSR);

        generateGraphPrintMessageWithtime("BUGGY IMPLEMENTATION UNCOMMENT IF YOU NEED IT", 0);
        break;

    case 1: // Grid
    case 5:
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

    case 6: // CSR

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


struct BetweennessCentralityStats *runBetweennessCentralityAlgorithm(struct Arguments *arguments, void *graph)
{
    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BetweennessCentralityStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = betweennessCentralityGraphCSR(arguments, graphCSR);
        break;

    case 1: // Grid
    case 5:
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("Betweenness Centrality NOT YET IMPLEMENTED", 0);
        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        generateGraphPrintMessageWithtime("Betweenness Centrality NOT YET IMPLEMENTED", 0);
        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        generateGraphPrintMessageWithtime("Betweenness Centrality NOT YET IMPLEMENTED", 0);
        break;

    case 6: // CSR
        generateGraphPrintMessageWithtime("Betweenness Centrality NOT YET IMPLEMENTED", 0);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = betweennessCentralityGraphCSR(arguments, graphCSR);
        break;
    }

    // if you want to output pageranks and rankins sorted use this
    if(stats)
    {
        stats->realRanks = radixSortEdgesByPageRank (stats->betweennessCentrality, stats->realRanks, stats->num_vertices);
        printRanksBetweennessCentralityStats(stats);
    }
    return stats;

}


struct PageRankStats *runPageRankAlgorithm(struct Arguments *arguments, void *graph)
{


    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct PageRankStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        stats = pageRankGraphCSR(arguments, graphCSR);

        break;

    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;
        stats = pageRankGraphGrid(arguments, graphGrid);

        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;
        stats = pageRankGraphAdjLinkedList(arguments, graphAdjLinkedList);

        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;
        stats = pageRankGraphAdjArrayList(arguments, graphAdjArrayList);

        break;

    case 6: // CSR
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;
        stats = pageRankGraphCSR(arguments, graphCSR);

        break;
    }


    // if you want to output pageranks and rankins sorted use this
    if(stats)
    {
        stats->realRanks = radixSortEdgesByPageRank (stats->pageRanks, stats->realRanks, stats->num_vertices);
    }

    return stats;


}

struct BellmanFordStats *runBellmanFordAlgorithm(struct Arguments *arguments, void *graph)
{

    struct GraphCSR *graphCSR = NULL;
    struct GraphGrid *graphGrid = NULL;
    struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    struct GraphAdjArrayList *graphAdjArrayList = NULL;
    struct BellmanFordStats *stats = NULL;

    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;

        stats = bellmanFordGraphCSR(arguments, graphCSR);
        break;

    case 1: // Grid
    case 5:
        graphGrid = (struct GraphGrid *)graph;

        stats = bellmanFordGraphGrid(arguments, graphGrid);
        break;

    case 2: // Adj Linked List
        graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        stats = bellmanFordGraphAdjLinkedList(arguments, graphAdjLinkedList);
        break;

    case 3: // Adj Array List
        graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        stats = bellmanFordGraphAdjArrayList(arguments, graphAdjArrayList);
        break;

    case 6: // CSR
        generateGraphPrintMessageWithtime("NOT YET IMPLEMENTED", 0);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = bellmanFordGraphCSR(arguments, graphCSR);
        break;
    }

    return stats;

}


struct SSSPStats *runSSSPAlgorithm(struct Arguments *arguments, void *graph)
{

    struct GraphCSR *graphCSR = NULL;
    // struct GraphGrid *graphGrid = NULL;
    // struct GraphAdjLinkedList *graphAdjLinkedList = NULL;
    // struct GraphAdjArrayList *graphAdjArrayList = NULL;

    struct SSSPStats *stats = NULL;
    switch (arguments->datastructure)
    {
    case 0: // CSR
    case 4:
        graphCSR = (struct GraphCSR *)graph;

        stats = SSSPGraphCSR(arguments, graphCSR);

        break;

    case 1: // Grid
    case 5:
        // graphGrid = (struct GraphGrid *)graph;
        generateGraphPrintMessageWithtime("SSSP NOT YET IMPLEMENTED", 0);
        break;

    case 2: // Adj Linked List
        // graphAdjLinkedList = (struct GraphAdjLinkedList *)graph;

        generateGraphPrintMessageWithtime("SSSP NOT YET IMPLEMENTED", 0);
        break;

    case 3: // Adj Array List
        // graphAdjArrayList = (struct GraphAdjArrayList *)graph;

        generateGraphPrintMessageWithtime("SSSP NOT YET IMPLEMENTED", 0);
        break;

    case 6: // CSR
        generateGraphPrintMessageWithtime("SSSP NOT YET IMPLEMENTED", 0);
        break;

    default:// CSR
        graphCSR = (struct GraphCSR *)graph;

        stats = SSSPGraphCSR(arguments, graphCSR);
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
    case 4:
        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
        break;

    case 1: // Grid
    case 5:
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

    case 6: // CSR
        graphCSR = (struct GraphCSR *)graph;
        Start(timer);
        graphCSRFree(graphCSR);
        Stop(timer);
        generateGraphPrintMessageWithtime("Free Graph CSR (Seconds)", Seconds(timer));
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
    case 2: // SSSP-Delta
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
    case 7: // Betweenness Centrality
    {
        struct BetweennessCentralityStats *freeStats = (struct BetweennessCentralityStats *)stats;
        freeBetweennessCentralityStats(freeStats);
    }
    break;
    case 8: // Triangle Counting
    {
        struct TCStats *freeStats = (struct TCStats *)stats;
        freeTCStats(freeStats);
    }
    break;
    case 9: // incremental Aggregation
    {
        struct IncrementalAggregationStats *freeStats = (struct IncrementalAggregationStats *)stats;
        freeIncrementalAggregationStats(freeStats);
    }
    break;
    default:// BFS
    {
        struct BFSStats *freeStatsBFS = (struct BFSStats *)stats;
        freeBFSStats(freeStatsBFS);
    }
    break;
    }

}

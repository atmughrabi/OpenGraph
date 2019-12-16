// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphGrid.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "timer.h"
#include "graphConfig.h"
#include "myMalloc.h"

#include "grid.h"
#include "graphGrid.h"
#include "edgeList.h"
#include "sortRun.h"
#include "reorder.h"







// #include "countsort.h"
// #include "radixsort.h"





void  graphGridReset(struct GraphGrid *graphGrid)
{

    graphGridResetActivePartitionsMap(graphGrid->grid);

}

void  graphGridPrint(struct GraphGrid *graphGrid)
{


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Graph Grid Properties");
    printf(" -----------------------------------------------------\n");
#if WEIGHTED
    printf("| %-51s | \n", "WEIGHTED");
#else
    printf("| %-51s | \n", "UN-WEIGHTED");
#endif

#if DIRECTED
    printf("| %-51s | \n", "DIRECTED");
#else
    printf("| %-51s | \n", "UN-DIRECTED");
#endif
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Vertices (V)");
    printf("| %-51u | \n", graphGrid->grid->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Edges (E)");
    printf("| %-51u | \n", graphGrid->grid->num_edges);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Partitions (P)");
    printf("| %-51u | \n", graphGrid->grid->num_partitions);
    printf(" -----------------------------------------------------\n");




    //   uint32_t i;
    //    for ( i = 0; i < ( graphGrid->grid->num_partitions*graphGrid->grid->num_partitions); ++i)
    //       {

    //       uint32_t x = i % graphGrid->grid->num_partitions;    // % is the "modulo operator", the remainder of i / width;
    // uint32_t y = i / graphGrid->grid->num_partitions;

    //      if(graphGrid->grid->partitions[i].num_edges){

    //       printf("| %-11s (%u,%u) \n", "Partition: ", y, x);
    //      printf("| %-11s %-40u   \n", "Edges: ", graphGrid->grid->partitions[i].num_edges);
    //      printf("| %-11s %-40u   \n", "Vertices: ", graphGrid->grid->partitions[i].num_vertices);
    //      edgeListPrint(graphGrid->grid->partitions[i].edgeList);
    //       }

    //       }


}


struct GraphGrid *graphGridNew(struct EdgeList *edgeList)
{


    struct GraphGrid *graphGrid = (struct GraphGrid *) my_malloc( sizeof(struct GraphGrid));

#if WEIGHTED
    graphGrid->max_weight =  edgeList->max_weight;
#endif

    graphGrid->num_edges = edgeList->num_edges;
    graphGrid->num_vertices = edgeList->num_vertices;

    graphGrid->grid = gridNew(edgeList);


    return graphGrid;

}

void   graphGridFree(struct GraphGrid *graphGrid)
{

    if(graphGrid->grid)
        gridFree(graphGrid->grid);

    if(graphGrid)
        free(graphGrid);

}



struct GraphGrid *graphGridPreProcessingStep (struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    printf("Filename : %s \n", arguments->fnameb);


    Start(timer);
    struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted);
    Stop(timer);
    // edgeListPrint(edgeList);
    graphGridPrintMessageWithtime("Read Edge List From File (Seconds)", Seconds(timer));



    if(arguments->lmode)
        edgeList = reorderGraphProcess(edgeList, arguments);

    // Start(timer);
    edgeList = sortRunAlgorithms(edgeList, arguments->sort);

    if(arguments->dflag)
    {
        Start(timer);
        edgeList = removeDulpicatesSelfLoopEdges(edgeList);
        Stop(timer);
        graphCSRPrintMessageWithtime("Removing duplicate edges (Seconds)", Seconds(timer));
    }
    // Stop(timer);
    // graphGridPrintMessageWithtime("Radix Sort Edges By Source (Seconds)",Seconds(timer));

    Start(timer);
    struct GraphGrid *graphGrid = graphGridNew(edgeList);
    Stop(timer);
    graphGridPrintMessageWithtime("Create Graph Grid (Seconds)", Seconds(timer));


    graphGridPrint(graphGrid);


    freeEdgeList(edgeList);
    free(timer);
    return graphGrid;


}



void graphGridPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}


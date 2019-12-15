// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : test_grid.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:29
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "libcxl.h"
#include "capienv.h"

#include "adjLinkedList.h"
#include "dynamicQueue.h"
#include "edgeList.h"

#include "countsort.h"
#include "radixsort.h"

#include "graphCSR.h"
#include "graphAdjLinkedList.h"
#include "graphAdjArrayList.h"
#include "graphGrid.h"

#include "vertex.h"
#include "timer.h"
#include "BFS.h"





void printMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}


int main()
{
    // create the graph given in above fugure
    // int V = 5;

    // const char * fname = "00_Graph_OpenMP/datasets/test/test.txt";
    // const char * fname = "00_Graph_OpenMP/datasets/wiki-vote/wiki-Vote.txt";
    // const char * fname = "00_Graph_OpenMP/datasets/twitter/twitter_rv.txt";
    // const char * fname = "00_Graph_OpenMP/datasets/facebook/facebook_combined.txt";

    // /
    const char *fnameb = "00_Graph_OpenMP/datasets/test/test.txt.bin";
    // const char * fnameb = "00_Graph_OpenMP/datasets/twitter/twitter_rv.txt.bin";
    // const char * fnameb = "00_Graph_OpenMP/datasets/twitter/twitter_rv.txt.bin8";
    // const char * fnameb = "00_Graph_OpenMP/datasets/facebook/facebook_combined.txt.bin";
    // const char * fnameb = "00_Graph_OpenMP/datasets/wiki-vote/wiki-Vote.txt.bin";


    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    printf("Filename : %s \n", fnameb);

    // Start(timer);
    // readEdgeListstxt(fname);
    // Stop(timer);
    // printf("Read Edge List From File converted to binary : %f Seconds \n",Seconds(timer));

    Start(timer);
    struct EdgeList *edgeList = readEdgeListsbin(fnameb, 0);
    Stop(timer);
    // edgeListPrint(edgeList);
    printMessageWithtime("Read Edge List From File (Seconds)", Seconds(timer));

    Start(timer);
    edgeList = radixSortEdgesBySourceAndDestination(edgeList);
    Stop(timer);
    printMessageWithtime("Radix Sort Edges By Source (Seconds)", Seconds(timer));

    Start(timer);
    struct GraphGrid *graphGrid = graphGridNew(edgeList);
    Stop(timer);
    printMessageWithtime("Create Graph Grid (Seconds)", Seconds(timer));

    freeEdgeList(edgeList);


    graphGridPrint(graphGrid);

    Start(timer);
    graphGridFree(graphGrid);
    Stop(timer);
    printMessageWithtime("Free Graph Grid (Seconds)", Seconds(timer));




    free(timer);
    return 0;
}
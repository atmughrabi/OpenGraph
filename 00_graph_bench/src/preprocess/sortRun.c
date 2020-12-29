// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : sortRun.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:35:52
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "countsort.h"
#include "radixsort.h"
#include "edgeList.h"
#include "timer.h"
#include "sortRun.h"

struct EdgeList *sortRunAlgorithms(struct EdgeList *edgeList, uint32_t sort)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    switch (sort)
    {
    case 0: // radix sort by source
        Start(timer);
        edgeList = radixSortEdgesBySource(edgeList);
        Stop(timer);
        sortRunPrintMessageWithtime("Radix Sort Edges By Source (Seconds)", Seconds(timer));
        break;
    case 1: // radix sort by source and destination
        Start(timer);
        edgeList = radixSortEdgesBySourceAndDestination(edgeList);
        Stop(timer);
        sortRunPrintMessageWithtime("Radix Sort Edges By Source/Dest (Seconds)", Seconds(timer));
        break;
    case 2: // count sort by source
        Start(timer);
        edgeList = countSortEdgesBySource(edgeList);
        Stop(timer);
        sortRunPrintMessageWithtime("Count Sort Edges By Source (Seconds)", Seconds(timer));
        break;
    case 3: // count sort by source and destination
        Start(timer);
        edgeList = countSortEdgesBySourceAndDestination(edgeList);
        Stop(timer);
        sortRunPrintMessageWithtime("Count Sort Edges By Source/Dest (Seconds)", Seconds(timer));
        break;
    default:// bfs file name root
        Start(timer);
        edgeList = radixSortEdgesBySource(edgeList);
        Stop(timer);
        sortRunPrintMessageWithtime("Radix Sort Edges By Source (Seconds)", Seconds(timer));
        break;
    }


    free(timer);
    return edgeList;

}

void sortRunPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}
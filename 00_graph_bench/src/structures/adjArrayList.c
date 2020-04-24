// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : adjArrayList.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "edgeList.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "adjArrayList.h"


void adjArrayListPrint(struct AdjArrayList *adjArrayList)
{

    if(adjArrayList->out_degree)
    {
        edgeListPrint(adjArrayList->outNodes);
    }

#if DIRECTED
    if(adjArrayList->in_degree)
    {
        edgeListPrint(adjArrayList->inNodes);
    }
#endif

}

struct AdjArrayList *adjArrayListNew()
{

    struct AdjArrayList *newNode = (struct AdjArrayList *) my_malloc(sizeof(struct AdjArrayList));

    newNode->out_degree = 0;
    newNode->outNodes = NULL;

#if DIRECTED
    newNode->in_degree = 0;
    newNode->inNodes = NULL;
#endif

    return newNode;

}

struct AdjArrayList *adjArrayListCreateNeighbourList(struct AdjArrayList *adjArrayList)
{

    if(adjArrayList->out_degree)
        adjArrayList->outNodes = newEdgeList(adjArrayList->out_degree);

#if DIRECTED
    if(adjArrayList->in_degree)
        adjArrayList->inNodes = newEdgeList(adjArrayList->in_degree);
#endif

    return adjArrayList;
}

struct AdjArrayList *adjArrayListCreateNeighbourListOutNodes(struct AdjArrayList *adjArrayList)
{

    adjArrayList->outNodes = newEdgeList(adjArrayList->out_degree);

    return adjArrayList;
}



struct AdjArrayList *adjArrayListCreateNeighbourListInNodes(struct AdjArrayList *adjArrayList)
{


#if DIRECTED
    adjArrayList->inNodes = newEdgeList(adjArrayList->in_degree);
#endif

    return adjArrayList;
}


void adjArrayListFree(struct AdjArrayList *adjArrayList)
{

    if(adjArrayList)
    {
        if(adjArrayList->outNodes)
            freeEdgeList(adjArrayList->outNodes);
#if DIRECTED
        if(adjArrayList->inNodes)
            freeEdgeList(adjArrayList->inNodes);
#endif
        free(adjArrayList);
    }
}

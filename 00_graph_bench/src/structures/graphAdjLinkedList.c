// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphAdjLinkedList.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <linux/types.h>
#include <omp.h>

#include "timer.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "edgeList.h"
#include "vertex.h"
#include "sortRun.h"
#include "reorder.h"

#include "adjLinkedList.h"
#include "graphAdjLinkedList.h"





// A utility function that creates a graphAdjLinkedList of V vertices
struct GraphAdjLinkedList *graphAdjLinkedListGraphNew(__u32 V)
{

    struct GraphAdjLinkedList *graphAdjLinkedList = (struct GraphAdjLinkedList *) my_malloc( sizeof(struct GraphAdjLinkedList));
    graphAdjLinkedList->num_vertices = V;
    graphAdjLinkedList->vertices = (struct AdjLinkedList *) my_malloc( V * sizeof(struct AdjLinkedList));


    __u32 i;
    #pragma omp parallel for
    for(i = 0; i < V; i++)
    {

        graphAdjLinkedList->vertices[i].outNodes = NULL;
        graphAdjLinkedList->vertices[i].out_degree = 0;

#if DIRECTED
        graphAdjLinkedList->vertices[i].inNodes = NULL;
        graphAdjLinkedList->vertices[i].in_degree = 0;
#endif

        graphAdjLinkedList->vertices[i].visited = 0;
    }
    // printf("\n Success!!! V: %d\n ", V);

    return graphAdjLinkedList;

}

struct GraphAdjLinkedList *graphAdjLinkedListEdgeListNew(struct EdgeList *edgeList)
{

    struct GraphAdjLinkedList *graphAdjLinkedList = (struct GraphAdjLinkedList *) my_malloc( sizeof(struct GraphAdjLinkedList));


    graphAdjLinkedList->num_vertices = edgeList->num_vertices;
    graphAdjLinkedList->num_edges = edgeList->num_edges;
    graphAdjLinkedList->vertices = (struct AdjLinkedList *) my_malloc( graphAdjLinkedList->num_vertices * sizeof(struct AdjLinkedList));

#if WEIGHTED
    graphAdjLinkedList->max_weight =  edgeList->max_weight;
#endif



    __u32 i;
    #pragma omp parallel for
    for(i = 0; i < graphAdjLinkedList->num_vertices; i++)
    {

        graphAdjLinkedList->vertices[i].outNodes = NULL;
        graphAdjLinkedList->vertices[i].out_degree = 0;

#if DIRECTED
        graphAdjLinkedList->vertices[i].inNodes = NULL;
        graphAdjLinkedList->vertices[i].in_degree = 0;
#endif

        graphAdjLinkedList->vertices[i].visited = 0;
    }

    omp_lock_t *vertex_lock  = (omp_lock_t *) my_malloc( graphAdjLinkedList->num_vertices * sizeof(omp_lock_t));



    #pragma omp parallel for
    for (i = 0; i < graphAdjLinkedList->num_vertices; i++)
    {
        omp_init_lock(&(vertex_lock[i]));
    }

    // #pragma omp parallel for
    for(i = 0; i < edgeList->num_edges; i++)
    {

        adjLinkedListAddEdge(graphAdjLinkedList, edgeList, i, vertex_lock);

    }

    #pragma omp parallel for
    for (i = 0; i < graphAdjLinkedList->num_vertices; i++)
    {
        omp_destroy_lock(&(vertex_lock[i]));
    }

    free(vertex_lock);

    return graphAdjLinkedList;

}


// A utility function to print the adjacency list
// representation of graphAdjLinkedList
void graphAdjLinkedListPrint(struct GraphAdjLinkedList *graphAdjLinkedList)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "GraphAdjLinkedList Properties");
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
    printf("| %-51u | \n", graphAdjLinkedList->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Edges (E)");
    printf("| %-51u | \n", graphAdjLinkedList->num_edges);
    printf(" -----------------------------------------------------\n");

}

void graphAdjLinkedListFree(struct GraphAdjLinkedList *graphAdjLinkedList)
{

    __u32 v;
    struct AdjLinkedListNode *pCrawl;
    struct AdjLinkedListNode *pFree;

    for (v = 0; v < graphAdjLinkedList->num_vertices; ++v)
    {
        pCrawl = graphAdjLinkedList->vertices[v].outNodes;
        pFree  = graphAdjLinkedList->vertices[v].outNodes;

        while (pCrawl)
        {

            pFree = pCrawl;
            pCrawl = pCrawl->next;
            if(pFree)
                free(pFree);

        }

#if DIRECTED
        pCrawl = graphAdjLinkedList->vertices[v].inNodes;
        pFree  = graphAdjLinkedList->vertices[v].inNodes;

        while (pCrawl)
        {
            pFree = pCrawl;
            pCrawl = pCrawl->next;
            if(pFree)
                free(pFree);
        }
#endif

    }

    if(graphAdjLinkedList->vertices)
        free(graphAdjLinkedList->vertices);

    if(graphAdjLinkedList)
        free(graphAdjLinkedList);


}

void adjLinkedListAddEdge(struct GraphAdjLinkedList *graphAdjLinkedList, struct EdgeList *edge, __u32 i, omp_lock_t *vertex_lock)
{

    // omp_set_lock(&(vertex_lock[edge->edges_array_src[i]]));
    // omp_unset_lock((&vertex_lock[edge->edges_array_src[i]]));

    // Add an edge from src to dest.  A new node is
    // added to the adjacency list of src.  The node
    // is added at the begining
    struct AdjLinkedListNode *newNode = newAdjLinkedListOutNode(edge->edges_array_dest[i]);
#if WEIGHTED
    newNode->weight = edge->edges_array_weight[i];
#endif
    // omp_set_lock(&(vertex_lock[edge->edges_array_src[i]]));
    newNode->next = graphAdjLinkedList->vertices[edge->edges_array_src[i]].outNodes;
    graphAdjLinkedList->vertices[edge->edges_array_src[i]].out_degree++;
    graphAdjLinkedList->vertices[edge->edges_array_src[i]].visited = 0;
    graphAdjLinkedList->vertices[edge->edges_array_src[i]].outNodes = newNode;

    // omp_unset_lock((&vertex_lock[edge->edges_array_src[i]]));

    // omp_set_lock(&(vertex_lock[edge->edges_array_dest[i]]));
    // omp_unset_lock((&vertex_lock[edge->edges_array_dest[i]]));
    // Since graphAdjLinkedList is undirected, add an edge from
    // dest to src also
    newNode = newAdjLinkedListInNode(edge->edges_array_src[i]);
#if WEIGHTED
    newNode->weight = edge->edges_array_weight[i];
#endif
    // omp_set_lock(&(vertex_lock[edge->edges_array_dest[i]]));
#if DIRECTED
    newNode->next = graphAdjLinkedList->vertices[edge->edges_array_dest[i]].inNodes;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].in_degree++;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].visited = 0;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].inNodes = newNode;
#else
    newNode->next = graphAdjLinkedList->vertices[edge->edges_array_dest[i]].outNodes;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].out_degree++;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].visited = 0;
    graphAdjLinkedList->vertices[edge->edges_array_dest[i]].outNodes = newNode;
#endif

    // omp_unset_lock((&vertex_lock[edge->edges_array_dest[i]]));

}



void   graphAdjLinkedListPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}
struct GraphAdjLinkedList *graphAdjLinkedListPreProcessingStep (struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    Start(timer);
    struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted);
    Stop(timer);
    // edgeListPrint(edgeList);

    if(arguments->lmode)
        edgeList = reorderGraphProcess(edgeList, arguments);

    edgeList = sortRunAlgorithms(edgeList, arguments->sort);

    if(arguments->dflag)
    {
        Start(timer);
        edgeList = removeDulpicatesSelfLoopEdges(edgeList);
        Stop(timer);
        graphCSRPrintMessageWithtime("Removing duplicate edges (Seconds)", Seconds(timer));
    }

    graphAdjLinkedListPrintMessageWithtime("Read Edge List From File (Seconds)", Seconds(timer));

    Start(timer);
    struct GraphAdjLinkedList *graphAdjLinkedList = graphAdjLinkedListEdgeListNew(edgeList);
    Stop(timer);
    graphAdjLinkedListPrintMessageWithtime("Create Adj Linked List from EdgeList (Seconds)", Seconds(timer));

    graphAdjLinkedListPrint(graphAdjLinkedList);

    freeEdgeList(edgeList);
    free(timer);
    return graphAdjLinkedList;


}
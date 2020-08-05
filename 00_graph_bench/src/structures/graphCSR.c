// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphCSR.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <err.h>
#include <string.h>

#include "myMalloc.h"
#include "timer.h"
#include "graphConfig.h"

#include "edgeList.h"
#include "sortRun.h"
#include "vertex.h"

#include "graphCSR.h"
#include "reorder.h"

//edgelist prerpcessing
// #include "countsort.h"
// #include "radixsort.h"






void graphCSRFree (struct GraphCSR *graphCSR)
{

    if(graphCSR)
    {
        if(graphCSR->vertices)
            freeVertexArray(graphCSR->vertices);
        if(graphCSR->sorted_edges_array)
            freeEdgeList(graphCSR->sorted_edges_array);

#if DIRECTED
        if(graphCSR->inverse_vertices)
            freeVertexArray(graphCSR->inverse_vertices);
        if(graphCSR->inverse_sorted_edges_array)
            freeEdgeList(graphCSR->inverse_sorted_edges_array);
#endif


        free(graphCSR);

    }

}

void graphCSRPrint(struct GraphCSR *graphCSR)
{


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "GraphCSR Properties");
    printf(" -----------------------------------------------------\n");
#if WEIGHTED
    printf("| %-51s | \n", "WEIGHTED");
    printf("| %-51s | \n", "MAX WEIGHT");
    printf("| %-51f | \n", graphCSR->max_weight);
#else
    printf("| %-51s | \n", "UN-WEIGHTED");
#endif

#if DIRECTED
    printf("| %-51s | \n", "DIRECTED");
#else
    printf("| %-51s | \n", "UN-DIRECTED");
#endif
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Average Degree (D)");
    printf("| %-51u | \n", graphCSR->avg_degree);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Vertices (V)");
    printf("| %-51u | \n", graphCSR->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Edges (E)");
    printf("| %-51u | \n", graphCSR->num_edges);
    printf(" -----------------------------------------------------\n");



    vertexArrayMaxOutdegree(graphCSR->vertices, graphCSR->num_vertices);
#if DIRECTED
    vertexArrayMaxInDegree(graphCSR->inverse_vertices, graphCSR->num_vertices);
#endif

}


struct GraphCSR *graphCSRNew(uint32_t V, uint32_t E, uint8_t inverse)
{

    struct GraphCSR *graphCSR = (struct GraphCSR *) my_malloc( sizeof(struct GraphCSR));
    graphCSR->num_vertices = V;
    graphCSR->num_edges = E;
    graphCSR->avg_degree = E / V;

#if WEIGHTED
    graphCSR->max_weight = 0;
#endif

    graphCSR->vertices = newVertexArray(V);

#if DIRECTED
    if (inverse)
    {
        graphCSR->inverse_vertices = newVertexArray(V);
    }
#endif

    return graphCSR;
}


struct GraphCSR *graphCSRAssignEdgeList (struct GraphCSR *graphCSR, struct EdgeList *edgeList, uint8_t inverse)
{


#if DIRECTED

    if(inverse)
        graphCSR->inverse_sorted_edges_array = edgeList;
    else
        graphCSR->sorted_edges_array = edgeList;

#else

    graphCSR->sorted_edges_array = edgeList;

#endif

#if WEIGHTED
    graphCSR->max_weight =  edgeList->max_weight;
#endif

    return mapVerticesWithInOutDegree (graphCSR, inverse);


}


struct GraphCSR *graphCSRPreProcessingStep (struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
   
    Start(timer);
    struct EdgeList *edgeList = readEdgeListsbin(arguments->fnameb, 0, arguments->symmetric, arguments->weighted); // read edglist from binary file
    Stop(timer);
    // edgeListPrint(edgeList);
    graphCSRPrintMessageWithtime("Read Edge List From File (Seconds)", Seconds(timer));


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

#if DIRECTED
    struct GraphCSR *graphCSR = graphCSRNew(edgeList->num_vertices, edgeList->num_edges, 1);
#else
    struct GraphCSR *graphCSR = graphCSRNew(edgeList->num_vertices, edgeList->num_edges, 0);
#endif



    // edgeListPrint(edgeList);
    Start(timer);
    graphCSR = graphCSRAssignEdgeList (graphCSR, edgeList, 0);
    Stop(timer);

    graphCSRPrintMessageWithtime("Mappign Vertices to CSR (Seconds)", Seconds(timer));

#if DIRECTED

    Start(timer);
    struct EdgeList *inverse_edgeList = readEdgeListsMem(edgeList, 1, 0, 0); // read edglist from memory since we pre loaded it
    Stop(timer);

    graphCSRPrintMessageWithtime("Read Inverse Edge List From Memory (Seconds)", Seconds(timer));

    inverse_edgeList = sortRunAlgorithms(inverse_edgeList, arguments->sort);

    Start(timer);
    graphCSR = graphCSRAssignEdgeList (graphCSR, inverse_edgeList, 1);
    Stop(timer);
    graphCSRPrintMessageWithtime("Process In/Out degrees of Inverse Nodes (Seconds)", Seconds(timer));

#endif


    graphCSRPrint(graphCSR);


    free(timer);

    return graphCSR;


}

void writeToBinFileGraphCSR (const char *fname, struct GraphCSR *graphCSR)
{


    FILE  *pBinary;
    uint32_t vertex_id;

    char *fname_txt = (char *) malloc((strlen(fname) + 10) * sizeof(char));
    char *fname_bin = (char *) malloc((strlen(fname) + 10) * sizeof(char));

    fname_txt = strcpy (fname_txt, fname);
    fname_bin = strcat (fname_txt, ".csr");

    pBinary = fopen(fname_bin, "wb");

    if (pBinary == NULL)
    {
        err(1, "open: %s", fname_bin);
        return ;
    }

    fwrite(&(graphCSR->num_edges), sizeof (graphCSR->num_edges), 1, pBinary);
    fwrite(&(graphCSR->num_vertices), sizeof (graphCSR->num_vertices), 1, pBinary);
#if WEIGHTED
    fwrite(&(graphCSR->max_weight), sizeof (graphCSR->max_weight), 1, pBinary);
#endif

    for(vertex_id = 0; vertex_id < graphCSR->num_vertices ; vertex_id++)
    {

        fwrite(&(graphCSR->vertices->out_degree[vertex_id]), sizeof (graphCSR->vertices->out_degree[vertex_id]), 1, pBinary);
        fwrite(&(graphCSR->vertices->in_degree[vertex_id]), sizeof (graphCSR->vertices->in_degree[vertex_id]), 1, pBinary);
        fwrite(&(graphCSR->vertices->edges_idx[vertex_id]), sizeof (graphCSR->vertices->edges_idx[vertex_id]), 1, pBinary);

#if DIRECTED
        if(graphCSR->inverse_vertices)
        {
            fwrite(&(graphCSR->inverse_vertices->out_degree[vertex_id]), sizeof (graphCSR->inverse_vertices->out_degree[vertex_id]), 1, pBinary);
            fwrite(&(graphCSR->inverse_vertices->in_degree[vertex_id]), sizeof (graphCSR->inverse_vertices->in_degree[vertex_id]), 1, pBinary);
            fwrite(&(graphCSR->inverse_vertices->edges_idx[vertex_id]), sizeof (graphCSR->inverse_vertices->edges_idx[vertex_id]), 1, pBinary);
        }
#endif
    }

    for(vertex_id = 0; vertex_id < graphCSR->num_edges ; vertex_id++)
    {

        fwrite(&(graphCSR->sorted_edges_array->edges_array_src[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_src[vertex_id]), 1, pBinary);
        fwrite(&(graphCSR->sorted_edges_array->edges_array_dest[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_dest[vertex_id]), 1, pBinary);

#if WEIGHTED
        fwrite(&(graphCSR->sorted_edges_array->edges_array_weight[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_weight[vertex_id]), 1, pBinary);
#endif

#if DIRECTED
        if(graphCSR->inverse_vertices)
        {
            fwrite(&(graphCSR->inverse_sorted_edges_array->edges_array_src[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_src[vertex_id]), 1, pBinary);
            fwrite(&(graphCSR->inverse_sorted_edges_array->edges_array_dest[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_dest[vertex_id]), 1, pBinary);
#if WEIGHTED
            fwrite(&(graphCSR->inverse_sorted_edges_array->edges_array_weight[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_weight[vertex_id]), 1, pBinary);
#endif
        }
#endif
    }

    fclose(pBinary);


}


struct GraphCSR *readFromBinFileGraphCSR (const char *fname)
{

    FILE  *pBinary;
    uint32_t vertex_id;
    uint32_t num_vertices;
    uint32_t num_edges;
#if WEIGHTED
    float max_weight;
#endif
    size_t ret;

    pBinary = fopen(fname, "rb");
    if (pBinary == NULL)
    {
        err(1, "open: %s", fname);
        return NULL;
    }


    ret = fread(&num_edges, sizeof(num_edges), 1, pBinary);
    ret = fread(&num_vertices, sizeof(num_vertices), 1, pBinary);
#if WEIGHTED
    ret = fread(&max_weight, sizeof(max_weight), 1, pBinary);
#endif


#if DIRECTED
    struct GraphCSR *graphCSR = graphCSRNew(num_vertices, num_edges, 1);
#else
    struct GraphCSR *graphCSR = graphCSRNew(num_vertices, num_edges, 0);
#endif


#if WEIGHTED
    graphCSR->max_weight = max_weight;
#endif

    struct EdgeList *sorted_edges_array = newEdgeList(num_edges);
    sorted_edges_array->num_vertices = num_vertices;
#if WEIGHTED
    sorted_edges_array->max_weight = max_weight;
#endif

    graphCSR->sorted_edges_array = sorted_edges_array;

#if DIRECTED
    struct EdgeList *inverse_sorted_edges_array = newEdgeList(num_edges);
    inverse_sorted_edges_array->num_vertices = num_vertices;
#if WEIGHTED
    inverse_sorted_edges_array->max_weight = max_weight;
#endif
    graphCSR->inverse_sorted_edges_array = inverse_sorted_edges_array;
#endif

    for(vertex_id = 0; vertex_id < graphCSR->num_vertices ; vertex_id++)
    {

        ret = fread(&(graphCSR->vertices->out_degree[vertex_id]), sizeof (graphCSR->vertices->out_degree[vertex_id]), 1, pBinary);
        ret = fread(&(graphCSR->vertices->in_degree[vertex_id]), sizeof (graphCSR->vertices->in_degree[vertex_id]), 1, pBinary);
        ret = fread(&(graphCSR->vertices->edges_idx[vertex_id]), sizeof (graphCSR->vertices->edges_idx[vertex_id]), 1, pBinary);

#if DIRECTED
        if(graphCSR->inverse_vertices)
        {
            ret = fread(&(graphCSR->inverse_vertices->out_degree[vertex_id]), sizeof (graphCSR->inverse_vertices->out_degree[vertex_id]), 1, pBinary);
            ret = fread(&(graphCSR->inverse_vertices->in_degree[vertex_id]), sizeof (graphCSR->inverse_vertices->in_degree[vertex_id]), 1, pBinary);
            ret = fread(&(graphCSR->inverse_vertices->edges_idx[vertex_id]), sizeof (graphCSR->inverse_vertices->edges_idx[vertex_id]), 1, pBinary);
        }
#endif
    }

    for(vertex_id = 0; vertex_id < graphCSR->num_edges ; vertex_id++)
    {

        ret = fread(&(graphCSR->sorted_edges_array->edges_array_src[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_src[vertex_id]), 1, pBinary);
        ret = fread(&(graphCSR->sorted_edges_array->edges_array_dest[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_dest[vertex_id]), 1, pBinary);

#if WEIGHTED
        ret = fread(&(graphCSR->sorted_edges_array->edges_array_weight[vertex_id]), sizeof (graphCSR->sorted_edges_array->edges_array_weight[vertex_id]), 1, pBinary);
#endif

#if DIRECTED
        if(graphCSR->inverse_vertices)
        {
            ret = fread(&(graphCSR->inverse_sorted_edges_array->edges_array_src[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_src[vertex_id]), 1, pBinary);
            ret = fread(&(graphCSR->inverse_sorted_edges_array->edges_array_dest[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_dest[vertex_id]), 1, pBinary);
#if WEIGHTED
            ret = fread(&(graphCSR->inverse_sorted_edges_array->edges_array_weight[vertex_id]), sizeof (graphCSR->inverse_sorted_edges_array->edges_array_weight[vertex_id]), 1, pBinary);
#endif
        }
#endif
    }




    if(ret)
    {
        graphCSRPrint(graphCSR);
    }

    fclose(pBinary);

    return graphCSR;

}

void graphCSRPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}

// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : reorder.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:35:52
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>


#include "timer.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "edgeList.h"
#include "fixedPoint.h"
#include "sortRun.h"
#include "quantization.h"

#include "graphCSR.h"
#include "reorder.h"
#include "epochReorder.h"

#include "pageRank.h"
#include "incrementalAggregation.h"

struct EdgeList *reorderGraphListPageRank(struct GraphCSR *graph)
{


    uint32_t v;
    double epsilon = 1e-6;
    uint32_t iterations = 100;
    struct PageRankStats  *stats = NULL;
    uint32_t *labelsInverse;
    uint32_t *labels;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    labelsInverse = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));
    struct EdgeList *edgeList = NULL;



    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting PageRank Reordering/Relabeling");
    printf(" -----------------------------------------------------\n");

    Start(timer);
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        labelsInverse[v] = v;
    }



    stats = pageRankDataDrivenPushGraphCSR(epsilon, iterations, graph);
    // stats = pageRankPulCacheAnalysisGraphCSR(epsilon, iterations, graph);


    // make sure that nodes with no in/out degrees have zero scores
    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        if(!(graph->vertices[v].out_degree || graph->vertices[v].in_degree))
        {
            stats->pageRanks[v] = 0;
        }
    }

    labelsInverse = radixSortEdgesByPageRank(stats->pageRanks, labelsInverse, graph->num_vertices);

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {

        labels[labelsInverse[v]] = v;
    }


    edgeList = graph->sorted_edges_array;

    edgeList = relabelEdgeList(edgeList, labels);

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "PageRank Reordering/Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(labelsInverse);
    free(labels);

    return edgeList;
}


struct EdgeList *reorderGraphListEpochPageRank(struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t *labelsInverse;
    uint32_t *labels;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    struct EdgeList *edgeList = NULL;
    labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting PageRank Epoch Reordering/Relabeling");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    labelsInverse = epochReorderPageRank(graph);

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        labels[labelsInverse[v]] = v;
    }

    edgeList = graph->sorted_edges_array;

    edgeList = relabelEdgeList(edgeList, labels);

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "PageRank Epoch Reordering/Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(labelsInverse);
    free(labels);

    return edgeList;
}

struct EdgeList *reorderGraphListEpochBFS(struct GraphCSR *graph)
{

    uint32_t v;
    uint32_t *labelsInverse;
    uint32_t *labels;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    struct EdgeList *edgeList = NULL;
    labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting BFS Epoch Reordering/Relabeling");
    printf(" -----------------------------------------------------\n");

    Start(timer);


    labelsInverse = epochReorderRecordBFS(graph);

    #pragma omp parallel for
    for(v = 0; v < graph->num_vertices; v++)
    {
        labels[labelsInverse[v]] = v;
    }


    edgeList = graph->sorted_edges_array;

    edgeList = relabelEdgeList(edgeList, labels);

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "BFS Epoch Reordering/Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    free(labelsInverse);
    free(labels);

    return edgeList;
}

struct EdgeList *reorderGraphListEpochRabbit(struct GraphCSR *graph)
{

    // uint32_t v;
    // uint32_t *labels;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    struct EdgeList *edgeList = NULL;
    // labels = (uint32_t *) my_malloc(graph->num_vertices * sizeof(uint32_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting RABBIT Reordering/Relabeling");
    printf(" -----------------------------------------------------\n");

    Start(timer);

    struct IncrementalAggregationStats *stats;

    stats = incrementalAggregationGraphCSR(graph);

    // #pragma omp parallel for
    // for(v = 0; v < graph->num_vertices; v++)
    // {
    //     labels[labelsInverse[v]] = v;
    // }


    edgeList = graph->sorted_edges_array;

    edgeList = relabelEdgeList(edgeList, stats->labels);

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "RABBIT Reordering/Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    free(timer);
    freeIncrementalAggregationStats(stats);

    return edgeList;
}


void radixSortCountSortEdgesByRanks (uint32_t **pageRanksFP, uint32_t **pageRanksFPTemp, uint32_t **labels, uint32_t **labelsTemp, uint32_t radix, uint32_t buckets, uint32_t *buckets_count, uint32_t num_vertices)
{

    uint32_t *tempPointer1 = NULL;
    uint32_t *tempPointer2 = NULL;
    uint32_t t = 0;
    uint32_t o = 0;
    uint32_t u = 0;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t P = omp_get_max_threads();  // 32/8 8 bit radix needs 4 iterations
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t base = 0;

    #pragma omp parallel default(none) num_threads(P) shared(pageRanksFP, pageRanksFPTemp,radix,labels,labelsTemp,buckets,buckets_count, num_vertices) firstprivate(t_id, P, offset_end,offset_start,base,i,j,t,u,o)
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id * (num_vertices / P);


        if(t_id == (P - 1))
        {
            offset_end = offset_start + (num_vertices / P) + (num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (num_vertices / P);
        }


        //HISTOGRAM-KEYS
        for(i = 0; i < buckets; i++)
        {
            buckets_count[(t_id * buckets) + i] = 0;
        }


        for (i = offset_start; i < offset_end; i++)
        {
            u = (*pageRanksFP)[i];
            t = (u >> (radix * 8)) & 0xff;
            buckets_count[(t_id * buckets) + t]++;
        }


        #pragma omp barrier


        // SCAN BUCKETS
        if(t_id == 0)
        {
            for(i = 0; i < buckets; i++)
            {
                for(j = 0 ; j < P; j++)
                {
                    t = buckets_count[(j * buckets) + i];
                    buckets_count[(j * buckets) + i] = base;
                    base += t;
                }
            }
        }


        #pragma omp barrier

        //RANK-AND-PERMUTE
        for (i = offset_start; i < offset_end; i++)         /* radix sort */
        {
            u = (*pageRanksFP)[i];
            t = (u >> (radix * 8)) & 0xff;
            o = buckets_count[(t_id * buckets) + t];
            (*pageRanksFPTemp)[o] = (*pageRanksFP)[i];
            (*labelsTemp)[o] = (*labels)[i];
            buckets_count[(t_id * buckets) + t]++;

        }

    }

    tempPointer1 = *labels;
    *labels = *labelsTemp;
    *labelsTemp = tempPointer1;


    tempPointer2 = *pageRanksFP;
    *pageRanksFP = *pageRanksFPTemp;
    *pageRanksFPTemp = tempPointer2;

}

uint32_t *radixSortEdgesByPageRank (float *pageRanks, uint32_t *labels, uint32_t num_vertices)
{


    // printf("*** START Radix Sort Edges By Source *** \n");

    // struct Graph* graph = graphNew(edgeList->num_vertices, edgeList->num_edges, inverse);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    uint32_t v;
    uint32_t radix = 4;  // 32/8 8 bit radix needs 4 iterations
    uint32_t P = omp_get_max_threads();  // 32/8 8 bit radix needs 4 iterations
    uint32_t buckets = 256; // 2^radix = 256 buckets
    uint32_t *buckets_count = NULL;

    // omp_set_num_threads(P);

    uint32_t j = 0; //1,2,3 iteration

    uint32_t *pageRanksFP = NULL;
    uint32_t *pageRanksFPTemp = NULL;
    uint32_t *labelsTemp = NULL;

    buckets_count = (uint32_t *) my_malloc(P * buckets * sizeof(uint32_t));
    pageRanksFP = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    pageRanksFPTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    labelsTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < num_vertices; v++)
    {
        pageRanksFP[v] = FLOAT_2_U(*(uint32_t *)&pageRanks[v]);
    }

    for(j = 0 ; j < radix ; j++)
    {
        radixSortCountSortEdgesByRanks (&pageRanksFP, &pageRanksFPTemp, &labels, &labelsTemp, j, buckets, buckets_count, num_vertices);
    }


    free(buckets_count);
    free(pageRanksFP);
    free(pageRanksFPTemp);
    free(labelsTemp);

    //  for(v = 0; v < num_vertices; v++)
    // {
    //     printf("rank %u label %u pr %.22f \n",v, labelsInternal[v], pageRanks[labelsInternal[v]]);
    // }

    return labels;

}

uint32_t *radixSortEdgesByDegree (uint32_t *degrees, uint32_t *labels, uint32_t num_vertices)
{


    // printf("*** START Radix Sort Edges By Source *** \n");

    // struct Graph* graph = graphNew(edgeList->num_vertices, edgeList->num_edges, inverse);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    uint32_t radix = 4;  // 32/8 8 bit radix needs 4 iterations
    uint32_t P = omp_get_max_threads();  // 32/8 8 bit radix needs 4 iterations
    uint32_t buckets = 256; // 2^radix = 256 buckets
    uint32_t *buckets_count = NULL;

    // omp_set_num_threads(P);

    uint32_t j = 0; //1,2,3 iteration
    uint32_t *degreesTemp = NULL;
    uint32_t *labelsTemp = NULL;

    buckets_count = (uint32_t *) my_malloc(P * buckets * sizeof(uint32_t));
    degreesTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    labelsTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));


    for(j = 0 ; j < radix ; j++)
    {
        radixSortCountSortEdgesByRanks (&degrees, &degreesTemp, &labels, &labelsTemp, j, buckets, buckets_count, num_vertices);
    }


    free(buckets_count);
    free(degreesTemp);
    free(labelsTemp);

    return labels;

}

struct EdgeList *reorderGraphProcessPageRank(struct EdgeList *edgeList, struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));



    // Start(timer);
    edgeList = sortRunAlgorithms(edgeList, 1);
    // edgeList = radixSortEdgesBySourceOptimized(edgeList);
    // edgeListPrint(edgeList);
    // Stop(timer);
    // graphCSRPrintMessageWithtime("Radix Sort Edges By Source (Seconds)",Seconds(timer));

    // edgeListPrint(edgeList);
    if(arguments->dflag)
    {
        Start(timer);
        edgeList = removeDulpicatesSelfLoopEdges(edgeList);
        Stop(timer);
        graphCSRPrintMessageWithtime("Removing duplicate edges (Seconds)", Seconds(timer));
    }
    // edgeListPrint(edgeList);

    if(arguments->lmode == 5)
    {
        arguments->lmode = (5 + 3);
        edgeList = reorderGraphProcess(edgeList, arguments);
    }
    // edgeListPrint(edgeList);
    arguments->lmode = arguments->lmode - 3;
    edgeList = sortRunAlgorithms(edgeList, 1);
    // edgeListPrint(edgeList);

#if DIRECTED
    struct GraphCSR *graph = graphCSRNew(edgeList->num_vertices, edgeList->num_edges, 1);
#else
    struct GraphCSR *graph = graphCSRNew(edgeList->num_vertices, edgeList->num_edges, 0);
#endif

    Start(timer);
    graph = graphCSRAssignEdgeList (graph, edgeList, 0);
    Stop(timer);


    graphCSRPrintMessageWithtime("Process In/Out degrees of Nodes (Seconds)", Seconds(timer));

#if DIRECTED

    Start(timer);
    // struct EdgeList* inverse_edgeList = readEdgeListsbin(fnameb,1);
    struct EdgeList *inverse_edgeList = readEdgeListsMem(edgeList, 1, 0, 0);
    Stop(timer);
    // edgeListPrint(inverse_edgeList);
    graphCSRPrintMessageWithtime("Read Inverse Edge List From File (Seconds)", Seconds(timer));


    // Start(timer);
    inverse_edgeList = sortRunAlgorithms(inverse_edgeList, 1);
    // inverse_edgeList = radixSortEdgesBySourceOptimized(inverse_edgeList);
    // Stop(timer);
    // graphCSRPrintMessageWithtime("Radix Sort Inverse Edges By Source (Seconds)",Seconds(timer));

    Start(timer);
    graph = graphCSRAssignEdgeList (graph, inverse_edgeList, 1);
    Stop(timer);
    graphCSRPrintMessageWithtime("Process In/Out degrees of Inverse Nodes (Seconds)", Seconds(timer));

#endif

    if(arguments->lmode == 1) // pageRank
        edgeList =  reorderGraphListPageRank(graph);
    else if(arguments->lmode == 5) //epoch RABBIT
        edgeList =  reorderGraphListEpochRabbit(graph);
    else if(arguments->lmode == 6) //epoch pagerank
        edgeList = reorderGraphListEpochPageRank(graph); // in-degree
    else if(arguments->lmode == 7) //epoch BFS
        edgeList = reorderGraphListEpochBFS(graph); // in-degree


    if(graph->vertices)
        freeVertexArray(graph->vertices);
    // if(graph->sorted_edges_array)
    //   freeEdgeArray(graph->sorted_edges_array);
#if DIRECTED
    if(graph->inverse_vertices)
        freeVertexArray(graph->inverse_vertices);
    if(graph->inverse_sorted_edges_array)
        freeEdgeList(graph->inverse_sorted_edges_array);
#endif


    free(timer);

    return edgeList;


}


struct EdgeList *reorderGraphProcessDegree( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode)
{


    uint32_t *degrees;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));


    degrees = reorderGraphProcessInOutDegrees( degrees, edgeList, lmode);

    edgeList = reorderGraphListDegree( edgeList, degrees, lmode);

    return edgeList;

}

uint32_t reorderGraphProcessVertexSize( struct EdgeList *edgeList)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;
    uint32_t num_vertices = 0;

    #pragma omp parallel for default(none) private(i,src,dest) shared(edgeList) reduction(max: num_vertices)
    for(i = 0; i < edgeList->num_edges; i++)
    {

        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];
        num_vertices = maxTwoIntegers(num_vertices, maxTwoIntegers(src, dest));

    }

    return num_vertices;
}


uint32_t *reorderGraphProcessInOutDegrees(uint32_t *degrees, struct EdgeList *edgeList, uint32_t lmode)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;

    #pragma omp parallel for default(none) private(i,src,dest) shared(edgeList,degrees,lmode)
    for(i = 0; i < edgeList->num_edges; i++)
    {
        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        if(lmode == 1)
        {
            #pragma omp atomic update
            degrees[src]++;
        }
        else if(lmode == 2)
        {
            #pragma omp atomic update
            degrees[dest]++;
        }
        else if(lmode == 3)
        {
            #pragma omp atomic update
            degrees[dest]++;
            #pragma omp atomic update
            degrees[src]++;
        }

    }

    return degrees;
}



struct EdgeList *reorderGraphProcess(struct EdgeList *edgeList, struct Arguments *arguments)
{




    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    // printf("Filename : %s \n",fnameb);

    printf(" *****************************************************\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Reorder Process");
    printf(" -----------------------------------------------------\n");
    Start(timer);

    switch(arguments->lmode)
    {
    case 1 :
        edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// in-degree
        break;
    case 2  :
        edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 3  :
        edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// in/out-degree
        break;
    case 4  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 5  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 6  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 7  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 8  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 9  :
        // edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degre
        break;
    case 10 :
        edgeList = relabelEdgeListFromFile(edgeList, arguments->fnameb, edgeList->num_vertices);// load from file
        break;

    default :
        edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// out-degree
    }

    Stop(timer);


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Total Reorder Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");
    printf(" *****************************************************\n");

    free(timer);

    return edgeList;

}


struct EdgeList *reorderGraphListDegree(struct EdgeList *edgeList, uint32_t *degrees, uint32_t lmode)
{

    uint32_t v;
    uint32_t *labelsInverse;
    uint32_t *labels;
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    labelsInverse = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Starting Degree Reordering/Relabeling");
    printf(" -----------------------------------------------------\n");
    if(lmode == 1)  // in-degree
    {
        printf("| %-51s | \n", "OUT-DEGREE");
    }
    else if(lmode == 2)
    {
        printf("| %-51s | \n", "IN-DEGREE");
    }
    else if(lmode == 3)
    {
        printf("| %-51s | \n", "(IN+OUT-DEGREE)");
    }
    printf(" -----------------------------------------------------\n");

    Start(timer);

    #pragma omp parallel for
    for(v = 0; v < edgeList->num_vertices; v++)
    {
        labelsInverse[v] = v;
    }

    labelsInverse = radixSortEdgesByDegree(degrees, labelsInverse, edgeList->num_vertices);

    #pragma omp parallel for
    for(v = 0; v < edgeList->num_vertices; v++)
    {
        labels[labelsInverse[v]] = edgeList->num_vertices - 1 - v;
    }

    edgeList = relabelEdgeList(edgeList, labels);

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Degree Reordering/Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");


    free(timer);
    free(labelsInverse);
    free(labels);

    return edgeList;
}

struct EdgeList *relabelEdgeList(struct EdgeList *edgeList, uint32_t *labels)
{

    uint32_t i;

    #pragma omp parallel for
    for(i = 0; i < edgeList->num_edges; i++)
    {
        uint32_t src;
        uint32_t dest;
        src = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        edgeList->edges_array_src[i] = labels[src];
        edgeList->edges_array_dest[i] = labels[dest];

    }



    return edgeList;

}


struct EdgeList *relabelEdgeListFromFile(struct EdgeList *edgeList, const char *fnameb, uint32_t size)
{

    FILE *pText;
    uint32_t i;
    uint32_t dest = 0;
    uint32_t x = 0;

    uint32_t *labels;

    labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));


    char *fname_txt = (char *) malloc((strlen(fnameb) + 10) * sizeof(char));

    fname_txt = strcpy (fname_txt, fnameb);
    fname_txt = strcat (fname_txt, ".labels");

    printf("%s\n", fname_txt );
    pText = fopen(fname_txt, "r");

    if (pText == NULL)
    {
        return NULL;
    }

    while (1)
    {

        i = fscanf(pText, "%u\n", &dest);
        labels[x] = dest;
        x++;

        if( i == EOF )
            break;

    }
    fclose(pText);



    edgeList = relabelEdgeList(edgeList, labels);


    free(labels);
    free(fname_txt);

    return edgeList;
}


void writeLabelsToFile(const char *fnameb, uint32_t *labels, uint32_t size)
{

    FILE *fptr;
    uint32_t x;
    fptr = fopen(fnameb, "w");
    for(x = 0; x < size; x++)
    {
        fprintf(fptr, "%u %u\n", x, labels[x]);
    }

    fclose(fptr);


}

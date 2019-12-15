// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : graphStats.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <linux/types.h>
#include <string.h>
#include <omp.h>

#include "timer.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "graphCSR.h"
#include "graphStats.h"








void collectStats(struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
    // printf("Filename : %s \n",fnameb);

    printf(" *****************************************************\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Collect Stats Process");
    printf(" -----------------------------------------------------\n");
    Start(timer);

    struct GraphCSR *graphStats = graphCSRPreProcessingStep (arguments);


    __u32 *histogram_in = (__u32 *) my_malloc(sizeof(__u32) * arguments->binSize);
    __u32 *histogram_out = (__u32 *) my_malloc(sizeof(__u32) * arguments->binSize);


    __u32 i = 0;
    #pragma omp parallel for
    for(i = 0 ; i < arguments->binSize; i++)
    {
        histogram_in[i] = 0;
        histogram_out[i] = 0;
    }

    char *fname_txt = (char *) malloc((strlen(arguments->fnameb) + 20) * sizeof(char));
    char *fname_stats_out = (char *) malloc((strlen(arguments->fnameb) + 20) * sizeof(char));
    char *fname_stats_in = (char *) malloc((strlen(arguments->fnameb) + 20) * sizeof(char));
    char *fname_adjMat = (char *) malloc((strlen(arguments->fnameb) + 20) * sizeof(char));


    fname_txt = strcpy (fname_txt, arguments->fnameb);
    fname_adjMat = strcpy (fname_adjMat, arguments->fnameb);


    fname_adjMat  = strcat (fname_adjMat, ".bin-adj-SM.dat");// out-degree

    if(arguments->lmode == 1)
    {
        fname_stats_in = strcat (fname_txt, ".in-degree.dat");// in-degree
        countHistogram(graphStats, histogram_in, arguments->binSize, arguments->inout_degree);
        printHistogram(fname_stats_in, histogram_in, arguments->binSize);
    }
    else if(arguments->lmode == 2)
    {
        fname_stats_out = strcat (fname_txt, ".out-degree.dat");// out-degree
        countHistogram(graphStats, histogram_out, arguments->binSize, arguments->inout_degree);
        printHistogram(fname_stats_out, histogram_out, arguments->binSize);
    }


    printSparseMatrixList(fname_adjMat,  graphStats, arguments->binSize);


    Stop(timer);


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Collect Stats Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");
    printf(" *****************************************************\n");

    free(timer);
    graphCSRFree(graphStats);
    free(histogram_in);
    free(histogram_out);
    free(fname_txt);
    free(fname_stats_out);
    free(fname_stats_in);
    free(fname_adjMat);

}


void countHistogram(struct GraphCSR *graphStats, __u32 *histogram, __u32 binSize, __u32 inout_degree)
{

    __u32 v;
    __u32 index;

    #pragma omp parallel for
    for(v = 0; v < graphStats->num_vertices; v++)
    {

        index = v / ((graphStats->num_vertices / binSize) + 1);

        if(inout_degree == 1)
        {
            #pragma omp atomic update
            histogram[index] += graphStats->vertices->in_degree[v];
        }
        else if(inout_degree == 2)
        {
            #pragma omp atomic update
            histogram[index] += graphStats->vertices->out_degree[v];
        }
    }

}


void printHistogram(const char *fname_stats, __u32 *histogram, __u32 binSize)
{

    __u32 index;
    FILE *fptr;
    fptr = fopen(fname_stats, "w");
    for(index = 0; index < binSize; index++)
    {
        fprintf(fptr, "%u %u \n", index, histogram[index]);
    }
    fclose(fptr);
}


void printSparseMatrixList(const char *fname_stats, struct GraphCSR *graphStats, __u32 binSize)
{


    __u32 *SparseMatrix = (__u32 *) my_malloc(sizeof(__u32) * binSize * binSize);


    __u32 x;
    __u32 y;
    #pragma omp parallel for private(y) shared(SparseMatrix)
    for(x = 0; x < binSize; x++)
    {
        for(y = 0; y < binSize; y++)
        {
            SparseMatrix[(binSize * y) + x] = 0;
        }
    }


    __u32 i;

    #pragma omp parallel for
    for(i = 0; i < graphStats->num_edges; i++)
    {
        __u32 src;
        __u32 dest;
        src = graphStats->sorted_edges_array->edges_array_src[i] / ((graphStats->num_vertices / binSize) + 1);
        dest = graphStats->sorted_edges_array->edges_array_dest[i] / ((graphStats->num_vertices / binSize) + 1);

        #pragma omp atomic update
        SparseMatrix[(binSize * dest) + src]++;

    }

    FILE *fptr;
    fptr = fopen(fname_stats, "w");
    for(x = 0; x < binSize; x++)
    {
        for(y = 0; y < binSize; y++)
        {
            fprintf(fptr, "%u %u %u\n", x, y, SparseMatrix[(binSize * y) + x]);
        }
    }

    fclose(fptr);
    free(SparseMatrix);

}


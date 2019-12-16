// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : grid.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <omp.h>

#include "grid.h"
#include "edgeList.h"
#include "vertex.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "bitmap.h"
#include "timer.h"

void gridPrint(struct Grid *grid)
{


    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Grid Properties");
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
    printf("| %-51u | \n", grid->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Edges (E)");
    printf("| %-51u | \n", grid->num_edges);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Partitions (P)");
    printf("| %-51u | \n", grid->num_partitions);
    printf(" -----------------------------------------------------\n");

    // _u32 i;
    //  for ( i = 0; i < grid->num_vertices; ++i)
    //     {

    //         uint32_t begin = getPartitionRangeBegin();
    //         uint32_t end = getPartitionRangeEnd();



    //     }

    //   uint32_t i;
    //    for ( i = 0; i < (grid->num_partitions*grid->num_partitions); ++i)
    //       {

    //       uint32_t x = i % grid->num_partitions;    // % is the "modulo operator", the remainder of i / width;
    // uint32_t y = i / grid->num_partitions;


    //       printf("| %-11s (%u,%u)   | \n", "Partition: ", y, x);
    //          printf("| %-11s %-40u  | \n", "Edges: ", grid->partitions[i].num_edges);
    //          printf("| %-11s %-40u  | \n", "Vertices: ", grid->partitions[i].num_vertices);
    //          edgeListPrint(grid->partitions[i].edgeList);
    //       }


}

void   graphGridResetActivePartitions(struct Grid *grid)
{

    uint32_t totalPartitions = 0;
    totalPartitions = grid->num_partitions * grid->num_partitions;
    uint32_t i;

    #pragma omp parallel for default(none) shared(grid,totalPartitions) private(i)
    for (i = 0; i < totalPartitions; ++i)
    {
        grid->activePartitions[i] = 0;
    }



}

void   graphGridResetActivePartitionsMap(struct Grid *grid)
{

    clearBitmap(grid->activePartitionsMap);

}

void   graphGridSetActivePartitionsMap(struct Grid *grid, uint32_t vertex)
{

    uint32_t row = getPartitionID(grid->num_vertices, grid->num_partitions, vertex);
    uint32_t Partition_idx = 0;
    uint32_t i;
    uint32_t totalPartitions = 0;
    totalPartitions = grid->num_partitions;

    // #pragma omp parallel for default(none) shared(grid,totalPartitions,row) private(i,Partition_idx)

    for ( i = 0; i < totalPartitions; ++i)
    {

        Partition_idx = (row * totalPartitions) + i;

        if(grid->partitions[Partition_idx].edgeList->num_edges)
        {
            if(!getBit(grid->activePartitionsMap, Partition_idx))
            {
                setBitAtomic(grid->activePartitionsMap, Partition_idx);
            }
        }
    }
}

void   graphGridSetActivePartitions(struct Grid *grid, uint32_t vertex)
{

    uint32_t row = getPartitionID(grid->num_vertices, grid->num_partitions, vertex);
    uint32_t Partition_idx = 0;
    uint32_t i;
    uint32_t totalPartitions = 0;
    totalPartitions = grid->num_partitions;

    // #pragma omp parallel for default(none) shared(grid,totalPartitions,row) private(i,Partition_idx)
    for ( i = 0; i < totalPartitions; ++i)
    {

        Partition_idx = (row * totalPartitions) + i;
        if(grid->partitions[Partition_idx].edgeList->num_edges)
        {
            grid->activePartitions[Partition_idx] = 1;
        }
    }
}



struct Grid *gridNew(struct EdgeList *edgeList)
{


    uint32_t totalPartitions = 0;


    struct Grid *grid = (struct Grid *) my_malloc( sizeof(struct Grid));


    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    grid->num_edges = edgeList->num_edges;
    grid->num_vertices = edgeList->num_vertices;
    grid->num_partitions = gridCalculatePartitions(edgeList);
    totalPartitions = grid->num_partitions * grid->num_partitions;

    grid->partitions = (struct Partition *) my_malloc(totalPartitions * sizeof(struct Partition));
    grid->activePartitions = (uint32_t *) my_malloc(totalPartitions * sizeof(uint32_t));
    grid->out_degree = (uint32_t *) my_malloc(grid->num_vertices * sizeof(uint32_t));
    grid->in_degree = (uint32_t *) my_malloc(grid->num_vertices * sizeof(uint32_t));

    // grid->activeVertices = newBitmap(grid->num_vertices);
    grid->activePartitionsMap = newBitmap(totalPartitions);

    uint32_t i;
    #pragma omp parallel for default(none) private(i) shared(totalPartitions,grid)
    for (i = 0; i < totalPartitions; ++i)
    {

        grid->partitions[i].num_edges = 0;
        grid->partitions[i].num_vertices = 0;   /* code */
        grid->activePartitions[i] = 0;


    }


    #pragma omp parallel for default(none) private(i) shared(grid)
    for (i = 0; i < grid->num_vertices ; ++i)
    {

        grid->out_degree[i] = 0;
        grid->in_degree[i] = 0;

    }


    Start(timer);
    grid = graphGridProcessInOutDegrees(grid, edgeList);
    Stop(timer);
    gridPrintMessageWithtime("Grid Process In Out Degrees (Seconds)", Seconds(timer));


    Start(timer);
    grid = gridPartitionEdgeListSizePreprocessing(grid, edgeList);
    Stop(timer);
    gridPrintMessageWithtime("Partition EdgeList Size (Seconds)", Seconds(timer));


    Start(timer);
    grid = gridPartitionsMemoryAllocations(grid);
    Stop(timer);
    gridPrintMessageWithtime("Partitions Memory Allocations (Seconds)", Seconds(timer));


    Start(timer);
    grid = gridPartitionEdgePopulation(grid, edgeList);
    Stop(timer);
    gridPrintMessageWithtime("Partition Edge Population (Seconds)", Seconds(timer));


    Start(timer);
    grid = gridPartitionVertexSizePreprocessing(grid);
    Stop(timer);
    gridPrintMessageWithtime("Partition Vertex Size (Seconds)", Seconds(timer));

    return grid;
}



void  gridFree(struct Grid *grid)
{
    

    if(grid)
    {
        uint32_t totalPartitions = grid->num_partitions * grid->num_partitions;
        uint32_t i;

        for (i = 0; i < totalPartitions; ++i)
        {

            freeEdgeList(grid->partitions[i].edgeList);
        }

        freeBitmap(grid->activePartitionsMap);

        if(grid->activePartitions)
            free(grid->activePartitions);

        if(grid->out_degree)
            free(grid->out_degree);

        if(grid->in_degree)
            free(grid->in_degree);

        if(grid->partitions)
            free(grid->partitions);

        free(grid);
    }
}


struct Grid *graphGridProcessInOutDegrees(struct Grid *grid, struct EdgeList *edgeList)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;

    #pragma omp parallel for default(none) private(i,src,dest) shared(edgeList,grid)
    for(i = 0; i < edgeList->num_edges; i++)
    {

        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        #pragma omp atomic update
        grid->out_degree[src]++;

        #pragma omp atomic update
        grid->in_degree[dest]++;

    }

    return grid;

}

struct Grid *gridPartitionVertexSizePreprocessing(struct Grid *grid)
{

    uint32_t i;
    uint32_t j;
    uint32_t src;
    uint32_t dest;
    uint32_t num_vertices = 0;
    uint32_t totalPartitions = grid->num_partitions * grid->num_partitions;

    // #pragma omp parallel for default(none) private(i) shared(totalPartitions,grid)
    #pragma omp parallel for default(none) private(i,src,dest,num_vertices) shared(totalPartitions,grid) schedule(dynamic,1024)
    for ( j = 0; j < totalPartitions; ++j)
    {
        num_vertices = 0;
        // #pragma omp parallel for default(none) private(i,src,dest) shared(j,grid) schedule(dynamic,1024) reduction(max:num_vertices)
        for(i = 0; i <  grid->partitions[j].edgeList->num_edges; i++)
        {

            src  =  grid->partitions[j].edgeList->edges_array_src[i];
            dest =  grid->partitions[j].edgeList->edges_array_dest[i];

            num_vertices = maxTwoIntegers(num_vertices, maxTwoIntegers(src, dest));

        }

        grid->partitions[j].num_vertices = num_vertices;
        grid->partitions[j].edgeList->num_vertices = num_vertices;
    }

    return grid;


}

struct Grid *gridPartitionEdgeListSizePreprocessing(struct Grid *grid, struct EdgeList *edgeList)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;
    uint32_t Partition_idx;

    uint32_t num_partitions = grid->num_partitions;
    uint32_t num_vertices = grid->num_vertices;


    uint32_t row;
    uint32_t col;

    #pragma omp parallel for default(none) private(i,row,col,src,dest,Partition_idx) shared(num_vertices, num_partitions,edgeList,grid)
    for(i = 0; i < edgeList->num_edges; i++)
    {

        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        // __sync_fetch_and_add(&grid->out_degree[src],1);
        // __sync_fetch_and_add(&grid->in_degree[dest],1);

        row = getPartitionID(num_vertices, num_partitions, src);
        col = getPartitionID(num_vertices, num_partitions, dest);
        Partition_idx = (row * num_partitions) + col;

        // __sync_fetch_and_add(&grid->partitions[Partition_idx].num_edges,1);

        #pragma omp atomic update
        grid->partitions[Partition_idx].num_edges++;

    }

    return grid;

}


struct Grid *gridPartitionEdgePopulation(struct Grid *grid, struct EdgeList *edgeList)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;
    uint32_t Partition_idx;
    uint32_t Edge_idx;

    uint32_t num_partitions = grid->num_partitions;
    uint32_t num_vertices = grid->num_vertices;

    uint32_t row;
    uint32_t col;




    #pragma omp parallel for default(none) private(Edge_idx,i,row,col,src,dest,Partition_idx) shared(num_vertices, num_partitions,edgeList,grid)
    for(i = 0; i < edgeList->num_edges; i++)
    {


        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        row = getPartitionID(num_vertices, num_partitions, src);
        col = getPartitionID(num_vertices, num_partitions, dest);
        Partition_idx = (row * num_partitions) + col;

        Edge_idx = __sync_fetch_and_add(&grid->partitions[Partition_idx].num_edges, 1);

        grid->partitions[Partition_idx].edgeList->edges_array_src[Edge_idx] = edgeList->edges_array_src[i];
        grid->partitions[Partition_idx].edgeList->edges_array_dest[Edge_idx] = edgeList->edges_array_dest[i];
#if WEIGHTED
        grid->partitions[Partition_idx].edgeList->edges_array_weight[Edge_idx] = edgeList->edges_array_weight[i];
#endif
    }



    return grid;

}


struct Grid *gridPartitionsMemoryAllocations(struct Grid *grid)
{

    uint32_t i;
    uint32_t totalPartitions = grid->num_partitions * grid->num_partitions;

    #pragma omp parallel for default(none) private(i) shared(totalPartitions,grid)
    for ( i = 0; i < totalPartitions; ++i)
    {

        grid->partitions[i].edgeList = newEdgeList(grid->partitions[i].num_edges);
        grid->partitions[i].edgeList->num_vertices = grid->partitions[i].num_vertices;
        grid->partitions[i].num_edges = 0;

    }

    return grid;


}

uint32_t gridCalculatePartitions(struct EdgeList *edgeList)
{
    //epfl everything graph
    uint32_t num_vertices  = edgeList->num_vertices;
    uint32_t num_Paritions = (num_vertices * 8 / 1024) / 20;
    if(num_Paritions > 512)
        num_Paritions = 256;
    if(num_Paritions == 0 )
        num_Paritions = 4;

    return num_Paritions;

}



inline uint32_t getPartitionID(uint32_t vertices, uint32_t partitions, uint32_t vertex_id)
{

    uint32_t partition_size = vertices / partitions;

    if (vertices % partitions == 0)
    {

        return vertex_id / partition_size;
    }

    partition_size += 1;

    uint32_t split_point = vertices % partitions * partition_size;

    return (vertex_id < split_point) ? vertex_id / partition_size : (vertex_id - split_point) / (partition_size - 1) + (vertices % partitions);
}


void gridPrintMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}
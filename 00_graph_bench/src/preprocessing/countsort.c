// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : countsort.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-11-09 10:34:42
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#include "countsort.h"
#include "edgeList.h"
#include "vertex.h"
#include "myMalloc.h"
#include "graphCSR.h"

struct EdgeList  *countSortEdgesBySource (struct EdgeList *edgeList)
{




    __u32 key = 0;
    __u32 pos = 0;
    __u32 num_vertices = edgeList->num_vertices;
    __u32 num_edges = edgeList->num_edges;
    __u32 i = 0;
    __u32 j = 0;
    __u32 P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    __u32 t_id = 0;
    __u32 offset_start = 0;
    __u32 offset_end = 0;
    __u32 base = 0;


    __u32 *vertex_count = (__u32 *) my_malloc( P * num_vertices * sizeof(__u32));



    struct EdgeList *sorted_edges_array = newEdgeList(num_edges);

    #pragma omp parallel default(none) shared(vertex_count,sorted_edges_array,edgeList,num_edges,num_vertices) firstprivate(t_id, P, offset_end,offset_start,base,i,j,key,pos)
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id * (num_edges / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (num_edges / P) + (num_edges % P) ;
        }
        else
        {
            offset_end = offset_start + (num_edges / P);
        }

        //HISTOGRAM-KEYS
        for(i = 0; i < num_vertices; i++)
        {
            vertex_count[(t_id * num_vertices) + i] = 0;
        }

        // count occurrence of key: id of the source vertex
        for(i = offset_start; i < offset_end; i++)
        {
            key = edgeList->edges_array_src[i];
            vertex_count[(t_id * num_vertices) + key]++;
        }

        #pragma omp barrier

        //SCAN BUCKETS
        if(t_id == 0)
        {
            for(i = 0; i < num_vertices; i++)
            {
                for(j = 0 ; j < P; j++)
                {
                    pos = vertex_count[(j * num_vertices) + i];
                    vertex_count[(j * num_vertices) + i] = base;
                    base += pos;
                }
            }
        }

        #pragma omp barrier

        //RANK-AND-PERMUTE

        for(i = offset_start; i < offset_end; i++)
        {

            key = edgeList->edges_array_src[i];
            pos = vertex_count[(t_id * num_vertices) + key];

            sorted_edges_array->edges_array_dest[pos] = edgeList->edges_array_dest[i];
            sorted_edges_array->edges_array_src[pos] = edgeList->edges_array_src[i];
#if WEIGHTED
            sorted_edges_array->edges_array_weight[pos] = edgeList->edges_array_weight[i];
#endif

            vertex_count[(t_id * num_vertices) + key]++;

        }

    }

    free(vertex_count);
    freeEdgeList(edgeList);

    edgeList = sorted_edges_array;

    return edgeList;

}


struct EdgeList *countSortEdgesByDestination (struct EdgeList *edgeList)
{

    __u32 key = 0;
    __u32 pos = 0;
    __u32 num_vertices = edgeList->num_vertices;
    __u32 num_edges = edgeList->num_edges;
    __u32 i = 0;
    __u32 j = 0;
    __u32 P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    __u32 t_id = 0;
    __u32 offset_start = 0;
    __u32 offset_end = 0;
    __u32 base = 0;


    __u32 *vertex_count = (__u32 *) my_malloc( P * num_vertices * sizeof(__u32));



    struct EdgeList *sorted_edges_array = newEdgeList(num_edges);

    #pragma omp parallel default(none) shared(vertex_count,sorted_edges_array,edgeList,num_edges,num_vertices) firstprivate(t_id, P, offset_end,offset_start,base,i,j,key,pos)
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id * (num_edges / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (num_edges / P) + (num_edges % P) ;
        }
        else
        {
            offset_end = offset_start + (num_edges / P);
        }

        //HISTOGRAM-KEYS
        for(i = 0; i < num_vertices; i++)
        {
            vertex_count[(t_id * num_vertices) + i] = 0;
        }

        // count occurrence of key: id of the source vertex
        for(i = offset_start; i < offset_end; i++)
        {
            key = edgeList->edges_array_dest[i];
            vertex_count[(t_id * num_vertices) + key]++;
        }

        #pragma omp barrier

        //SCAN BUCKETS
        if(t_id == 0)
        {
            for(i = 0; i < num_vertices; i++)
            {
                for(j = 0 ; j < P; j++)
                {
                    pos = vertex_count[(j * num_vertices) + i];
                    vertex_count[(j * num_vertices) + i] = base;
                    base += pos;
                }
            }
        }

        #pragma omp barrier

        //RANK-AND-PERMUTE

        for(i = offset_start; i < offset_end; i++)
        {

            key = edgeList->edges_array_dest[i];
            pos = vertex_count[(t_id * num_vertices) + key];

            sorted_edges_array->edges_array_dest[pos] = edgeList->edges_array_dest[i];
            sorted_edges_array->edges_array_src[pos] = edgeList->edges_array_src[i];
#if WEIGHTED
            sorted_edges_array->edges_array_weight[pos] = edgeList->edges_array_weight[i];
#endif

            vertex_count[(t_id * num_vertices) + key]++;

        }

    }

    free(vertex_count);
    freeEdgeList(edgeList);

    edgeList = sorted_edges_array;

    return edgeList;




}


struct EdgeList  *countSortEdgesBySourceAndDestination (struct EdgeList *edgeList)
{


    edgeList = countSortEdgesByDestination (edgeList);
    edgeList = countSortEdgesBySource (edgeList);

    return edgeList;
}







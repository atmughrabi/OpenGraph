// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : vertex.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <linux/types.h>
#include <omp.h>

#include "graphCSR.h"
#include "vertex.h"
#include "myMalloc.h"

struct Vertex *newVertexArray(__u32 num_vertices)
{


    struct Vertex *vertex_array = (struct Vertex *) my_malloc(sizeof(struct Vertex));

    vertex_array->out_degree = (__u32 *) my_malloc( num_vertices * sizeof(__u32));
    vertex_array->in_degree = (__u32 *) my_malloc( num_vertices * sizeof(__u32));
    vertex_array->edges_idx = (__u32 *) my_malloc( num_vertices * sizeof(__u32));

    __u32 i;

    for(i = 0; i < num_vertices; i++)
    {

        vertex_array->edges_idx[i]  = 0;
        vertex_array->out_degree[i] = 0;
        vertex_array->in_degree[i] = 0;

    }

    return vertex_array;

}

struct GraphCSR *mapVertices (struct GraphCSR *graph, __u8 inverse)
{

    __u32 i;
    __u32 vertex_id;

    struct Vertex *vertices;
    struct EdgeList *sorted_edges_array;


#if DIRECTED
    if(inverse)
    {
        sorted_edges_array = graph->inverse_sorted_edges_array;
        vertices = graph->inverse_vertices; // sorted edge array
    }
    else
    {
        sorted_edges_array = graph->sorted_edges_array;
        vertices = graph->vertices;
    }

#else
    sorted_edges_array = graph->sorted_edges_array;
    vertices = graph->vertices;
#endif


    vertex_id = sorted_edges_array->edges_array_src[0];
    vertices->edges_idx[vertex_id] = 0;

    for(i = 1; i < graph->num_edges; i++)
    {

        if(sorted_edges_array->edges_array_src[i] != sorted_edges_array->edges_array_src[i - 1])
        {

            vertex_id = sorted_edges_array->edges_array_src[i];
            vertices->edges_idx[vertex_id] = i;

        }
    }

    return graph;

}

void partitionEdgeListOffsetStartEnd(struct GraphCSR *graph, struct EdgeList *sorted_edges_array, __u32 *offset_start, __u32 *offset_end)
{

    __u32 i;
    __u32 j;
    __u32 P = numThreads;

    if(P >  graph->num_edges && graph->num_edges != 0)
        P = graph->num_edges;


    for(i = 0 ; i < P ; i++)
    {

        offset_start[i] = 0;
        offset_end[i] = 0;

    }

    offset_start[0] = 0;
    offset_end[0] = offset_start[0] + (graph->num_edges / P);



    if(1 == (P))
    {
        offset_end[0] = graph->num_edges;
    }

    for(i = 1 ; i < P ; i++)
    {

        j = offset_end[i - 1];

        if(j == graph->num_edges)
        {
            offset_start[i] = graph->num_edges;
            offset_end[i] = graph->num_edges;
            continue;
        }


        for(; j < graph->num_edges; j++)
        {

            if(sorted_edges_array->edges_array_src[j] != sorted_edges_array->edges_array_src[j - 1])
            {
                offset_start[i] = j;
                offset_end[i - 1] = j;

                if(i == (P - 1))
                {
                    offset_end[i] = i * (graph->num_edges / P) + (graph->num_edges / P) + (graph->num_edges % P) ;
                }
                else
                {
                    offset_end[i] =  offset_start[i] + (graph->num_edges / P);
                }

                if(offset_end[i] > graph->num_edges && offset_start[i] < graph->num_edges)
                {
                    offset_end[i] = graph->num_edges;
                    // printf("3-%u %u\n", offset_start[i], offset_end[i] );
                }

                break;
            }
            else if(sorted_edges_array->edges_array_src[j] == sorted_edges_array->edges_array_src[j - 1] && j == (graph->num_edges - 1))
            {
                offset_start[i] = graph->num_edges;
                offset_end[i] = graph->num_edges;
                offset_end[i - 1] = graph->num_edges;

            }

        }


    }
    // for(i=0 ; i < P ; i++){

    //    printf("%u %u\n", offset_start[i], offset_end[i] );

    // }


}

struct GraphCSR *mapVerticesWithInOutDegree (struct GraphCSR *graph, __u8 inverse)
{

    __u32 i;
    __u32 vertex_id;
    // __u32 vertex_id_dest;
    __u32 P = numThreads;
    struct Vertex *vertices;
    struct EdgeList *sorted_edges_array;

    __u32 *offset_start_arr = (__u32 *) my_malloc( P * sizeof(__u32));
    __u32 *offset_end_arr = (__u32 *) my_malloc( P * sizeof(__u32));

    // for(vertex_id = 0; vertex_id < graph->num_vertices; vertex_id++){

    //     printf("-->v %u out_degree %u\n",vertex_id, graph->vertices->out_degree[vertex_id] );
    // }


#if DIRECTED

    if(inverse)
    {
        sorted_edges_array = graph->inverse_sorted_edges_array;
        vertices = graph->inverse_vertices; // sorted edge array
    }
    else
    {
        sorted_edges_array = graph->sorted_edges_array;
        vertices = graph->vertices;
    }

#else
    sorted_edges_array = graph->sorted_edges_array;
    vertices = graph->vertices;
#endif

    //edge list must be sorted
    partitionEdgeListOffsetStartEnd(graph, sorted_edges_array, offset_start_arr, offset_end_arr);


    __u32 t_id = 0;
    __u32 offset_start = 0;
    __u32 offset_end = 0;



    #pragma omp parallel default(none) private(i,vertex_id) shared(graph,vertices,sorted_edges_array,offset_start_arr,offset_end_arr) firstprivate(t_id, offset_end,offset_start)
    {

        t_id = omp_get_thread_num();

        offset_start = offset_start_arr[t_id];
        offset_end = offset_end_arr[t_id];

        // printf("t_id %u start %u end %u \n",t_id,offset_start, offset_end);

        if(offset_start < graph->num_edges)
        {

            vertex_id = sorted_edges_array->edges_array_src[offset_start];
            vertices->edges_idx[vertex_id] = offset_start;
            vertices->out_degree[vertex_id]++;

            for(i = offset_start + 1; i < offset_end; i++)
            {

                vertex_id = sorted_edges_array->edges_array_src[i];
                vertices->out_degree[vertex_id]++;

                if(sorted_edges_array->edges_array_src[i] != sorted_edges_array->edges_array_src[i - 1])
                {
                    vertices->edges_idx[vertex_id] = i;
                }
            }



        }

    }

#if DIRECTED
    if(!inverse)
    {

        #pragma omp parallel for default(none) private(vertex_id) shared(vertices,graph)
        for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
        {
            graph->inverse_vertices->in_degree[vertex_id] = vertices->out_degree[vertex_id];
        }

    }
    else
    {
        #pragma omp parallel for default(none) private(vertex_id) shared(vertices,graph)
        for(vertex_id = 0; vertex_id < graph->num_vertices ; vertex_id++)
        {
            graph->vertices->in_degree[vertex_id] = vertices->out_degree[vertex_id];
        }

    }
#endif


    free(offset_start_arr);
    free(offset_end_arr);

    // for(vertex_id = 0; vertex_id < graph->num_vertices; vertex_id++){

    //     printf("<--v %u out_degree %u\n",vertex_id, graph->vertices->out_degree[vertex_id] );
    // }

    return graph;

}

void vertexArrayMaxOutdegree(struct Vertex *vertex_array, __u32 num_vertices)
{


    __u32 i;
    __u32 out_degree = 0;
    __u32 index = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Node", "*max_out_degree", "in_degree");
    printf(" -----------------------------------------------------\n");



    for(i = 0; i < num_vertices; i++)
    {


        out_degree = maxTwoIntegers(out_degree, vertex_array->out_degree[i]);
        if(vertex_array->out_degree[i] == out_degree)
            index = i;


    }


    printf("| %-15u | %-15u | %-15u | \n", index,  vertex_array->out_degree[index], vertex_array->in_degree[index]);
    printf(" -----------------------------------------------------\n");

}

void vertexArrayMaxInDegree(struct Vertex *vertex_array, __u32 num_vertices)
{


    __u32 i;
    __u32 in_degree = 0;
    __u32 index = 0;
    printf(" -----------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | \n", "Node", "out_degree", "*max_in_degree");
    printf(" -----------------------------------------------------\n");

    for(i = 0; i < num_vertices; i++)
    {


        in_degree = maxTwoIntegers(in_degree, vertex_array->out_degree[i]);
        if(vertex_array->out_degree[i] == in_degree)
            index = i;


    }


    printf("| %-15u | %-15u | %-15u | \n", index,  vertex_array->in_degree[index], vertex_array->out_degree[index]);
    printf(" -----------------------------------------------------\n");

}

void printVertexArray(struct Vertex *vertex_array, __u32 num_vertices)
{


    __u32 i;

    printf("| %-15s | %-15s | %-15s |\n", "Node", "out_degree", "in_degree");

    for(i = 0; i < num_vertices; i++)
    {

        if((vertex_array->out_degree[i] > 0) )
            printf("| %-15u | %-15u | %-15u | \n", i,  vertex_array->out_degree[i], vertex_array->in_degree[i]);

    }

}

void freeVertexArray(struct Vertex *vertices)
{
    if(vertices)
    {
        if(vertices->edges_idx)
            free(vertices->edges_idx);
        if(vertices->out_degree)
            free(vertices->out_degree);
        if(vertices->in_degree)
            free(vertices->in_degree);

        free(vertices);
    }
}


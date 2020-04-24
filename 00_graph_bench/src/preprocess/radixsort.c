// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : radixsort.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:35:52
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <omp.h>

#include "radixsort.h"
#include "edgeList.h"
#include "vertex.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "timer.h"

// A function to do counting sort of edgeList according to
// the digit represented by exp
// The parallel version has the following pseudo code
// parallel_for part in 0..K-1
//   for i in indexes(part)
//     bucket = compute_bucket(a[i])
//     Cnt[part][bucket]++

// base = 0
// for bucket in 0..R-1
//   for part in 0..K-1
//     Cnt[part][bucket] += base
//     base = Cnt[part][bucket]

// parallel_for part in 0..K-1
//   for i in indexes(part)
//     bucket = compute_bucket(a[i])
//     out[Cnt[part][bucket]++] = a[i]
void radixSortCountSortEdgesBySource (struct EdgeList** sorted_edges_array, struct EdgeList** edgeList, uint32_t radix, uint32_t buckets, uint32_t* buckets_count){

	struct EdgeList* temp_edges_array = NULL; 
    uint32_t num_edges = (*edgeList)->num_edges;
    uint32_t t = 0;
    uint32_t o = 0;
    uint32_t u = 0;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t base = 0;

    #pragma omp parallel default(none) shared(sorted_edges_array,edgeList,radix,buckets,buckets_count,num_edges) firstprivate(t_id, P, offset_end,offset_start,base,i,j,t,u,o) 
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id*(num_edges/P);


        if(t_id == (P-1)){
            offset_end = offset_start+(num_edges/P) + (num_edges%P) ;
        }
        else{
            offset_end = offset_start+(num_edges/P);
        }
        

        //HISTOGRAM-KEYS 
        for(i=0; i < buckets; i++){ 
            buckets_count[(t_id*buckets)+i] = 0;
        }

       
        for (i = offset_start; i < offset_end; i++) {      
            u = (*edgeList)->edges_array_src[i];
            t = (u >> (radix*8)) & 0xff;
            buckets_count[(t_id*buckets)+t]++;
        }


        #pragma omp barrier

       
        // SCAN BUCKETS
        if(t_id == 0){
            for(i=0; i < buckets; i++){
                 for(j=0 ; j < P; j++){
                 t = buckets_count[(j*buckets)+i];
                 buckets_count[(j*buckets)+i] = base;
                 base += t;
                }
            }
        }


        #pragma omp barrier

        //RANK-AND-PERMUTE
        for (i = offset_start; i < offset_end; i++) {       /* radix sort */
            u = (*edgeList)->edges_array_src[i];
            t = (u >> (radix*8)) & 0xff;
            o = buckets_count[(t_id*buckets)+t];
            (*sorted_edges_array)->edges_array_dest[o] = (*edgeList)->edges_array_dest[i];
            (*sorted_edges_array)->edges_array_src[o] = (*edgeList)->edges_array_src[i];        
            #if WEIGHTED
               (*sorted_edges_array)->edges_array_weight[o]= (*edgeList)->edges_array_weight[i];
            #endif
            buckets_count[(t_id*buckets)+t]++;

        }

    }

    temp_edges_array = *sorted_edges_array;
    *sorted_edges_array = *edgeList;
    *edgeList = temp_edges_array;
    
}

void radixSortCountSortEdgesByDestination (struct EdgeList** sorted_edges_array, struct EdgeList** edgeList, uint32_t radix, uint32_t buckets, uint32_t* buckets_count){

    struct EdgeList* temp_edges_array = NULL; 
    uint32_t num_edges = (*edgeList)->num_edges;
    uint32_t t = 0;
    uint32_t o = 0;
    uint32_t u = 0;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t base = 0;

    #pragma omp parallel default(none) shared(sorted_edges_array,edgeList,radix,buckets,buckets_count,num_edges) firstprivate(t_id, P, offset_end,offset_start,base,i,j,t,u,o) 
    {
        P = omp_get_num_threads();
        t_id = omp_get_thread_num();
        offset_start = t_id*(num_edges/P);


        if(t_id == (P-1)){
            offset_end = offset_start+(num_edges/P) + (num_edges%P) ;
        }
        else{
            offset_end = offset_start+(num_edges/P);
        }
        

        //HISTOGRAM-KEYS 
        for(i=0; i < buckets; i++){ 
            buckets_count[(t_id*buckets)+i] = 0;
        }

       
        for (i = offset_start; i < offset_end; i++) {      
            u = (*edgeList)->edges_array_dest[i];
            t = (u >> (radix*8)) & 0xff;
            buckets_count[(t_id*buckets)+t]++;
        }


        #pragma omp barrier


        //SCAN BUCKETS
        if(t_id == 0){

        for(i=0; i < buckets; i++){
             for(j=0 ; j < P; j++){
             t = buckets_count[(j*buckets)+i];
             buckets_count[(j*buckets)+i] = base;
             base += t;
         }
        }

        }

        #pragma omp barrier

        //RANK-AND-PERMUTE
        for (i = offset_start; i < offset_end; i++) {       /* radix sort */
            u = (*edgeList)->edges_array_dest[i];
            t = (u >> (radix*8)) & 0xff;
            o = buckets_count[(t_id*buckets)+t];
            (*sorted_edges_array)->edges_array_dest[o] = (*edgeList)->edges_array_dest[i];
            (*sorted_edges_array)->edges_array_src[o] = (*edgeList)->edges_array_src[i];        
            #if WEIGHTED
               (*sorted_edges_array)->edges_array_weight[o]= (*edgeList)->edges_array_weight[i];
            #endif
            buckets_count[(t_id*buckets)+t]++;

        }

    }

    temp_edges_array = *sorted_edges_array;
    *sorted_edges_array = *edgeList;
    *edgeList = temp_edges_array;
    
}

// This algorithm coded in accordance to Zagha et al paper 1991

struct EdgeList* radixSortEdgesBySource (struct EdgeList* edgeList){

	    // printf("*** START Radix Sort Edges By Source *** \n");

    // struct Graph* graph = graphNew(edgeList->num_vertices, edgeList->num_edges, inverse);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number

    uint32_t radix = 4;  // 32/8 8 bit radix needs 4 iterations
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t buckets = 256; // 2^radix = 256 buckets
    uint32_t num_edges = edgeList->num_edges;
    uint32_t* buckets_count = NULL;

    // omp_set_num_threads(P);
   
    uint32_t j = 0; //1,2,3 iteration

    struct EdgeList* sorted_edges_array = newEdgeList(num_edges);

    sorted_edges_array->num_vertices = edgeList->num_vertices;

    buckets_count = (uint32_t*) my_malloc(P * buckets * sizeof(uint32_t));
    

    for(j=0 ; j < radix ; j++){
        radixSortCountSortEdgesBySource (&sorted_edges_array, &edgeList, j, buckets, buckets_count);
    }

    free(buckets_count);
    freeEdgeList(sorted_edges_array);

    return edgeList;

}

// This algorithm coded in accordance to Zagha et al paper 1991

struct EdgeList* radixSortEdgesBySourceAndDestination (struct EdgeList* edgeList){

        // printf("*** START Radix Sort Edges By Source *** \n");

    // struct Graph* graph = graphNew(edgeList->num_vertices, edgeList->num_edges, inverse);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number

    uint32_t radix = 4;  // 32/8 8 bit radix needs 4 iterations
    uint32_t P = numThreads;  // 32/8 8 bit radix needs 4 iterations
    uint32_t buckets = 256; // 2^radix = 256 buckets
    uint32_t num_edges = edgeList->num_edges;
    uint32_t* buckets_count = NULL;

    // omp_set_num_threads(P);
   
    uint32_t j = 0; //1,2,3 iteration

    struct EdgeList* sorted_edges_array = newEdgeList(num_edges);

    sorted_edges_array->num_vertices = edgeList->num_vertices;
   
    buckets_count = (uint32_t*) my_malloc(P * buckets * sizeof(uint32_t));
   

    for(j=0 ; j < radix ; j++){
        radixSortCountSortEdgesByDestination (&sorted_edges_array, &edgeList, j, buckets, buckets_count);
    }
    for(j=0 ; j < radix ; j++){
        radixSortCountSortEdgesBySource (&sorted_edges_array, &edgeList, j, buckets, buckets_count);
    }
    

    free(buckets_count);
    freeEdgeList(sorted_edges_array);

    return edgeList;

}


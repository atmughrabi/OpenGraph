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
#include "mt19937.h"
#include "vc_vector.h"

#include "graphCSR.h"
#include "reorder.h"


uint32_t RegionAtomicDecrement(uint32_t *region)
{

    uint32_t oldValue;
    uint32_t flag = 0;

    do
    {

        oldValue = *region;
        if(oldValue > 0)
        {
            if(__sync_bool_compare_and_swap(region, oldValue, (oldValue - 1)))
            {
                flag = 1;
            }
        }
        else
        {
            return 0;
        }
    }
    while(!flag);

    return 1;
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
    uint32_t P = 1;  // 32/8 8 bit radix needs 4 iterations
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t base = 0;

    #pragma omp parallel default(none) shared(P,pageRanksFP, pageRanksFPTemp,radix,labels,labelsTemp,buckets,buckets_count, num_vertices) firstprivate(t_id, offset_end,offset_start,base,i,j,t,u,o)
    {

        t_id = omp_get_thread_num();

        if(t_id == 0)
        {
            P = omp_get_num_threads();
        }

        #pragma omp barrier

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

    buckets_count   = (uint32_t *) my_malloc(P * buckets * sizeof(uint32_t));
    pageRanksFP     = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    pageRanksFPTemp = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));
    labelsTemp      = (uint32_t *) my_malloc(num_vertices * sizeof(uint32_t));

    #pragma omp parallel for
    for(v = 0; v < num_vertices; v++)
    {
        pageRanksFP[v] = FloatToFixed32SORT(pageRanks[v]);
        pageRanksFPTemp[v] = 0;
        labelsTemp[v] = 0;
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

    #pragma omp parallel
    for (j = 0; j < num_vertices; ++j)
    {
        labelsTemp[j] = 0;
        degreesTemp[j] = 0;
    }


    for(j = 0 ; j < radix ; j++)
    {
        radixSortCountSortEdgesByRanks (&degrees, &degreesTemp, &labels, &labelsTemp, j, buckets, buckets_count, num_vertices);
    }


    free(buckets_count);
    free(degreesTemp);
    free(labelsTemp);

    return labels;

}

// ********************************************************************************************
// ***************                  Degree relabel                               **************
// ********************************************************************************************

struct EdgeList *reorderGraphProcessDegree( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode)
{
    uint32_t i;
    uint32_t *degrees;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));

    #pragma omp parallel
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        degrees[i] = 0;
    }

    degrees = reorderGraphGenerateInOutDegrees( degrees, edgeList, lmode);

    edgeList = reorderGraphListDegree( edgeList, degrees, lmode);

    free(degrees);
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


    switch(lmode)
    {
    case 1  :
        printf("| %-51s | \n", "OUT-DEGREE");
        break;
    case 2  :
        printf("| %-51s | \n", "IN-DEGREE");
        break;
    case 3  :
        printf("| %-51s | \n", "(IN+OUT)-DEGREE");
        break;
    case 10  :
        printf("| %-51s | \n", "RANDOM-DEGREE");
        break;
    default :
        printf("| %-51s | \n", "OUT-DEGREE");
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

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Reordering Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    Start(timer);
    edgeList = relabelEdgeList(edgeList, labels);
    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    #pragma omp parallel for
    for (v = 0; v < edgeList->num_vertices; ++v)
    {
        edgeList->label_array[v] = labels[edgeList->label_array[v]];
        edgeList->inverse_label_array[edgeList->label_array[v]] = v;
    }

    free(timer);
    free(labelsInverse);
    free(labels);

    return edgeList;
}


// ********************************************************************************************
// ***************                  DBG relabel                                  **************
// ********************************************************************************************

struct EdgeList *reorderGraphProcessDBG( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode)
{

    // UINT32_MAX
    uint32_t  i;
    uint32_t *degrees;
    uint32_t *thresholds;
    uint32_t  num_buckets = 11;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    thresholds = (uint32_t *) my_malloc(num_buckets * sizeof(uint32_t));

    #pragma omp parallel
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        degrees[i] = 0;
    }

    // START initialize thresholds
    if(edgeList->avg_degree <= 1)
        thresholds[0] = 1;
    else
        thresholds[0] = (edgeList->avg_degree / 2);
    for ( i = 1; i < (num_buckets - 1); ++i)
    {
        thresholds[i] = thresholds[i - 1] * 2;
    }
    thresholds[num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds

    switch(lmode)
    {
    case 4  :
        printf("| %-51s | \n", "DBG OUT-DEGREE");
        break;
    case 5  :
        printf("| %-51s | \n", "DBG IN-DEGREE");
        break;
    default :
        printf("| %-51s | \n", "DBG OUT-DEGREE");
    }

    degrees = reorderGraphGenerateInOutDegrees(degrees, edgeList, lmode);

    edgeList = reorderGraphListDBG(edgeList, degrees, thresholds, num_buckets, lmode);

    free(thresholds);
    free(degrees);
    return edgeList;

}

struct EdgeList *reorderGraphListDBG(struct EdgeList *edgeList, uint32_t *degrees, uint32_t *thresholds, uint32_t num_buckets, uint32_t lmode)
{

    uint32_t  i = 0;
    int32_t  j = 0;
    int32_t  k = 0;
    void  *iter = 0;
    uint32_t  v = 0;
    uint32_t  t = 0;
    uint32_t  temp_idx = 0;
    uint32_t P = 1;
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;

    uint32_t *start_idx = NULL;
    vc_vector **buckets = NULL;
    uint32_t *labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));


    Start(timer);
    #pragma omp parallel default(none) shared(P,labels,buckets,edgeList,num_buckets,degrees,thresholds,start_idx) firstprivate(iter,temp_idx,k,offset_start,offset_end,t_id,i,j,v,t)
    {

        t_id = omp_get_thread_num();

        if(t_id == 0)
        {
            P = omp_get_num_threads();
            start_idx = (uint32_t *) my_malloc(P * num_buckets * sizeof(uint32_t));
            buckets = (vc_vector **) malloc(P * num_buckets * sizeof(vc_vector *));
        }

        #pragma omp barrier

        for (i = 0; i < num_buckets; ++i)
        {
            buckets[(t_id * num_buckets) + i] = vc_vector_create(0, sizeof(uint32_t), NULL);
        }

        offset_start = t_id * (edgeList->num_vertices / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (edgeList->num_vertices / P) + (edgeList->num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (edgeList->num_vertices / P);
        }

        for (v = offset_start; v < offset_end; ++v)
        {
            for ( i = 0; i < num_buckets; ++i)
            {
                if(degrees[v] <= thresholds[i])
                {
                    vc_vector_push_back(buckets[(t_id * num_buckets) + i], &v);
                    break;
                }
            }
        }

        #pragma omp barrier

        if(t_id == 0)
        {
            for ( j = num_buckets - 1; j >= 0; --j)
            {
                for (t = 0; t < P; ++t)
                {
                    start_idx[(t * num_buckets) + j] = temp_idx;
                    temp_idx += vc_vector_count(buckets[(t * num_buckets) + j]);
                }
            }
        }

        #pragma omp barrier

        for ( j = num_buckets - 1 ; j >= 0 ; --j)
        {
            k = start_idx[(t_id * num_buckets) + j];
            for (   iter = vc_vector_begin(buckets[(t_id * num_buckets) + j]);
                    iter != vc_vector_end(buckets[(t_id * num_buckets) + j]);
                    iter = vc_vector_next(buckets[(t_id * num_buckets) + j], iter))
            {
                labels[(*(uint32_t *)iter)] = k++;
            }
        }

    }

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Reordering Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    Start(timer);
    edgeList = relabelEdgeList(edgeList, labels);
    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    #pragma omp parallel for
    for (v = 0; v < edgeList->num_vertices; ++v)
    {
        edgeList->label_array[v] = labels[edgeList->label_array[v]];
        edgeList->inverse_label_array[edgeList->label_array[v]] = v;
    }

    for (i = 0; i < (P * num_buckets); ++i)
    {
        vc_vector_release(buckets[i]);
    }

    free(timer);
    free(buckets);
    free(start_idx);
    free(labels);
    return edgeList;
}

// ********************************************************************************************
// ***************                  HUBSort relabel                              **************
// ********************************************************************************************

struct EdgeList *reorderGraphProcessHUBSort( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode)
{

    // UINT32_MAX
    uint32_t  i;
    uint32_t *degrees;
    uint32_t *thresholds;
    uint32_t  num_buckets = 2;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    thresholds = (uint32_t *) my_malloc(num_buckets * sizeof(uint32_t));

    #pragma omp parallel for
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        degrees[i] = 0;
    }

    // START initialize thresholds
    if(edgeList->avg_degree <= 1)
        thresholds[0] = 1;
    else
        thresholds[0] = (edgeList->avg_degree / 2);
    for ( i = 1; i < (num_buckets - 1); ++i)
    {
        thresholds[i] = thresholds[i - 1] * 2;
    }
    thresholds[num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds

    switch(lmode)
    {
    case 6  :
        printf("| %-51s | \n", "HUBSort OUT-DEGREE");
        break;
    case 7  :
        printf("| %-51s | \n", "HUBSort IN-DEGREE");
        break;
    default :
        printf("| %-51s | \n", "HUBSort OUT-DEGREE");
    }

    degrees = reorderGraphGenerateInOutDegrees(degrees, edgeList, lmode);

    edgeList = reorderGraphListHUBSort(edgeList, degrees, thresholds, num_buckets, lmode);

    free(thresholds);
    free(degrees);
    return edgeList;

}

struct EdgeList *reorderGraphListHUBSort(struct EdgeList *edgeList, uint32_t *degrees, uint32_t *thresholds, uint32_t num_buckets, uint32_t lmode)
{

    uint32_t  i = 0;
    int32_t  j = 0;
    int32_t  k = 0;
    void  *iter = 0;
    uint32_t  v = 0;
    uint32_t  t = 0;
    uint32_t  temp_idx = 0;
    uint32_t P = 1;
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;

    uint32_t *start_idx = NULL;
    vc_vector **buckets = NULL;
    uint32_t *labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    uint32_t *sizeHot = (uint32_t *) my_malloc(num_buckets * sizeof(uint32_t));
    uint32_t **degreesHot  = (uint32_t **) my_malloc(num_buckets * sizeof(uint32_t *));
    uint32_t **verticesHot  = (uint32_t **) my_malloc(num_buckets * sizeof(uint32_t *));

    Start(timer);
    #pragma omp parallel default(none) shared(P,verticesHot,degreesHot,sizeHot,labels,buckets,edgeList,num_buckets,degrees,thresholds,start_idx) firstprivate(iter,temp_idx,k,offset_start,offset_end,t_id,i,j,v,t)
    {
        t_id = omp_get_thread_num();

        if(t_id == 0)
        {
            P = omp_get_num_threads();
            start_idx = (uint32_t *) my_malloc(P * num_buckets * sizeof(uint32_t));
            buckets = (vc_vector **) malloc(P * num_buckets * sizeof(vc_vector *));
        }

        #pragma omp barrier

        for (i = 0; i < num_buckets; ++i)
        {
            buckets[(t_id * num_buckets) + i] = vc_vector_create(0, sizeof(uint32_t), NULL);
        }

        offset_start = t_id * (edgeList->num_vertices / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (edgeList->num_vertices / P) + (edgeList->num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (edgeList->num_vertices / P);
        }

        for (v = offset_start; v < offset_end; ++v)
        {
            for ( i = 0; i < num_buckets; ++i)
            {
                if(degrees[v] <= thresholds[i])
                {
                    vc_vector_push_back(buckets[(t_id * num_buckets) + i], &v);
                    break;
                }
            }
        }

        #pragma omp barrier

        if(t_id == 0)
        {
            for ( j = num_buckets - 1; j >= 0; --j)
            {
                temp_idx = 0;
                for (t = 0; t < P; ++t)
                {
                    start_idx[(t * num_buckets) + j] = temp_idx;
                    temp_idx += vc_vector_count(buckets[(t * num_buckets) + j]);
                }

                sizeHot[j] = temp_idx;
                degreesHot[j]  = (uint32_t *) my_malloc(sizeHot[j] * sizeof(uint32_t));
                verticesHot[j]  = (uint32_t *) my_malloc(sizeHot[j] * sizeof(uint32_t));
            }
        }

        #pragma omp barrier

        for ( j = num_buckets - 1 ; j >= 0 ; --j)
        {
            k = start_idx[(t_id * num_buckets) + j];
            for (   iter = vc_vector_begin(buckets[(t_id * num_buckets) + j]);
                    iter != vc_vector_end(buckets[(t_id * num_buckets) + j]);
                    iter = vc_vector_next(buckets[(t_id * num_buckets) + j], iter))
            {

                verticesHot[j][k] = (*(uint32_t *)iter);
                degreesHot[j][k] = degrees[(*(uint32_t *)iter)];
                k++;
            }
        }

    }

    verticesHot[num_buckets - 1] = radixSortEdgesByDegree(degreesHot[num_buckets - 1], verticesHot[num_buckets - 1], sizeHot[num_buckets - 1]);

    #pragma omp parallel for
    for(v = 0; v < sizeHot[1]; v++)
    {
        labels[verticesHot[1][v]] = sizeHot[1] - 1 - v;
    }

    #pragma omp parallel for
    for(v = 0; v < sizeHot[0]; v++)
    {
        labels[verticesHot[0][v]] = sizeHot[1] + (v);
    }

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Reordering Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    Start(timer);
    edgeList = relabelEdgeList(edgeList, labels);
    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    #pragma omp parallel for
    for (v = 0; v < edgeList->num_vertices; ++v)
    {
        edgeList->label_array[v] = labels[edgeList->label_array[v]];
        edgeList->inverse_label_array[edgeList->label_array[v]] = v;
    }

    for (i = 0; i < (P * num_buckets); ++i)
    {
        vc_vector_release(buckets[i]);
    }

    for (i = 0; i <  num_buckets; ++i)
    {
        free(degreesHot[i]);
        free(verticesHot[i]);
    }

    free(degreesHot);
    free(verticesHot);
    free(sizeHot);

    free(timer);
    free(buckets);
    free(start_idx);
    free(labels);
    return edgeList;
}



// ********************************************************************************************
// ***************                  HUBCluster relabel                           **************
// ********************************************************************************************

struct EdgeList *reorderGraphProcessHUBCluster( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode)
{


    // UINT32_MAX
    uint32_t  i;
    uint32_t *degrees;
    uint32_t *thresholds;
    uint32_t  num_buckets = 2;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    thresholds = (uint32_t *) my_malloc(num_buckets * sizeof(uint32_t));

    #pragma omp parallel for
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        degrees[i] = 0;
    }

    // START initialize thresholds
    thresholds[0] = (edgeList->avg_degree);
    thresholds[num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds

    switch(lmode)
    {
    case 8  :
        printf("| %-51s | \n", "HUBCluster OUT-DEGREE");
        break;
    case 9  :
        printf("| %-51s | \n", "HUBCluster IN-DEGREE");
        break;
    default :
        printf("| %-51s | \n", "HUBCluster OUT-DEGREE");
    }

    degrees = reorderGraphGenerateInOutDegrees(degrees, edgeList, lmode);

    edgeList = reorderGraphListHUBCluster(edgeList, degrees, thresholds, num_buckets, lmode);

    free(thresholds);
    free(degrees);
    return edgeList;

}

struct EdgeList *reorderGraphListHUBCluster(struct EdgeList *edgeList, uint32_t *degrees, uint32_t *thresholds, uint32_t num_buckets, uint32_t lmode)
{

    uint32_t i = 0;
    int32_t  j = 0;
    int32_t  k = 0;
    void  *iter = 0;
    uint32_t v = 0;
    uint32_t t = 0;
    uint32_t temp_idx = 0;
    uint32_t P = 1;
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;

    uint32_t *start_idx = NULL;
    vc_vector **buckets = NULL;

    uint32_t *labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    Start(timer);
    #pragma omp parallel default(none) shared(P,labels,buckets,edgeList,num_buckets,degrees,thresholds,start_idx) firstprivate(iter,temp_idx,k,offset_start,offset_end,t_id,i,j,v,t)
    {
        t_id = omp_get_thread_num();

        if(t_id == 0)
        {
            P = omp_get_num_threads();
            start_idx = (uint32_t *) my_malloc(P * num_buckets * sizeof(uint32_t));
            buckets = (vc_vector **) malloc(P * num_buckets * sizeof(vc_vector *));
        }

        #pragma omp barrier

        for (i = 0; i < num_buckets; ++i)
        {
            buckets[(t_id * num_buckets) + i] = vc_vector_create(0, sizeof(uint32_t), NULL);
        }

        offset_start = t_id * (edgeList->num_vertices / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (edgeList->num_vertices / P) + (edgeList->num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (edgeList->num_vertices / P);
        }

        for (v = offset_start; v < offset_end; ++v)
        {
            for ( i = 0; i < num_buckets; ++i)
            {
                if(degrees[v] <= thresholds[i])
                {
                    vc_vector_push_back(buckets[(t_id * num_buckets) + i], &v);
                    break;
                }
            }
        }

        #pragma omp barrier

        if(t_id == 0)
        {
            for ( j = num_buckets - 1; j >= 0; --j)
            {
                for (t = 0; t < P; ++t)
                {
                    start_idx[(t * num_buckets) + j] = temp_idx;
                    temp_idx += vc_vector_count(buckets[(t * num_buckets) + j]);
                }
            }
        }

        #pragma omp barrier

        for ( j = num_buckets - 1 ; j >= 0 ; --j)
        {
            k = start_idx[(t_id * num_buckets) + j];
            for (   iter = vc_vector_begin(buckets[(t_id * num_buckets) + j]);
                    iter != vc_vector_end(buckets[(t_id * num_buckets) + j]);
                    iter = vc_vector_next(buckets[(t_id * num_buckets) + j], iter))
            {
                labels[(*(uint32_t *)iter)] = k++;
            }
        }

    }

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Reordering Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    Start(timer);
    edgeList = relabelEdgeList(edgeList, labels);
    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Relabeling Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    #pragma omp parallel for
    for (v = 0; v < edgeList->num_vertices; ++v)
    {
        edgeList->label_array[v] = labels[edgeList->label_array[v]];
        edgeList->inverse_label_array[edgeList->label_array[v]] = v;
    }

    for (i = 0; i < (P * num_buckets); ++i)
    {
        vc_vector_release(buckets[i]);
    }

    free(timer);
    free(buckets);
    free(start_idx);
    free(labels);
    return edgeList;
}

// ********************************************************************************************
// ***************                  AccelGraph label-Masking                     **************
// ********************************************************************************************

struct EdgeList *maskGraphProcess(struct EdgeList *edgeList, struct Arguments *arguments)
{

    struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));

    printf(" *****************************************************\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Mask Process");
    printf(" -----------------------------------------------------\n");
    Start(timer);

    uint32_t cache_size = (32768 >> 2);

    switch(arguments->lmode)
    {
    case 1  :
    case 2  :
    case 3  :
    case 4 :
        edgeList = maskGraphProcessDegree( edgeList, arguments->mmode, cache_size); // degree
        break;
    default :
        edgeList = maskGraphProcessDegree( edgeList, arguments->mmode, cache_size); // out-degree
    }

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Total Mask Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");
    printf(" *****************************************************\n");

    free(timer);

    return edgeList;
}

struct EdgeList *maskGraphProcessDegree( struct EdgeList *edgeList, uint32_t mmode, uint32_t cache_size)
{

    // UINT32_MAX
    uint32_t  i;
    uint32_t *degrees;
    uint32_t *thresholds;
    uint32_t  num_buckets = 11;

    degrees = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    thresholds = (uint32_t *) my_malloc(num_buckets * sizeof(uint32_t));
    #pragma omp parallel
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        degrees[i] = 0;
    }

    // START initialize thresholds
    if(edgeList->avg_degree <= 1)
        thresholds[0] = 1;
    else
        thresholds[0] = (edgeList->avg_degree / 2);
    for ( i = 1; i < (num_buckets - 1); ++i)
    {
        thresholds[i] = thresholds[i - 1] * 2;
    }
    thresholds[num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds

    switch(mmode)
    {
    case 1  :
        printf("| %-51s | \n", "Vertex Property OUT-DEGREE");
        break;
    case 2  :
        printf("| %-51s | \n", "Vertex Structure IN-DEGREE");
        break;
    case 3  :
        printf("| %-51s | \n", "Vertex Property OUT-DEGREE");
        break;
    case 4  :
        printf("| %-51s | \n", "Vertex Structure IN-DEGREE");
        break;
    default :
        printf("| %-51s | \n", "Vertex Property OUT-DEGREE");
    }

    degrees = maskGraphProcessGenerateInOutDegrees(degrees, edgeList, mmode);

    edgeList = maskGraphProcessGenerateMaskArray(edgeList, degrees, thresholds, num_buckets, mmode, cache_size);

    free(thresholds);
    free(degrees);
    return edgeList;

}

uint32_t *maskGraphProcessGenerateInOutDegrees(uint32_t *degrees, struct EdgeList *edgeList, uint32_t mmode)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;

    #pragma omp parallel for default(none) private(i,src,dest) shared(edgeList,degrees,mmode)
    for(i = 0; i < edgeList->num_edges; i++)
    {
        src  = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        switch(mmode)
        {
        case 1 :
        case 3 :
        {
            #pragma omp atomic update
            degrees[src]++;
        }
        break;
        case 2 :
        case 4 :
        {
            #pragma omp atomic update
            degrees[dest]++;
        }
        break;
        case 5 :
        case 6 :
        {
            #pragma omp atomic update
            degrees[dest]++;
            #pragma omp atomic update
            degrees[src]++;
        }
        break;
        default :
        {
            #pragma omp atomic update
            degrees[src]++;
        }// out-degree
        }
    }
    return degrees;
}

struct EdgeList *maskGraphProcessGenerateMaskArray(struct EdgeList *edgeList, uint32_t *degrees, uint32_t *thresholds, uint32_t num_buckets, uint32_t mmode, uint32_t cache_size)
{

    uint32_t i = 0;
    int32_t  j = 0;
    void  *iter = 0;
    uint32_t  v = 0;
    uint32_t  t = 0;
    uint32_t  temp_idx = 0;
    uint32_t P = 1;
    uint32_t t_id = 0;
    uint32_t offset_start = 0;
    uint32_t offset_end = 0;
    uint32_t  num_masks = 4;

    uint32_t *start_idx      = NULL;
    vc_vector **buckets      = NULL;

    uint32_t *labels         = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    uint32_t *mask_array     = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));
    struct Timer *timer      = (struct Timer *) malloc(sizeof(struct Timer));
    uint32_t *cache_regions  = (uint32_t *) my_malloc(num_masks * sizeof(uint32_t));

    int diff = (int)edgeList->num_vertices - (int)cache_size;

    if(diff < (2 * (int)cache_size))
    {
        cache_regions[0] = edgeList->num_vertices / 3; // VERTEX_VALUE_HOT_U32
        cache_regions[1] = edgeList->num_vertices / 3; // VERTEX_CACHE_WARM_U32
        cache_regions[2] = edgeList->num_vertices / 3; // VERTEX_VALUE_LUKEWARM_U32
    }
    else
    {
        cache_regions[0] = cache_size*8; // VERTEX_VALUE_HOT_U32
        cache_regions[1] = cache_regions[0]*4; // VERTEX_CACHE_WARM_U32
        cache_regions[2] = cache_regions[1]*4; // VERTEX_VALUE_LUKEWARM_U32
    }

    cache_regions[3] = UINT32_MAX; // VERTEX_CACHE_COLD_U32

    #pragma omp parallel for
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        mask_array[i] = VERTEX_CACHE_COLD_U32;
    }

    Start(timer);
    #pragma omp parallel default(none) shared(P,mask_array,mmode,cache_regions,labels,buckets,edgeList,num_buckets,degrees,thresholds,start_idx) firstprivate(iter,temp_idx,offset_start,offset_end,t_id,i,j,v,t)
    {
        t_id = omp_get_thread_num();

        if(t_id == 0)
        {
            P = omp_get_num_threads();
            start_idx = (uint32_t *) my_malloc(P * num_buckets * sizeof(uint32_t));
            buckets = (vc_vector **) malloc(P * num_buckets * sizeof(vc_vector *));
        }

        #pragma omp barrier

        for (i = 0; i < num_buckets; ++i)
        {
            buckets[(t_id * num_buckets) + i] = vc_vector_create(0, sizeof(uint32_t), NULL);
        }

        offset_start = t_id * (edgeList->num_vertices / P);

        if(t_id == (P - 1))
        {
            offset_end = offset_start + (edgeList->num_vertices / P) + (edgeList->num_vertices % P) ;
        }
        else
        {
            offset_end = offset_start + (edgeList->num_vertices / P);
        }

        for (v = offset_start; v < offset_end; ++v)
        {
            for ( i = 0; i < num_buckets; ++i)
            {
                if(degrees[v] <= thresholds[i])
                {
                    vc_vector_push_back(buckets[(t_id * num_buckets) + i], &v);
                    break;
                }
            }
        }

        #pragma omp barrier

        if(t_id == 0)
        {
            for ( j = num_buckets - 1; j >= 0; --j)
            {
                for (t = 0; t < P; ++t)
                {
                    start_idx[(t * num_buckets) + j] = temp_idx;
                    temp_idx += vc_vector_count(buckets[(t * num_buckets) + j]);
                }
            }
        }

        #pragma omp barrier

        for ( j = num_buckets - 1 ; j >= 0 ; --j)
        {
            for (   iter = vc_vector_begin(buckets[(t_id * num_buckets) + j]);
                    iter != vc_vector_end(buckets[(t_id * num_buckets) + j]);
                    iter = vc_vector_next(buckets[(t_id * num_buckets) + j], iter))
            {
                if(RegionAtomicDecrement(&(cache_regions[0])))
                {
                    mask_array[(*(uint32_t *)iter)] = VERTEX_VALUE_HOT_U32;
                }
                else if(RegionAtomicDecrement(&(cache_regions[1])))
                {
                    mask_array[(*(uint32_t *)iter)] = VERTEX_CACHE_WARM_U32;
                }
                else if(RegionAtomicDecrement(&(cache_regions[2])))
                {
                    mask_array[(*(uint32_t *)iter)] = VERTEX_VALUE_LUKEWARM_U32;
                }
                else
                {
                    mask_array[(*(uint32_t *)iter)] = VERTEX_CACHE_COLD_U32;
                }
            }
        }

    }

    Stop(timer);

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Mask Complete");
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", Seconds(timer));
    printf(" -----------------------------------------------------\n");

    if(mmode == 1 || mmode == 2)
    {
        Start(timer);
        edgeList = maskEdgeList(edgeList, mask_array);
        Stop(timer);
        printf(" -----------------------------------------------------\n");
        printf("| %-51s | \n", "Relabeling Complete");
        printf(" -----------------------------------------------------\n");
        printf("| %-51f | \n", Seconds(timer));
        printf(" -----------------------------------------------------\n");
    }

    #pragma omp parallel for
    for (i = 0; i < edgeList->num_vertices; ++i)
    {
        edgeList->mask_array[i] = mask_array[i];
    }

    for (i = 0; i < (P * num_buckets); ++i)
    {
        vc_vector_release(buckets[i]);
    }

    free(mask_array);
    free(timer);
    free(buckets);
    free(start_idx);
    free(labels);
    free(cache_regions);
    return edgeList;
}

// ********************************************************************************************
// ***************                  generic functions                            **************
// ********************************************************************************************

uint32_t *reorderGraphGenerateInOutDegrees(uint32_t *degrees, struct EdgeList *edgeList, uint32_t lmode)
{

    uint32_t i;
    uint32_t src;
    uint32_t dest;



    if(lmode != 10)
    {
        #pragma omp parallel for default(none) private(i,src,dest) shared(edgeList,degrees,lmode)
        for(i = 0; i < edgeList->num_edges; i++)
        {
            src  = edgeList->edges_array_src[i];
            dest = edgeList->edges_array_dest[i];

            switch(lmode)
            {
            case 1 :
            case 4 :
            case 6 :
            case 8 :
            {
                #pragma omp atomic update
                degrees[src]++;
            } // degree
            break;
            case 2 :
            case 5 :
            case 7 :
            case 9 :
            {
                #pragma omp atomic update
                degrees[dest]++;
            }
            break;
            case 3 :
            {
                #pragma omp atomic update
                degrees[dest]++;
                #pragma omp atomic update
                degrees[src]++;
            }
            break;
            default :
            {
                #pragma omp atomic update
                degrees[src]++;
            }// out-degree
            }
        }
    }


    if(lmode == 10)
    {
        mt19937state *mt19937var = (mt19937state *) my_malloc(sizeof(mt19937state));
        initializeMersenneState (mt19937var, 27491095);
        #pragma omp parallel for firstprivate(mt19937var)
        for (i = 0; i < edgeList->num_vertices; ++i)
        {
            degrees[i] = (generateRandInt(mt19937var) % edgeList->num_vertices) + omp_get_thread_num();
        }
        free(mt19937var);
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
    case 1  :
    case 2  :
    case 3  :
    case 10 :
        edgeList = reorderGraphProcessDegree( arguments->sort, edgeList, arguments->lmode);// degree
        break;
    case 4  :
    case 5  :
        edgeList = reorderGraphProcessDBG( arguments->sort, edgeList, arguments->lmode);// DBG
        break;
    case 6  :
    case 7  :
        edgeList = reorderGraphProcessHUBSort( arguments->sort, edgeList, arguments->lmode);// HUBSort
        break;
    case 8  :
    case 9  :
        edgeList = reorderGraphProcessHUBCluster( arguments->sort, edgeList, arguments->lmode);// HUBCluster
        break;
    case 11 :
        edgeList = relabelEdgeListFromFile(edgeList, arguments->fnamel, edgeList->num_vertices);// load from file
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

struct EdgeList *maskEdgeList(struct EdgeList *edgeList, uint32_t *mask_array)
{
    uint32_t i;

    #pragma omp parallel for
    for(i = 0; i < edgeList->num_edges; i++)
    {
        uint32_t src;
        uint32_t dest;
        src = edgeList->edges_array_src[i];
        dest = edgeList->edges_array_dest[i];

        edgeList->edges_array_src[i] = src | mask_array[src];
        edgeList->edges_array_dest[i] = dest | mask_array[dest];
    }

    return edgeList;
}

// ********************************************************************************************
// ***************                  File relabel                                 **************
// ********************************************************************************************

struct EdgeList *relabelEdgeListFromFile(struct EdgeList *edgeList, const char *fnameb, uint32_t size)
{

    FILE *pText;
    uint32_t i;
    uint32_t v = 0;
    uint32_t dest = 0;
    uint32_t x = 0;
    uint32_t *labels = (uint32_t *) my_malloc(edgeList->num_vertices * sizeof(uint32_t));

    // char *fname_txt = (char *) malloc((strlen(fnameb) + 10) * sizeof(char));

    // fname_txt = strcpy (fname_txt, fnameb);

    printf("%s\n", fnameb );
    pText = fopen(fnameb, "r");

    if (pText == NULL)
    {
        return NULL;
    }

    while (1)
    {

        i = fscanf(pText, "%u\n", &dest);
        labels[x] = dest;
        x++;

        if( x == edgeList->num_vertices )
            break;

        if( i == EOF )
            break;

    }
    fclose(pText);



    edgeList = relabelEdgeList(edgeList, labels);

    #pragma omp parallel for
    for (v = 0; v < edgeList->num_vertices; ++v)
    {
        edgeList->label_array[v] = labels[edgeList->label_array[v]];
        edgeList->inverse_label_array[edgeList->label_array[v]] = v;
    }

    free(labels);
    // free(fname_txt);

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

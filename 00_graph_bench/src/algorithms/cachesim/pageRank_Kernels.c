
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "myMalloc.h"
#include "fixedPoint.h"
#include "quantization.h"

#include "graphGrid.h"

#include "cache.h"
#include "pageRank_Kernels.h"

float movingAvg(float *ptrArrNumbers, double *ptrSum, int pos, int len, float nextNum)
{
    //Subtract the oldest number from the prev sum, add the new number
    *ptrSum = *ptrSum - ptrArrNumbers[pos] + nextNum;
    //Assign the nextNum to the position in the array
    ptrArrNumbers[pos] = nextNum;
    //return the average
    return *ptrSum / len;
}

// ********************************************************************************************
// ***************          GRID DataStructure               **************
// ********************************************************************************************

void pageRankPullRowGraphGridKernel(float *riDividedOnDiClause_pull_grid, float *pageRanksNext_pull_grid,  struct Partition *partitions, uint32_t totalPartitions)
{

    uint32_t i;
    for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
    {
        uint32_t j;

        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t k;
            uint32_t src;
            uint32_t dest;

            for (k = 0; k < partitions[(i * totalPartitions) + j].num_edges; ++k)
            {
                src  = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_src[k];
                dest = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_dest[k];

                pageRanksNext_pull_grid[dest] +=  riDividedOnDiClause_pull_grid[src];
            }
        }
    }
}

// ********************************************************************************************

void pageRankPullRowFixedPointGraphGridKernel(uint64_t *riDividedOnDiClause_pull_grid_fp, uint64_t *pageRanksNext_pull_grid_fp,  struct Partition *partitions, uint32_t totalPartitions)
{

    uint32_t i;
    for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
    {
        uint32_t j;

        for (j = 0; j < totalPartitions; ++j)
        {
            uint32_t k;
            uint32_t src;
            uint32_t dest;

            for (k = 0; k < partitions[(i * totalPartitions) + j].num_edges; ++k)
            {
                src  = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_src[k];
                dest = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_dest[k];

                pageRanksNext_pull_grid_fp[dest] +=  riDividedOnDiClause_pull_grid_fp[src];
            }
        }
    }
}

// ********************************************************************************************

void pageRankPushColumnGraphGridKernel(float *riDividedOnDiClause_push_grid, float *pageRanksNext_push_grid,  struct Partition *partitions, uint32_t totalPartitions)
{

    uint32_t j;

    for (j = 0; j < totalPartitions; ++j)
    {
        uint32_t i;
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t k;
            uint32_t src;
            uint32_t dest;

            for (k = 0; k < partitions[(i * totalPartitions) + j].num_edges; ++k)
            {
                src  = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_src[k];
                dest = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_dest[k];

                pageRanksNext_push_grid[dest] +=  riDividedOnDiClause_push_grid[src];
            }
        }
    }
}

// ********************************************************************************************

void pageRankPushColumnFixedPointGraphGridKernel(uint64_t *riDividedOnDiClause_push_grid_fp, uint64_t *pageRanksNext_push_grid_fp,  struct Partition *partitions, uint32_t totalPartitions)
{

    uint32_t j;

    for (j = 0; j < totalPartitions; ++j)
    {
        uint32_t i;
        for (i = 0; i < totalPartitions; ++i)  // iterate over partitions rowwise
        {
            uint32_t k;
            uint32_t src;
            uint32_t dest;

            for (k = 0; k < partitions[(i * totalPartitions) + j].num_edges; ++k)
            {
                src  = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_src[k];
                dest = (partitions[(i * totalPartitions) + j]).edgeList->edges_array_dest[k];

                pageRanksNext_push_grid_fp[dest] +=  riDividedOnDiClause_push_grid_fp[src];
            }
        }
    }
}

// ********************************************************************************************
// ***************          CSR DataStructure                                    **************
// ********************************************************************************************


void pageRankPullGraphCSRKernel(float *riDividedOnDiClause_pull_csr, float *pageRanksNext_pull_csr, uint32_t *out_degree_pull_csr, uint32_t *edges_idx_pull_csr, uint32_t *sorted_edges_array_pull_csr, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;



    for(v = 0; v < num_vertices; v++)
    {
        float nodeIncomingPR = 0.0f;
        degree = out_degree_pull_csr[v];
        edge_idx = edges_idx_pull_csr[v];

        for(j = edge_idx ; j <  (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array_pull_csr[j];
            nodeIncomingPR += riDividedOnDiClause_pull_csr[u]; // pageRanks[v]/graph->vertices[v].out_degree;
        }
        pageRanksNext_pull_csr[v] = nodeIncomingPR;
    }
}

void pageRankPullGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause, float *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices)
{


    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;

    // #pragma omp parallel for private(v,j,u,degree,edge_idx) schedule(static, 8)
    for(v = 0; v < num_vertices; v++)
    {

        float nodeIncomingPR = 0.0f;
        degree = out_degree[v];
        edge_idx = edges_idx[v];

        // Access(cache->ref_cache, (uint64_t) & (out_degree[v]), 'r', v);
        // Access(cache->ref_cache, (uint64_t) & (edges_idx[v]), 'r', v);
        // AccessDoubleTaggedCacheFloat(cache, (uint64_t) & (out_degree[v]), 'r', v, out_degree[v]);
        // AccessDoubleTaggedCacheFloat(cache, (uint64_t) & (edges_idx[v]), 'r', v, edges_idx[v]);

        for(j = edge_idx ; j <  (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array[j];

            // Access(cache->accel_graph->cold_cache, (uint64_t) & (sorted_edges_array[j]), 'r', u);
            // Access(cache->ref_cache, (uint64_t) & (sorted_edges_array[j]), 'r', u);
            // AccessDoubleTaggedCacheFloat(cache, (uint64_t) & (sorted_edges_array[j]), 'r', u, sorted_edges_array[j]);

            nodeIncomingPR += riDividedOnDiClause[u]; // pageRanks[v]/graph->vertices[v].out_degree;
            // #pragma omp critical
            // {
            AccessDoubleTaggedCacheFloat(cache, (uint64_t) & (riDividedOnDiClause[u]), 'r', u, riDividedOnDiClause[u]);
            // }
        }
        pageRanksNext[v] = nodeIncomingPR;

        // #pragma omp critical
        // {
        AccessDoubleTaggedCacheFloat(cache, (uint64_t) & (pageRanksNext[v]), 'w', v, pageRanksNext[v]);
        // }
    }


}


// ********************************************************************************************

void pageRankPushGraphCSRKernel(float *riDividedOnDiClause_push_csr, float *pageRanksNext_push_csr, uint32_t *out_degree_push_csr, uint32_t *edges_idx_push_csr, uint32_t *sorted_edges_array_push_csr, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;



    for(v = 0; v < num_vertices; v++)
    {
        degree = out_degree_push_csr[v];
        edge_idx = edges_idx_push_csr[v];

        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array_push_csr[j];
            pageRanksNext_push_csr[u] += riDividedOnDiClause_push_csr[v];
        }
    }
}


void pageRankPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause, float *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if((v + 1) < num_vertices)
        {
            edge_idx = edges_idx[v + 1];
            for(j = edge_idx ; j < (edge_idx + out_degree[v + 1]) ; j++)
            {
                u = sorted_edges_array[j];
                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[u])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
                }

            }

            if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[v + 1])))
            {
                Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[v + 1]), 's', (v + 1));
            }
        }
#endif

        degree = out_degree[v];
        edge_idx = edges_idx[v];

        // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree[v]), 'r', v);
        // Access(cache->accel_graph->cold_cache, (uint64_t) & (edges_idx[v]), 'r', v);

        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array[j];

            Access(cache->accel_graph->cold_cache, (uint64_t) & (sorted_edges_array[j]), 'r', u);

            pageRanksNext[u] += riDividedOnDiClause[v];

            Access(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[v]), 'r', v);

            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'w', u);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
        }
    }
}

// ********************************************************************************************



void pageRankPullFixedPointGraphCSRKernel(uint64_t *riDividedOnDiClause_pull_csr_fp, uint64_t *pageRanksNext_pull_csr_fp, uint32_t *out_degree_pull_csr_fp, uint32_t *edges_idx_pull_csr_fp, uint32_t *sorted_edges_array_pull_csr_fp, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;



    for(v = 0; v < num_vertices; v++)
    {

        uint64_t nodeIncomingPR = 0;
        degree = out_degree_pull_csr_fp[v];
        edge_idx = edges_idx_pull_csr_fp[v];

        for(j = edge_idx ; j <  (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array_pull_csr_fp[j];
            nodeIncomingPR += riDividedOnDiClause_pull_csr_fp[u]; // pageRanks[v]/graph->vertices[v].out_degree;
        }
        pageRanksNext_pull_csr_fp[v] = nodeIncomingPR;
    }

}

void pageRankPullFixedPointGraphCSRKernelCache(struct DoubleTaggedCache *cache, uint64_t *riDividedOnDiClause, uint64_t *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices)
{


    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if((v + 1) < num_vertices)
        {
            edge_idx = edges_idx[v + 1];
            for(j = edge_idx ; j < (edge_idx + out_degree[v + 1]) ; j++)
            {
                u = sorted_edges_array[j];
                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[u])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[u]), 's', u);
                }
            }
            if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[v + 1])))
            {
                Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[v + 1]), 'r', (v + 1));
            }
        }
#endif

        uint64_t nodeIncomingPR = 0;
        degree = out_degree[v];
        edge_idx = edges_idx[v];

        // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree[v]), 'r', v);
        // Access(cache->accel_graph->cold_cache, (uint64_t) & (edges_idx[v]), 'r', v);

        for(j = edge_idx ; j <  (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array[j];

            Access(cache->accel_graph->cold_cache, (uint64_t) & (sorted_edges_array[j]), 'r', u);

            nodeIncomingPR += riDividedOnDiClause[u]; // pageRanks[v]/graph->vertices[v].out_degree;

            Access(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[u]), 'r', u);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[u]), 'r', u);

        }

        pageRanksNext[v] = nodeIncomingPR;
        Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[v]), 'r', v);
        Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[v]), 'w', v);
        Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[v]), 'r', v);
    }


}


// ********************************************************************************************

void pageRankPushFixedPointGraphCSRKernel(uint64_t *riDividedOnDiClause_push_csr_fp, uint64_t *pageRanksNext_push_csr_fp, uint32_t *out_degree_push_csr_fp, uint32_t *edges_idx_push_csr_fp, uint32_t *sorted_edges_array_push_csr_fp, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;



    for(v = 0; v < num_vertices; v++)
    {
        degree = out_degree_push_csr_fp[v];
        edge_idx = edges_idx_push_csr_fp[v];

        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array_push_csr_fp[j];
            pageRanksNext_push_csr_fp[u] += riDividedOnDiClause_push_csr_fp[v];
        }
    }
}


void pageRankPushFixedPointGraphCSRKernelCache(struct DoubleTaggedCache *cache, uint64_t *riDividedOnDiClause, uint64_t *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if((v + 1) < num_vertices)
        {
            edge_idx = edges_idx[v + 1];
            for(j = edge_idx ; j < (edge_idx + out_degree[v + 1]) ; j++)
            {
                u = sorted_edges_array[j];
                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[u])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
                }

            }

            if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[v + 1])))
            {
                Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[v + 1]), 's', (v + 1));
            }
        }
#endif

        degree = out_degree[v];
        edge_idx = edges_idx[v];

        // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree[v]), 'r', v);
        // Access(cache->accel_graph->cold_cache, (uint64_t) & (edges_idx[v]), 'r', v);

        for(j = edge_idx ; j < (edge_idx + degree) ; j++)
        {
            u = sorted_edges_array[j];

            Access(cache->accel_graph->cold_cache, (uint64_t) & (sorted_edges_array[j]), 'r', u);

            pageRanksNext[u] += riDividedOnDiClause[v];

            Access(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause[v]), 'r', v);

            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanksNext[u]), 'w', u);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanksNext[u]), 'r', u);
        }
    }
}

// ********************************************************************************************

uint32_t pageRankDataDrivenPullGraphCSRKernel(float *riDividedOnDiClause_dd_pull_csr, float *pageRanks_dd_pull_csr,
        uint32_t *in_degree_dd_pull_csr, uint32_t *in_edges_idx_dd_pull_csr, uint32_t *in_sorted_edges_array_dd_pull_csr,
        uint32_t *out_degree_dd_pull_csr, uint32_t *out_edges_idx_dd_pull_csr, uint32_t *out_sorted_edges_array_dd_pull_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;
    double base_pr = 1 - damp;


    for(v = 0; v < num_vertices; v++)
    {
        if(workListCurr[v])
        {
            double error = 0;
            float nodeIncomingPR = 0;


            degree = in_degree_dd_pull_csr[v]; // when directed we use inverse graph out degree means in degree
            edge_idx = in_edges_idx_dd_pull_csr[v];


            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = in_sorted_edges_array_dd_pull_csr[j];
                nodeIncomingPR += riDividedOnDiClause_dd_pull_csr[u]; // sum (PRi/outDegree(i))
            }

            float oldPageRank =  pageRanks_dd_pull_csr[v];
            float newPageRank =  base_pr + (damp * nodeIncomingPR);
            error = fabs(newPageRank - oldPageRank);
            (*error_total) += error / num_vertices;

            if(error >= epsilon)
            {
                pageRanks_dd_pull_csr[v] = newPageRank;
                degree = out_degree_dd_pull_csr[v];
                edge_idx = out_edges_idx_dd_pull_csr[v];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = out_sorted_edges_array_dd_pull_csr[j];
                    workListNext[u] = 1;
                }

                activeVertices++;
            }
        }
    }


    return activeVertices;
}


// ********************************************************************************************

uint32_t pageRankDataDrivenPullGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause_dd_pull_csr, float *pageRanks_dd_pull_csr,
        uint32_t *in_degree_dd_pull_csr, uint32_t *in_edges_idx_dd_pull_csr, uint32_t *in_sorted_edges_array_dd_pull_csr,
        uint32_t *out_degree_dd_pull_csr, uint32_t *out_edges_idx_dd_pull_csr, uint32_t *out_sorted_edges_array_dd_pull_csr,
        uint8_t  *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;
    double base_pr = 1 - damp;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if(workListCurr[v + 1])
            if((v + 1) < num_vertices)
            {
                degree = in_degree_dd_pull_csr[v + 1]; // when directed we use inverse graph out degree means in degree
                edge_idx = in_edges_idx_dd_pull_csr[v + 1];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = in_sorted_edges_array_dd_pull_csr[j];
                    if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause_dd_pull_csr[u])))
                    {
                        Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause_dd_pull_csr[u]), 'r', u);
                    }

                }

                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_pull_csr[v + 1])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pull_csr[v + 1]), 's', (v + 1));
                }
            }
#endif

        Access(cache->accel_graph->cold_cache, (uint64_t) & (workListCurr[v]), 'r', v);
        if(workListCurr[v])
        {

            double error = 0;
            float nodeIncomingPR = 0;

            degree = in_degree_dd_pull_csr[v]; // when directed we use inverse graph out degree means in degree
            edge_idx = in_edges_idx_dd_pull_csr[v];
            // Access(cache->accel_graph->cold_cache, (uint64_t) & (in_degree_dd_pull_csr[v]), 'r', v);
            // Access(cache->accel_graph->cold_cache, (uint64_t) & (in_edges_idx_dd_pull_csr[v]), 'r', v);


            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = in_sorted_edges_array_dd_pull_csr[j];
                Access(cache->accel_graph->cold_cache, (uint64_t) & (in_sorted_edges_array_dd_pull_csr[j]), 'r', u);


                nodeIncomingPR += riDividedOnDiClause_dd_pull_csr[u]; // sum (PRi/outDegree(i))
                Access(cache->accel_graph->cold_cache, (uint64_t) & (riDividedOnDiClause_dd_pull_csr[u]), 'r', u);
                Access(cache->accel_graph->warm_cache, (uint64_t) & (riDividedOnDiClause_dd_pull_csr[u]), 'r', u);

            }

            float oldPageRank =  pageRanks_dd_pull_csr[v];
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pull_csr[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_pull_csr[v]), 'r', v);

            float newPageRank =  base_pr + (damp * nodeIncomingPR);
            error = fabs(newPageRank - oldPageRank);

            (*error_total) += error / num_vertices;
            Access(cache->accel_graph->cold_cache, (uint64_t) & ((*error_total)), 'r', v);
            Access(cache->accel_graph->cold_cache, (uint64_t) & ((*error_total)), 'w', v);

            if(error >= epsilon)
            {
                pageRanks_dd_pull_csr[v] = newPageRank;
                Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pull_csr[v]), 'w', v);

                degree = out_degree_dd_pull_csr[v];
                edge_idx = out_edges_idx_dd_pull_csr[v];
                // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree_dd_pull_csr[v]), 'r', v);
                // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_edges_idx_dd_pull_csr[v]), 'r', v);

                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = out_sorted_edges_array_dd_pull_csr[j];
                    Access(cache->accel_graph->cold_cache, (uint64_t) & (out_sorted_edges_array_dd_pull_csr[j]), 'r', u);

                    workListNext[u] = 1;
                    Access(cache->accel_graph->cold_cache, (uint64_t) & (workListNext[u]), 'w', u);
                }

                activeVertices++;
            }
        }
    }


    return activeVertices;
}

// ********************************************************************************************

uint32_t pageRankDataDrivenPushGraphCSRKernel(float *aResiduals_dd_push_csr, float *pageRanks_dd_push_csr,
        uint32_t *out_degree_dd_push_csr, uint32_t *out_edges_idx_dd_push_csr, uint32_t *out_sorted_edges_array_dd_push_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;


    for(v = 0; v < num_vertices; v++)
    {
        if(workListCurr[v])
        {
            float oldPageRank =  pageRanks_dd_push_csr[v];
            float newPageRank =  aResiduals_dd_push_csr[v] + pageRanks_dd_push_csr[v];
            (*error_total) += fabs(newPageRank / num_vertices - oldPageRank / num_vertices);

            pageRanks_dd_push_csr[v] = newPageRank;

            degree = out_degree_dd_push_csr[v];
            float delta = damp * (aResiduals_dd_push_csr[v] / degree);

            edge_idx = out_edges_idx_dd_push_csr[v];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = out_sorted_edges_array_dd_push_csr[j];
                float prevResidual = 0.0f;

                prevResidual = aResiduals_dd_push_csr[u];

                aResiduals_dd_push_csr[u] += delta;

                if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
                {
                    activeVertices++;
                    if(!workListNext[u])
                    {
                        workListNext[u] = 1;
                    }
                }
            }
            aResiduals_dd_push_csr[v] = 0.0f;
        }
    }
    return activeVertices;
}

// ********************************************************************************************

uint32_t pageRankDataDrivenPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *aResiduals_dd_push_csr, float *pageRanks_dd_push_csr,
        uint32_t *out_degree_dd_push_csr, uint32_t *out_edges_idx_dd_push_csr, uint32_t *out_sorted_edges_array_dd_push_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if(workListCurr[v + 1])
            if((v + 1) < num_vertices)
            {
                degree = out_degree_dd_push_csr[v + 1]; // when directed we use inverse graph out degree means in degree
                edge_idx = out_edges_idx_dd_push_csr[v + 1];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = out_sorted_edges_array_dd_push_csr[j];
                    if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_push_csr[u])))
                    {
                        Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[u]), 'r', u);
                    }

                }

                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_push_csr[v + 1])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[v + 1]), 's', (v + 1));
                }
            }
#endif

        Access(cache->accel_graph->cold_cache, (uint64_t) & (workListCurr[v]), 'r', v);
        if(workListCurr[v])
        {
            float oldPageRank =  pageRanks_dd_push_csr[v];
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_push_csr[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_push_csr[v]), 'r', v);

            float newPageRank =  aResiduals_dd_push_csr[v] + pageRanks_dd_push_csr[v];
            Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_push_csr[v]), 'r', v);

            (*error_total) += fabs(newPageRank / num_vertices - oldPageRank / num_vertices);

            pageRanks_dd_push_csr[v] = newPageRank;

            degree = out_degree_dd_push_csr[v];
            float delta = damp * (aResiduals_dd_push_csr[v] / degree);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[v]), 'r', v);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree_dd_push_csr[v]), 'r', v);

            edge_idx = out_edges_idx_dd_push_csr[v];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = out_sorted_edges_array_dd_push_csr[j];
                float prevResidual = 0.0f;

                prevResidual = aResiduals_dd_push_csr[u];
                Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[u]), 'r', u);
                Access(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_push_csr[u]), 'r', u);

                aResiduals_dd_push_csr[u] += delta;
                Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_push_csr[u]), 'w', u);

                if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
                {
                    activeVertices++;
                    if(!workListNext[u])
                    {
                        workListNext[u] = 1;
                        Access(cache->accel_graph->cold_cache, (uint64_t) & (workListNext[u]), 'w', u);
                        Access(cache->accel_graph->warm_cache, (uint64_t) & (workListNext[u]), 'w', u);

                    }
                }
            }
            aResiduals_dd_push_csr[v] = 0.0f;
        }
    }
    return activeVertices;
}

// ********************************************************************************************

uint32_t pageRankDataDrivenPullPushGraphCSRKernel(float *aResiduals_dd_pullpush_csr, float *pageRanks_dd_pullpush_csr,
        uint32_t *in_degree_dd_pullpush_csr, uint32_t *in_edges_idx_dd_pullpush_csr, uint32_t *in_sorted_edges_array_dd_pullpush_csr,
        uint32_t *out_degree_dd_pullpush_csr, uint32_t *out_edges_idx_dd_pullpush_csr, uint32_t *out_sorted_edges_array_dd_pullpush_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{
    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;
    double base_pr = 1.0 - damp;


    for(v = 0; v < num_vertices; v++)
    {
        if(workListCurr[v])
        {
            float nodeIncomingPR = 0.0f;
            degree = in_degree_dd_pullpush_csr[v];
            edge_idx = in_edges_idx_dd_pullpush_csr[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = in_sorted_edges_array_dd_pullpush_csr[j];
                nodeIncomingPR += pageRanks_dd_pullpush_csr[u] / out_degree_dd_pullpush_csr[u];
            }

            float newPageRank = base_pr + (damp * nodeIncomingPR);
            float oldPageRank =  pageRanks_dd_pullpush_csr[v];

            (*error_total) += fabs(newPageRank / num_vertices - oldPageRank / num_vertices);

            pageRanks_dd_pullpush_csr[v] = newPageRank;

            degree = out_degree_dd_pullpush_csr[v];
            float delta = damp * (aResiduals_dd_pullpush_csr[v] / degree);
            edge_idx = out_edges_idx_dd_pullpush_csr[v];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = out_sorted_edges_array_dd_pullpush_csr[j];
                float prevResidual = 0.0f;

                prevResidual = aResiduals_dd_pullpush_csr[u];

                aResiduals_dd_pullpush_csr[u] += delta;

                if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
                {
                    activeVertices++;
                    if(!workListNext[u])
                    {
                        workListNext[u] = 1;
                    }
                }
            }
            aResiduals_dd_pullpush_csr[v] = 0.0f;
        }
    }
    return activeVertices;
}

// ********************************************************************************************

uint32_t pageRankDataDrivenPullPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *aResiduals_dd_pullpush_csr, float *pageRanks_dd_pullpush_csr,
        uint32_t *in_degree_dd_pullpush_csr, uint32_t *in_edges_idx_dd_pullpush_csr, uint32_t *in_sorted_edges_array_dd_pullpush_csr,
        uint32_t *out_degree_dd_pullpush_csr, uint32_t *out_edges_idx_dd_pullpush_csr, uint32_t *out_sorted_edges_array_dd_pullpush_csr,
        uint8_t  *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices)
{

    uint32_t j;
    uint32_t v;
    uint32_t u;
    uint32_t degree;
    uint32_t edge_idx;
    uint32_t activeVertices = 0;
    double damp = 0.85;
    double base_pr = 1 - damp;

    for(v = 0; v < num_vertices; v++)
    {

#ifdef PREFETCH
        if(workListCurr[v + 1])
            if((v + 1) < num_vertices)
            {
                degree = in_degree_dd_pullpush_csr[v + 1]; // when directed we use inverse graph out degree means in degree
                edge_idx = in_edges_idx_dd_pullpush_csr[v + 1];
                for(j = edge_idx ; j < (edge_idx + degree) ; j++)
                {
                    u = in_sorted_edges_array_dd_pullpush_csr[j];
                    if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[u])))
                    {
                        Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[u]), 'r', u);
                    }

                    if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (out_degree_dd_pullpush_csr[u])))
                    {
                        Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (out_degree_dd_pullpush_csr[u]), 'r', u);
                    }

                }

                if(checkInCache(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v + 1])))
                {
                    Prefetch(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v + 1]), 's', (v + 1));
                }
            }
#endif

        Access(cache->accel_graph->cold_cache, (uint64_t) & (workListCurr[v]), 'r', v);
        if(workListCurr[v])
        {
            float nodeIncomingPR = 0.0f;
            degree = in_degree_dd_pullpush_csr[v];
            edge_idx = in_edges_idx_dd_pullpush_csr[v];
            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = in_sorted_edges_array_dd_pullpush_csr[j];
                // Access(cache->accel_graph->cold_cache, (uint64_t) & (in_sorted_edges_array_dd_pullpush_csr[j]), 'r', u);

                nodeIncomingPR += pageRanks_dd_pullpush_csr[u] / out_degree_dd_pullpush_csr[u];
                Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[u]), 'r', u);
                Access(cache->accel_graph->cold_cache, (uint64_t) & (out_degree_dd_pullpush_csr[u]), 'r', u);

                Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[u]), 'r', u);
                Access(cache->accel_graph->warm_cache, (uint64_t) & (out_degree_dd_pullpush_csr[u]), 'r', u);
            }

            float newPageRank = base_pr + (damp * nodeIncomingPR);
            float oldPageRank =  pageRanks_dd_pullpush_csr[v];
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[v]), 'r', v);

            (*error_total) += fabs(newPageRank / num_vertices - oldPageRank / num_vertices);

            pageRanks_dd_pullpush_csr[v] = newPageRank;
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[v]), 'r', v);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (pageRanks_dd_pullpush_csr[v]), 'w', v);

            degree = out_degree_dd_pullpush_csr[v];
            float delta = damp * (aResiduals_dd_pullpush_csr[v] / degree);
            Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v]), 'r', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v]), 'r', v);

            edge_idx = out_edges_idx_dd_pullpush_csr[v];

            for(j = edge_idx ; j < (edge_idx + degree) ; j++)
            {
                u = out_sorted_edges_array_dd_pullpush_csr[j];
                // Access(cache->accel_graph->cold_cache, (uint64_t) & (out_sorted_edges_array_dd_pullpush_csr[j]), 'r', u);

                float prevResidual = 0.0f;

                prevResidual = aResiduals_dd_pullpush_csr[u];
                Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[u]), 'r', u);
                Access(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[u]), 'r', u);

                aResiduals_dd_pullpush_csr[u] += delta;
                Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[u]), 'w', u);


                if ((fabs(prevResidual + delta) >= epsilon) && (prevResidual <= epsilon))
                {
                    activeVertices++;

                    Access(cache->accel_graph->cold_cache, (uint64_t) & (workListNext[u]), 'r', u);
                    Access(cache->accel_graph->warm_cache, (uint64_t) & (workListNext[u]), 'r', u);
                    if(!workListNext[u])
                    {
                        workListNext[u] = 1;
                        Access(cache->accel_graph->cold_cache, (uint64_t) & (workListNext[u]), 'w', u);
                    }
                }
            }
            aResiduals_dd_pullpush_csr[v] = 0.0f;
            Access(cache->accel_graph->cold_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v]), 'w', v);
            Access(cache->accel_graph->warm_cache, (uint64_t) & (aResiduals_dd_pullpush_csr[v]), 'w', v);
        }
    }

    return activeVertices;

}

// ********************************************************************************************
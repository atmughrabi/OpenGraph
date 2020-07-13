#ifndef PAGERANK_KERNELS_H
#define PAGERANK_KERNELS_H

#include <stdint.h>
#include "cache.h"

// ********************************************************************************************
// ***************          GRID DataStructure               **************
// ********************************************************************************************

void pageRankPullRowGraphGridKernel(float *riDividedOnDiClause_pull_grid, float *pageRanksNext_pull_grid,  struct Partition *partitions, uint32_t totalPartitions);
// ********************************************************************************************
void pageRankPushColumnGraphGridKernel(float *riDividedOnDiClause_push_grid, float *pageRanksNext_push_grid,  struct Partition *partitions, uint32_t totalPartitions);
// ********************************************************************************************
void pageRankPullRowFixedPointGraphGridKernel(uint64_t *riDividedOnDiClause_push_grid_fp, uint64_t *pageRanksNext_push_grid_fp,  struct Partition *partitions, uint32_t totalPartitions);
// ********************************************************************************************
void pageRankPushColumnFixedPointGraphGridKernel(uint64_t *riDividedOnDiClause_push_grid_fp, uint64_t *pageRanksNext_push_grid_fp,  struct Partition *partitions, uint32_t totalPartitions);

// ********************************************************************************************
// ***************          CSR DataStructure                                    **************
// ********************************************************************************************

void pageRankPullGraphCSRKernel(float *riDividedOnDiClause_pull_csr, float *pageRanksNext_pull_csr, uint32_t *out_degree_pull_csr, uint32_t *edges_idx_pull_csr, uint32_t *sorted_edges_array_pull_csr, uint32_t num_vertices);
void pageRankPullGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause, float *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices);
// ********************************************************************************************
void pageRankPushGraphCSRKernel(float *riDividedOnDiClause_push_csr, float *pageRanksNext_push_csr, uint32_t *out_degree_push_csr, uint32_t *edges_idx_push_csr, uint32_t *sorted_edges_array_push_csr, uint32_t num_vertices);
void pageRankPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause, float *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices);
// ********************************************************************************************
void pageRankPullFixedPointGraphCSRKernel(uint64_t *riDividedOnDiClause_pull_csr_fp, uint64_t *pageRanksNext_pull_csr_fp, uint32_t *out_degree_pull_csr_fp, uint32_t *edges_idx_pull_csr_fp, uint32_t *sorted_edges_array_pull_csr_fp, uint32_t num_vertices);
void pageRankPullFixedPointGraphCSRKernelCache(struct DoubleTaggedCache *cache, uint64_t *riDividedOnDiClause, uint64_t *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices);
// ********************************************************************************************
void pageRankPushFixedPointGraphCSRKernel(uint64_t *riDividedOnDiClause_push_csr_fp, uint64_t *pageRanksNext_push_csr_fp, uint32_t *out_degree_push_csr_fp, uint32_t *edges_idx_push_csr_fp, uint32_t *sorted_edges_array_push_csr_fp, uint32_t num_vertices);
void pageRankPushFixedPointGraphCSRKernelCache(struct DoubleTaggedCache *cache, uint64_t *riDividedOnDiClause, uint64_t *pageRanksNext, uint32_t *out_degree, uint32_t *edges_idx, uint32_t *sorted_edges_array, uint32_t num_vertices);
// ********************************************************************************************
uint32_t pageRankDataDrivenPullGraphCSRKernel(float *riDividedOnDiClause_dd_pull_csr, float *pageRanks_dd_pull_csr,
        uint32_t *in_degree_dd_pull_csr, uint32_t *in_edges_idx_dd_pull_csr, uint32_t *in_sorted_edges_array_dd_pull_csr,
        uint32_t *out_degree_dd_pull_csr, uint32_t *out_edges_idx_dd_pull_csr, uint32_t *out_sorted_edges_array_dd_pull_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
uint32_t pageRankDataDrivenPullGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause_dd_pull_csr, float *pageRanks_dd_pull_csr,
        uint32_t *in_degree_dd_pull_csr, uint32_t *in_edges_idx_dd_pull_csr, uint32_t *in_sorted_edges_array_dd_pull_csr,
        uint32_t *out_degree_dd_pull_csr, uint32_t *out_edges_idx_dd_pull_csr, uint32_t *out_sorted_edges_array_dd_pull_csr,
        uint8_t  *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
// ********************************************************************************************
uint32_t pageRankDataDrivenPushGraphCSRKernel(float *aResiduals_dd_push_csr, float *pageRanks_dd_push_csr,
        uint32_t *out_degree_dd_push_csr, uint32_t *out_edges_idx_dd_push_csr, uint32_t *out_sorted_edges_array_dd_push_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
uint32_t pageRankDataDrivenPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *aResiduals_dd_push_csr, float *pageRanks_dd_push_csr,
        uint32_t *out_degree_dd_push_csr, uint32_t *out_edges_idx_dd_push_csr, uint32_t *out_sorted_edges_array_dd_push_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
// ********************************************************************************************
uint32_t pageRankDataDrivenPullPushGraphCSRKernel(float *riDividedOnDiClause_dd_pullpush_csr, float *pageRanks_dd_pullpush_csr,
        uint32_t *in_degree_dd_pullpush_csr, uint32_t *in_edges_idx_dd_pullpush_csr, uint32_t *in_sorted_edges_array_dd_pullpush_csr,
        uint32_t *out_degree_dd_pullpush_csr, uint32_t *out_edges_idx_dd_pullpush_csr, uint32_t *out_sorted_edges_array_dd_pullpush_csr,
        uint8_t *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
uint32_t pageRankDataDrivenPullPushGraphCSRKernelCache(struct DoubleTaggedCache *cache, float *riDividedOnDiClause_dd_pullpush_csr, float *pageRanks_dd_pullpush_csr,
        uint32_t *in_degree_dd_pullpush_csr, uint32_t *in_edges_idx_dd_pullpush_csr, uint32_t *in_sorted_edges_array_dd_pullpush_csr,
        uint32_t *out_degree_dd_pullpush_csr, uint32_t *out_edges_idx_dd_pullpush_csr, uint32_t *out_sorted_edges_array_dd_pullpush_csr,
        uint8_t  *workListCurr, uint8_t *workListNext, double *error_total, double epsilon, uint32_t num_vertices);
// ********************************************************************************************


#endif
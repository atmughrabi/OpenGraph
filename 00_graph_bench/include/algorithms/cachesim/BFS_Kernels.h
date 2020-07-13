#ifndef BFS_KERNELS_H
#define BFS_KERNELS_H

#include <stdint.h>

#include "bitmap.h"
#include "cache.h"



// ********************************************************************************************
// ***************          CSR DataStructure                                    **************
// ********************************************************************************************

// uint32_t topDownStepGraphCSRKernelAladdin(struct GraphCSR *graph, struct ArrayQueue *sharedFrontierQueue,  struct ArrayQueue **localFrontierQueues, struct BFSStats *stats);
void bottomUpStepGraphCSRKernel(  uint32_t *nf, int *parents,  uint32_t *distances, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, uint32_t *out_degree_pull_csr, uint32_t *edges_idx_pull_csr, uint32_t *sorted_edges_array_pull_csr, uint32_t num_vertices);

#endif
// -----------------------------------------------------------------------------
//
//		"00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : BFS_Kernels.c
// Create : 2019-10-11 16:26:36
// Revise : 2019-10-11 21:39:42
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "myMalloc.h"

#include "bitmap.h"
#include "cache.h"
#include "BFS_Kernels.h"

// you should add these to Aladdin as an extern "aladdin_sys_constants.h"
// unsigned ACCELGRAPH = 0x300;

void bottomUpStepGraphCSRKernel( uint32_t *nf, int *parents,  uint32_t *distances, struct Bitmap *bitmapCurr, struct Bitmap *bitmapNext, uint32_t *out_degree_pull_csr, uint32_t *edges_idx_pull_csr, uint32_t *sorted_edges_array_pull_csr, uint32_t num_vertices)
{

	uint32_t j;
 	uint32_t v;
 	uint32_t u;
    uint32_t edge_idx;
    uint32_t out_degree;

iter :	  
    for(v = 0 ; v < num_vertices ; v++)
    {
        out_degree = out_degree_pull_csr[v];
        if(parents[v] < 0)  // optmization
        {
            edge_idx = edges_idx_pull_csr[v];

            for(j = edge_idx ; j < (edge_idx + out_degree) ; j++)
            {
                u = sorted_edges_array_pull_csr[j];
                if(getBit(bitmapCurr, u))
                {
                    parents[v] = u;
                    distances[v] = distances[u] + 1;
                    setBit(bitmapNext, v);
                    (*nf)++;
                    break;
                }
            }

        }

    }
 
}
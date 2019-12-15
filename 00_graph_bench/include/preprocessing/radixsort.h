#ifndef RADIXSORT_H
#define RADIXSORT_H

#include <linux/types.h>
#include "edgeList.h"



struct EdgeList* radixSortEdgesBySource (struct EdgeList* edgeList);
struct EdgeList* radixSortEdgesBySourceAndDestination (struct EdgeList* edgeList);
// struct EdgeList* radixSortEdgesBySourceOptimized (struct EdgeList* edgeList);
// struct EdgeList* radixSortEdgesBySourceOptimizedParallel (struct EdgeList* edgeList);
void radixSortCountSortEdgesBySource (struct EdgeList** sorted_edges_array, struct EdgeList** edgeList, __u32 radix, __u32 buckets, __u32* buckets_count);
void radixSortCountSortEdgesByDestination (struct EdgeList** sorted_edges_array, struct EdgeList** edgeList, __u32 radix, __u32 buckets, __u32* buckets_count);

#endif
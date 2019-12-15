#ifndef REORDER_H
#define REORDER_H

#include <linux/types.h>
#include "edgeList.h"
#include "graphCSR.h"

struct EdgeList *relabelEdgeListFromFile(struct EdgeList *edgeList, const char *fnameb, __u32 size);
void writeLabelsToFile(const char *fnameb, __u32 *labels, __u32 size);
struct EdgeList *relabelEdgeList(struct EdgeList *edgeList, __u32 *labels);
struct EdgeList *reorderGraphProcess(struct EdgeList *edgeList, struct Arguments *arguments);
struct EdgeList *reorderGraphProcessPageRank(struct EdgeList *edgeList, struct Arguments *arguments);
struct EdgeList *reorderGraphProcessDegree( __u32 sort, struct EdgeList *edgeList, __u32 lmode);
struct EdgeList *reorderGraphListDegree(struct EdgeList *edgeList, __u32 *degrees, __u32 lmode);
struct EdgeList *reorderGraphListPageRank(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochPageRank(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochBFS(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochRabbit(struct GraphCSR *graph);

__u32 *reorderGraphProcessInOutDegrees(__u32 *degrees, struct EdgeList *edgeList, __u32 lmode);
__u32 reorderGraphProcessVertexSize( struct EdgeList *edgeList);
__u32 *radixSortEdgesByDegree (__u32 *degrees, __u32 *labels, __u32 num_vertices);
__u32 *radixSortEdgesByPageRank (float *pageRanks, __u32 *labels, __u32 num_vertices);

void radixSortCountSortEdges(__u32 **pageRanksFP, __u32 **pageRanksFPTemp, __u32 **labels, __u32 **labelsTemp, __u32 radix, __u32 buckets, __u32 *buckets_count, __u32 num_vertices);

#endif
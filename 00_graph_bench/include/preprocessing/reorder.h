#ifndef REORDER_H
#define REORDER_H

#include <stdint.h>
#include "edgeList.h"
#include "graphCSR.h"

struct EdgeList *relabelEdgeListFromFile(struct EdgeList *edgeList, const char *fnameb, uint32_t size);
void writeLabelsToFile(const char *fnameb, uint32_t *labels, uint32_t size);
struct EdgeList *relabelEdgeList(struct EdgeList *edgeList, uint32_t *labels);
struct EdgeList *reorderGraphProcess(struct EdgeList *edgeList, struct Arguments *arguments);
struct EdgeList *reorderGraphProcessPageRank(struct EdgeList *edgeList, struct Arguments *arguments);
struct EdgeList *reorderGraphProcessDegree( uint32_t sort, struct EdgeList *edgeList, uint32_t lmode);
struct EdgeList *reorderGraphListDegree(struct EdgeList *edgeList, uint32_t *degrees, uint32_t lmode);
struct EdgeList *reorderGraphListPageRank(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochPageRank(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochBFS(struct GraphCSR *graph);
struct EdgeList *reorderGraphListEpochRabbit(struct GraphCSR *graph);

uint32_t *reorderGraphProcessInOutDegrees(uint32_t *degrees, struct EdgeList *edgeList, uint32_t lmode);
uint32_t reorderGraphProcessVertexSize( struct EdgeList *edgeList);
uint32_t *radixSortEdgesByDegree (uint32_t *degrees, uint32_t *labels, uint32_t num_vertices);
uint32_t *radixSortEdgesByPageRank (float *pageRanks, uint32_t *labels, uint32_t num_vertices);

void radixSortCountSortEdges(uint32_t **pageRanksFP, uint32_t **pageRanksFPTemp, uint32_t **labels, uint32_t **labelsTemp, uint32_t radix, uint32_t buckets, uint32_t *buckets_count, uint32_t num_vertices);

#endif
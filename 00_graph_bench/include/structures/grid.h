#ifndef GRID_H
#define GRID_H

#include <stdint.h>
#include "bitmap.h"
#include "graphConfig.h"
#include "edgeList.h"
#include "vertex.h"



// A structure to represent an adjacency list
struct  Partition
{

    // uint32_t vertex_start;
    // uint32_t vertex_end;
    uint32_t num_edges;
    uint32_t num_vertices;
    struct EdgeList *edgeList;

};


// A structure to represent an adjacency list
struct  Grid
{

    uint32_t num_edges;
    uint32_t num_vertices;
    uint32_t num_partitions;
    struct Partition *partitions;
    uint32_t *activePartitions;
    uint32_t *out_degree;
    uint32_t *in_degree;
    struct Bitmap *activePartitionsMap;
    // struct Bitmap* activeVertices;
};


void gridPrint(struct Grid *grid);
struct Grid *gridNew(struct EdgeList *edgeList, uint32_t cache_size);
void  gridFree(struct Grid *grid);
void gridPrintMessageWithtime(const char *msg, double time);

struct Grid *gridPartitionEdgeListSizePreprocessing(struct Grid *grid, struct EdgeList *edgeList);
struct Grid *gridPartitionVertexSizePreprocessing(struct Grid *grid);
uint32_t gridCalculatePartitions(struct EdgeList *edgeList, uint32_t cache_size);
struct Grid *gridPartitionsMemoryAllocations(struct Grid *grid);
struct Grid *gridPartitionEdgePopulation(struct Grid *grid, struct EdgeList *edgeList);
struct Grid *graphGridProcessInOutDegrees(struct Grid *grid, struct EdgeList *edgeList);
void   graphGridSetActivePartitions(struct Grid *grid, uint32_t vertex);
void   graphGridResetActivePartitions(struct Grid *grid);

void   graphGridResetActivePartitionsMap(struct Grid *grid);
void   graphGridSetActivePartitionsMap(struct Grid *grid, uint32_t vertex);

uint32_t getPartitionID(uint32_t vertices, uint32_t partitions, uint32_t vertex_id);

#endif
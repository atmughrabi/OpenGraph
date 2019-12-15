#ifndef GRID_H
#define GRID_H

#include <linux/types.h>
#include "bitmap.h"
#include "graphConfig.h"
#include "edgeList.h"
#include "vertex.h"



// A structure to represent an adjacency list
struct  Partition
{

    // __u32 vertex_start;
    // __u32 vertex_end;
    __u32 num_edges;
    __u32 num_vertices;
    struct EdgeList *edgeList;

};


// A structure to represent an adjacency list
struct  Grid
{

    __u32 num_edges;
    __u32 num_vertices;
    __u32 num_partitions;
    struct Partition *partitions;
    __u32 *activePartitions;
    __u32 *out_degree;
    __u32 *in_degree;
    struct Bitmap *activePartitionsMap;
    // struct Bitmap* activeVertices;
};


void gridPrint(struct Grid *grid);
struct Grid *gridNew(struct EdgeList *edgeList);
void  gridFree(struct Grid *grid);
void gridPrintMessageWithtime(const char *msg, double time);

struct Grid *gridPartitionEdgeListSizePreprocessing(struct Grid *grid, struct EdgeList *edgeList);
struct Grid *gridPartitionVertexSizePreprocessing(struct Grid *grid);
__u32 gridCalculatePartitions(struct EdgeList *edgeList);
struct Grid *gridPartitionsMemoryAllocations(struct Grid *grid);
struct Grid *gridPartitionEdgePopulation(struct Grid *grid, struct EdgeList *edgeList);
struct Grid *graphGridProcessInOutDegrees(struct Grid *grid, struct EdgeList *edgeList);
void   graphGridSetActivePartitions(struct Grid *grid, __u32 vertex);
void   graphGridResetActivePartitions(struct Grid *grid);

void   graphGridResetActivePartitionsMap(struct Grid *grid);
void   graphGridSetActivePartitionsMap(struct Grid *grid, __u32 vertex);

__u32 getPartitionID(__u32 vertices, __u32 partitions, __u32 vertex_id);

#endif
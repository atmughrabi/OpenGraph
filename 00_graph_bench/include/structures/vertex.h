#ifndef VERTEX_H
#define VERTEX_H

#include <linux/types.h>
#include "graphCSR.h"

#define NO_OUTGOING_EDGES -1
#define NO_INCOMING_EDGES 0
#define NOT_VISITED 0



struct  Vertex
{
    __u32 *out_degree;
    __u32 *in_degree;
    __u32 *edges_idx;
};


struct GraphCSR *mapVertices (struct GraphCSR *graph, __u8 inverse);
struct GraphCSR *mapVerticesWithInOutDegree (struct GraphCSR *graph, __u8 inverse);

struct Vertex *newVertexArray(__u32 num_vertices);
void freeVertexArray(struct Vertex *vertices);
void printVertexArray(struct Vertex *vertex_array, __u32 num_vertices);
void vertexArrayMaxOutdegree(struct Vertex *vertex_array, __u32 num_vertices);
void vertexArrayMaxInDegree(struct Vertex *vertex_array, __u32 num_vertices);
void partitionEdgeListOffsetStartEnd(struct GraphCSR *graph, struct EdgeList *sorted_edges_array, __u32 *offset_start, __u32 *offset_end);


#endif
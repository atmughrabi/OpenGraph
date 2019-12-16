#ifndef VERTEX_H
#define VERTEX_H

#include <stdint.h>
#include "graphCSR.h"

#define NO_OUTGOING_EDGES -1
#define NO_INCOMING_EDGES 0
#define NOT_VISITED 0



struct  Vertex
{
    uint32_t *out_degree;
    uint32_t *in_degree;
    uint32_t *edges_idx;
};


struct GraphCSR *mapVertices (struct GraphCSR *graph, uint8_t inverse);
struct GraphCSR *mapVerticesWithInOutDegree (struct GraphCSR *graph, uint8_t inverse);

struct Vertex *newVertexArray(uint32_t num_vertices);
void freeVertexArray(struct Vertex *vertices);
void printVertexArray(struct Vertex *vertex_array, uint32_t num_vertices);
void vertexArrayMaxOutdegree(struct Vertex *vertex_array, uint32_t num_vertices);
void vertexArrayMaxInDegree(struct Vertex *vertex_array, uint32_t num_vertices);
void partitionEdgeListOffsetStartEnd(struct GraphCSR *graph, struct EdgeList *sorted_edges_array, uint32_t *offset_start, uint32_t *offset_end);


#endif
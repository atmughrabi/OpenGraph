#ifndef GRAPHCSR_H
#define GRAPHCSR_H

#include <stdint.h>

#include "graphConfig.h"
#include "edgeList.h"
#include "vertex.h"



struct GraphCSR
{

    uint32_t num_edges;
    uint32_t num_vertices;
    uint32_t avg_degree;
#if WEIGHTED
    float max_weight;
#endif

    struct Vertex *vertices;
    struct EdgeList *sorted_edges_array; // sorted edge array

#if DIRECTED
    struct Vertex *inverse_vertices;
    struct EdgeList *inverse_sorted_edges_array; // sorted edge array
#endif
};

void graphCSRFree (struct GraphCSR *graphCSR);
void graphCSRPrint (struct GraphCSR *graphCSR);
struct GraphCSR *graphCSRAssignEdgeList (struct GraphCSR *graphCSR, struct EdgeList *edgeList, uint8_t inverse);
struct GraphCSR *graphCSRNew(uint32_t V, uint32_t E,  uint8_t inverse);
struct GraphCSR *graphCSRPreProcessingStep (struct Arguments *arguments);
struct GraphCSR *graphCSRPreProcessingStepDualOrder (struct Arguments *arguments);
void graphCSRVertexLabelRemappingDualOrder (struct GraphCSR *graphCSR);
void graphCSRPrintMessageWithtime(const char *msg, double time);
struct GraphCSR *readFromBinFileGraphCSR (const char *fname);
void writeToBinFileGraphCSR (const char *fname, struct GraphCSR *graph);

#endif
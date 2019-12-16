#ifndef GRAPHSTATS_H
#define GRAPHSTATS_H

#include <stdint.h>
#include "graphConfig.h"
#include "graphCSR.h"

void collectStats( struct Arguments *arguments);
void countHistogram(struct GraphCSR *graphStats, uint32_t *histogram, uint32_t binSize, uint32_t inout_degree);
void printHistogram(const char *fname_stats, uint32_t *histogram, uint32_t binSize);
void printSparseMatrixList(const char *fname_stats, struct GraphCSR *graphStats, uint32_t binSize);

#endif


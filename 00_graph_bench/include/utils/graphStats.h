#ifndef GRAPHSTATS_H
#define GRAPHSTATS_H

#include <linux/types.h>
#include "graphConfig.h"
#include "graphCSR.h"

void collectStats( struct Arguments *arguments);
void countHistogram(struct GraphCSR *graphStats, __u32 *histogram, __u32 binSize, __u32 inout_degree);
void printHistogram(const char *fname_stats, __u32 *histogram, __u32 binSize);
void printSparseMatrixList(const char *fname_stats, struct GraphCSR *graphStats, __u32 binSize);

#endif


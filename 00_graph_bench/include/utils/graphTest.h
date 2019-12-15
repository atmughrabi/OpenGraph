#ifndef GRAPHTEST_H
#define GRAPHTEST_H


#include <linux/types.h>
#include "graphConfig.h"


__u32 compareRealRanks(__u32 *arr1, __u32 *arr2, __u32 arr1_size, __u32 arr2_size);
__u32 cmpGraphAlgorithmsTestStats(void *ref_stats, void *cmp_stats, __u32 algorithm);
__u32 compareDistanceArrays(__u32 *arr1, __u32 *arr2, __u32 arr1_size, __u32 arr2_size);
void *runGraphAlgorithmsTest(void *graph, struct Arguments *arguments);
__u32 equalFloat(float a, float b, float epsilon);
__u32 compareFloatArrays(float *arr1, float *arr2, __u32 arr1_size, __u32 arr2_size);

#endif
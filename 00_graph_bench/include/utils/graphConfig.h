#ifndef GRAPHCONFIG_H
#define GRAPHCONFIG_H

#include <linux/types.h>
#include "mt19937.h"


#define WEIGHTED 1
#define DIRECTED 1

extern int numThreads;
extern mt19937state *mt19937var;

/* Used by main to communicate with parse_opt. */
struct Arguments
{
    int wflag;
    int xflag;
    int sflag;
    int dflag;
    __u32 binSize;
    __u32 inout_degree;
    __u32 iterations;
    __u32 trials;
    double epsilon;
    int root;
    __u32 algorithm;
    __u32 datastructure;
    __u32 pushpull;
    __u32 sort;
    __u32 lmode;
    __u32 symmetric;
    __u32 weighted;
    __u32 delta;
    __u32 numThreads;
    char *fnameb;
    __u32 fnameb_format;
    __u32 convert_format;
};

#endif
#ifndef SORTRUN_H
#define SORTRUN_H

#include <linux/types.h>
#include "edgeList.h"

struct EdgeList *sortRunAlgorithms(struct EdgeList *edgeList, __u32 sort);
void sortRunPrintMessageWithtime(const char *msg, double time);

#endif
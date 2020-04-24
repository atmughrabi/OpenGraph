#ifndef SORTRUN_H
#define SORTRUN_H

#include <stdint.h>
#include "edgeList.h"

struct EdgeList *sortRunAlgorithms(struct EdgeList *edgeList, uint32_t sort);
void sortRunPrintMessageWithtime(const char *msg, double time);

#endif
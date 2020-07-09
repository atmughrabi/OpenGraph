#ifndef WORKLIST_H
#define WORKLIST_H

#include <stdint.h>
#include "graphConfig.h"


void swapWorkLists (uint8_t **workList1, uint8_t **workList2);
void resetWorkList(uint8_t *workList, uint32_t size);
void setWorkList(uint8_t *workList,  uint32_t size);

// int main()
// {
//     char mybits[(BITARRAY_BITS + 7) / 8];
//     memset(mybits, 0, sizeof(mybits));

//     ba_setbit(mybits, 33, 1);
//     if (!ba_get(33))
//         return 1;
//     return 0;
// };


#endif
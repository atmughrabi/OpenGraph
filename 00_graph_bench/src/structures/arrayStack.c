// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : arrayStack.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <linux/types.h>
#include <omp.h>

#include "myMalloc.h"
#include "arrayStack.h"
#include "bitmap.h"

struct ArrayStack *newArrayStack(__u32 size)
{

    struct ArrayStack *arrayStack = (struct ArrayStack *) my_malloc( sizeof(struct ArrayStack));

    arrayStack->head = 0;
    arrayStack->tail = 0;
    arrayStack->tail_next = 0;
    arrayStack->size = size;

    arrayStack->Stack = (__u32 *) my_malloc(size * sizeof(__u32));


    arrayStack->q_bitmap = newBitmap(size);

    arrayStack->q_bitmap_next = newBitmap(size);

    return arrayStack;

}


void resetArrayStack(struct ArrayStack *q)
{

    q->head = 0;
    q->tail = 0;
    q->tail_next = 0;
    clearBitmap(q->q_bitmap);

}

void freeArrayStack(struct ArrayStack *q)
{

    if(q)
    {
        if(q->q_bitmap_next)
            freeBitmap(q->q_bitmap_next);
        if(q->q_bitmap)
            freeBitmap(q->q_bitmap);
        if(q->Stack)
            free(q->Stack);
        free(q);
    }


}

void pushArrayStack (struct ArrayStack *q, __u32 k)
{

    q->Stack[q->tail] = k;
    q->tail = (q->tail + 1) % q->size;
    q->tail_next = q->tail;

}


void pushArrayStackWithBitmap (struct ArrayStack *q, __u32 k)
{

    q->Stack[q->tail] = k;
    setBit(q->q_bitmap, k);
    q->tail = q->tail_next;
    q->tail++;
    q->tail_next++;

}


void pushArrayStackAtomic (struct ArrayStack *q, __u32 k)
{

    __u32 local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->Stack[local_q_tail] = k;

}


void pushArrayStackWithBitmapAtomic (struct ArrayStack *q, __u32 k)
{

    __u32 local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->Stack[local_q_tail] = k;
    setBitAtomic(q->q_bitmap, k);

}


void pushArrayStackDelayed (struct ArrayStack *q, __u32 k)
{

    q->Stack[q->tail_next] = k;
    q->tail_next++;

}

void pushArrayStackDelayedWithBitmap (struct ArrayStack *q, __u32 k)
{

    q->Stack[q->tail_next] = k;
    setBit(q->q_bitmap_next, k);
    q->tail_next++;

}

void pushArrayStackDelayedWithBitmapAtomic (struct ArrayStack *q, __u32 k)
{

    __u32 local_q_tail_next = __sync_fetch_and_add(&q->tail_next, 1);
    setBitAtomic(q->q_bitmap, k);
    q->Stack[local_q_tail_next] = k;

}


void slideWindowArrayStack (struct ArrayStack *q)
{

    q->head = q->tail;
    q->tail = q->tail_next;

}

void slideWindowArrayStackBitmap (struct ArrayStack *q)
{

    q->head = q->tail;
    q->tail = q->tail_next;
    swapBitmaps(&q->q_bitmap, &q->q_bitmap_next);
    clearBitmap(q->q_bitmap_next);

}

__u32 popArrayStack(struct ArrayStack *q)
{

    __u32 k = q->Stack[q->tail - 1];
    clearBit(q->q_bitmap, k);
    q->tail = q->tail - 1;

    return k;

}


__u32 frontArrayStack (struct ArrayStack *q)
{

    __u32 k = q->Stack[q->head];

    return k;

}

__u8 isEmptyArrayStackCurr (struct ArrayStack *q)
{

    if((q->tail > q->head))
        return 0;
    else
        return 1;

}

__u8 isEmptyArrayStack (struct ArrayStack *q)
{

    if(!isEmptyArrayStackCurr(q) || !isEmptyArrayStackNext(q))
        return 0;
    else
        return 1;

}

__u8 isEmptyArrayStackNext (struct ArrayStack *q)
{

    if((q->tail_next > q->head))
        return 0;
    else
        return 1;

}

__u8  ispushArrayStack  (struct ArrayStack *q, __u32 k)
{


    return getBit(q->q_bitmap, k);

}

__u8  ispushArrayStackNext  (struct ArrayStack *q, __u32 k)
{


    return getBit(q->q_bitmap_next, k);

}

__u32 sizeArrayStackCurr(struct ArrayStack *q)
{

    return q->tail - q->head;

}

__u32 sizeArrayStackNext(struct ArrayStack *q)
{

    return q->tail_next - q->tail;
}


__u32 sizeArrayStack(struct ArrayStack *q)
{

    return q->tail_next - q->head;

}

void flushArrayStackToShared(struct ArrayStack *local_q, struct ArrayStack *shared_q)
{

    __u32 shared_q_tail_next = __sync_fetch_and_add(&shared_q->tail_next, local_q->tail);
    __u32 local_q_size = local_q->tail - local_q->head;

    memcpy(&shared_q->Stack[shared_q_tail_next], &local_q->Stack[local_q->head], local_q_size * (sizeof(__u32)));

    local_q->head = 0;
    local_q->tail = 0;
    local_q->tail_next = 0;

}



void arrayStackGenerateBitmap(struct ArrayStack *q)
{

    __u32 v;
    __u32 i;

    #pragma omp parallel for
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->Stack[i];
        setBitAtomic(q->q_bitmap, v);
    }

}


void arrayStackToBitmap(struct ArrayStack *q, struct Bitmap *b)
{

    __u32 v;
    __u32 i;

    #pragma omp parallel for default(none) shared(q,b) private(v,i)
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->Stack[i];
        setBitAtomic(b, v);
    }

    // b->numSetBits = q->q_bitmap->numSetBits;
    q->head = q->tail;
    q->tail_next = q->tail;


}

void bitmapToArrayStack(struct Bitmap *b, struct ArrayStack *q, struct ArrayStack **localFrontierStacks)
{

    #pragma omp parallel default(none) shared(b,localFrontierStacks,q)
    {
        __u32 i;

        __u32 t_id = omp_get_thread_num();
        struct ArrayStack *localFrontierStack = localFrontierStacks[t_id];

        #pragma omp for
        for(i = 0 ; i < (b->size); i++)
        {
            if(getBit(b, i))
            {
                localFrontierStack->Stack[localFrontierStack->tail] = i;
                localFrontierStack->tail++;
            }

        }

        flushArrayStackToShared(localFrontierStack, q);

    }

    slideWindowArrayStack(q);

}

// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : arrayQueue.c
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
#include "arrayQueue.h"
#include "bitmap.h"

struct ArrayQueue *newArrayQueue(__u32 size)
{

    struct ArrayQueue *arrayQueue = (struct ArrayQueue *) my_malloc( sizeof(struct ArrayQueue));


    arrayQueue->head = 0;
    arrayQueue->tail = 0;
    arrayQueue->tail_next = 0;
    arrayQueue->size = size;

    arrayQueue->queue = (__u32 *) my_malloc(size * sizeof(__u32));

    arrayQueue->q_bitmap = newBitmap(size);

    arrayQueue->q_bitmap_next = newBitmap(size);

    return arrayQueue;

}

void softResetArrayQueue(struct ArrayQueue *q)
{

    q->head = 0;
    q->tail = 0;
    q->tail_next = 0;
    // clearBitmap(q->q_bitmap);
    // clearBitmap(q->q_bitmap_next);

}

void resetArrayQueue(struct ArrayQueue *q)
{

    q->head = 0;
    q->tail = 0;
    q->tail_next = 0;
    clearBitmap(q->q_bitmap);
    // clearBitmap(q->q_bitmap_next);

}

void freeArrayQueue(struct ArrayQueue *q)
{
    if(q)
    {
        if(q->q_bitmap_next)
            freeBitmap(q->q_bitmap_next);
        if(q->q_bitmap)
            freeBitmap(q->q_bitmap);
        if(q->queue)
            free(q->queue);
        free(q);
    }
}

void enArrayQueue (struct ArrayQueue *q, __u32 k)
{

    q->queue[q->tail] = k;
    q->tail = (q->tail + 1) % q->size;
    q->tail_next = q->tail;

}


void enArrayQueueWithBitmap (struct ArrayQueue *q, __u32 k)
{

    q->queue[q->tail] = k;
    setBit(q->q_bitmap, k);
    q->tail = q->tail_next;
    q->tail++;
    q->tail_next++;

}


void enArrayQueueAtomic (struct ArrayQueue *q, __u32 k)
{

    __u32 local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->queue[local_q_tail] = k;

}


void enArrayQueueWithBitmapAtomic (struct ArrayQueue *q, __u32 k)
{

    __u32 local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->queue[local_q_tail] = k;
    setBitAtomic(q->q_bitmap, k);

}


void enArrayQueueDelayed (struct ArrayQueue *q, __u32 k)
{

    q->queue[q->tail_next] = k;
    q->tail_next++;

}

void enArrayQueueDelayedWithBitmap (struct ArrayQueue *q, __u32 k)
{

    q->queue[q->tail_next] = k;
    setBit(q->q_bitmap_next, k);
    q->tail_next++;

}

void enArrayQueueDelayedWithBitmapAtomic (struct ArrayQueue *q, __u32 k)
{

    __u32 local_q_tail_next = __sync_fetch_and_add(&q->tail_next, 1);
    setBitAtomic(q->q_bitmap, k);
    q->queue[local_q_tail_next] = k;

}


void slideWindowArrayQueue (struct ArrayQueue *q)
{

    q->head = q->tail;
    q->tail = q->tail_next;

}

void slideWindowArrayQueueBitmap (struct ArrayQueue *q)
{

    q->head = q->tail;
    q->tail = q->tail_next;
    swapBitmaps(&q->q_bitmap, &q->q_bitmap_next);
    clearBitmap(q->q_bitmap_next);

}

__u32 deArrayQueue(struct ArrayQueue *q)
{

    __u32 k = q->queue[q->head];
    clearBit(q->q_bitmap, k);
    q->head = (q->head + 1) % q->size;

    return k;

}


__u32 frontArrayQueue (struct ArrayQueue *q)
{

    __u32 k = q->queue[q->head];

    return k;

}

__u8 isEmptyArrayQueueCurr (struct ArrayQueue *q)
{

    if((q->tail > q->head))
        return 0;
    else
        return 1;

}

__u8 isEmptyArrayQueue (struct ArrayQueue *q)
{

    if(!isEmptyArrayQueueCurr(q) || !isEmptyArrayQueueNext(q))
        return 0;
    else
        return 1;

}

__u8 isEmptyArrayQueueNext (struct ArrayQueue *q)
{

    if((q->tail_next > q->head))
        return 0;
    else
        return 1;

}

__u8  isEnArrayQueued   (struct ArrayQueue *q, __u32 k)
{


    return getBit(q->q_bitmap, k);

}

__u8  isEnArrayQueuedNext   (struct ArrayQueue *q, __u32 k)
{


    return getBit(q->q_bitmap_next, k);

}

__u32 sizeArrayQueueCurr(struct ArrayQueue *q)
{

    return q->tail - q->head;

}

__u32 sizeArrayQueueNext(struct ArrayQueue *q)
{

    return q->tail_next - q->tail;
}


__u32 sizeArrayQueue(struct ArrayQueue *q)
{

    return q->tail_next - q->head;

}

void flushArrayQueueToShared(struct ArrayQueue *local_q, struct ArrayQueue *shared_q)
{

    __u32 shared_q_tail_next = __sync_fetch_and_add(&shared_q->tail_next, local_q->tail);
    __u32 local_q_size = local_q->tail - local_q->head;

    memcpy(&shared_q->queue[shared_q_tail_next], &local_q->queue[local_q->head], local_q_size * (sizeof(__u32)));

    local_q->head = 0;
    local_q->tail = 0;
    local_q->tail_next = 0;

}



void arrayQueueGenerateBitmap(struct ArrayQueue *q)
{

    __u32 v;
    __u32 i;

    #pragma omp parallel for
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->queue[i];
        setBitAtomic(q->q_bitmap, v);
    }

}


void arrayQueueToBitmap(struct ArrayQueue *q, struct Bitmap *b)
{

    __u32 v;
    __u32 i;

    #pragma omp parallel for default(none) shared(q,b) private(v,i)
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->queue[i];
        setBitAtomic(b, v);
    }

    // b->numSetBits = q->q_bitmap->numSetBits;
    q->head = q->tail;
    q->tail_next = q->tail;


}

void bitmapToArrayQueue(struct Bitmap *b, struct ArrayQueue *q, struct ArrayQueue **localFrontierQueues)
{

    #pragma omp parallel default(none) shared(b,localFrontierQueues,q)
    {
        __u32 i;

        __u32 t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];

        #pragma omp for
        for(i = 0 ; i < (b->size); i++)
        {
            if(getBit(b, i))
            {
                localFrontierQueue->queue[localFrontierQueue->tail] = i;
                localFrontierQueue->tail++;
            }

        }

        flushArrayQueueToShared(localFrontierQueue, q);

    }

    slideWindowArrayQueue(q);

}

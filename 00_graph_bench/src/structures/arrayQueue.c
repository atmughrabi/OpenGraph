// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
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
#include <stdint.h>
#include <omp.h>

#include "myMalloc.h"
#include "arrayQueue.h"
#include "bitmap.h"

struct ArrayQueue *newArrayQueue(uint32_t size)
{

    struct ArrayQueue *arrayQueue = (struct ArrayQueue *) my_malloc( sizeof(struct ArrayQueue));


    arrayQueue->head = 0;
    arrayQueue->tail = 0;
    arrayQueue->tail_next = 0;
    arrayQueue->size = size;

    arrayQueue->queue = (uint32_t *) my_malloc(size * sizeof(uint32_t));

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

void enArrayQueue (struct ArrayQueue *q, uint32_t k)
{

    q->queue[q->tail] = k;
    q->tail = (q->tail + 1) % q->size;
    q->tail_next = q->tail;

}


void enArrayQueueWithBitmap (struct ArrayQueue *q, uint32_t k)
{

    q->queue[q->tail] = k;
    setBit(q->q_bitmap, k);
    q->tail = q->tail_next;
    q->tail++;
    q->tail_next++;

}


void enArrayQueueAtomic (struct ArrayQueue *q, uint32_t k)
{

    uint32_t local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->queue[local_q_tail] = k;

}


void enArrayQueueWithBitmapAtomic (struct ArrayQueue *q, uint32_t k)
{

    uint32_t local_q_tail = __sync_fetch_and_add(&q->tail, 1);
    q->queue[local_q_tail] = k;
    setBitAtomic(q->q_bitmap, k);

}


void enArrayQueueDelayed (struct ArrayQueue *q, uint32_t k)
{

    q->queue[q->tail_next] = k;
    q->tail_next++;

}

void enArrayQueueDelayedWithBitmap (struct ArrayQueue *q, uint32_t k)
{

    q->queue[q->tail_next] = k;
    setBit(q->q_bitmap_next, k);
    q->tail_next++;

}

void enArrayQueueDelayedWithBitmapAtomic (struct ArrayQueue *q, uint32_t k)
{

    uint32_t local_q_tail_next = __sync_fetch_and_add(&q->tail_next, 1);
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

uint32_t deArrayQueue(struct ArrayQueue *q)
{

    uint32_t k = q->queue[q->head];
    clearBit(q->q_bitmap, k);
    q->head = (q->head + 1) % q->size;

    return k;

}


uint32_t frontArrayQueue (struct ArrayQueue *q)
{

    uint32_t k = q->queue[q->head];

    return k;

}

uint8_t isEmptyArrayQueueCurr (struct ArrayQueue *q)
{

    if((q->tail > q->head))
        return 0;
    else
        return 1;

}

uint8_t isEmptyArrayQueue (struct ArrayQueue *q)
{

    if(!isEmptyArrayQueueCurr(q) || !isEmptyArrayQueueNext(q))
        return 0;
    else
        return 1;

}

uint8_t isEmptyArrayQueueNext (struct ArrayQueue *q)
{

    if((q->tail_next > q->head))
        return 0;
    else
        return 1;

}

uint8_t  isEnArrayQueued   (struct ArrayQueue *q, uint32_t k)
{


    return getBit(q->q_bitmap, k);

}

uint8_t  isEnArrayQueuedNext   (struct ArrayQueue *q, uint32_t k)
{


    return getBit(q->q_bitmap_next, k);

}

uint32_t sizeArrayQueueCurr(struct ArrayQueue *q)
{

    return q->tail - q->head;

}

uint32_t sizeArrayQueueNext(struct ArrayQueue *q)
{

    return q->tail_next - q->tail;
}


uint32_t sizeArrayQueue(struct ArrayQueue *q)
{

    return q->tail_next - q->head;

}

void flushArrayQueueToShared(struct ArrayQueue *local_q, struct ArrayQueue *shared_q)
{

    uint32_t shared_q_tail_next = __sync_fetch_and_add(&shared_q->tail_next, local_q->tail);
    uint32_t local_q_size = local_q->tail - local_q->head;

    memcpy(&shared_q->queue[shared_q_tail_next], &local_q->queue[local_q->head], local_q_size * (sizeof(uint32_t)));

    local_q->head = 0;
    local_q->tail = 0;
    local_q->tail_next = 0;

}



void arrayQueueGenerateBitmap(struct ArrayQueue *q)
{

    uint32_t v;
    uint32_t i;

    #pragma omp parallel for
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->queue[i];
        setBitAtomic(q->q_bitmap, v);
    }

}


void arrayQueueToBitmap(struct ArrayQueue *q, struct Bitmap *b)
{

    uint32_t v;
    uint32_t i;

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
        uint32_t i;

        uint32_t t_id = omp_get_thread_num();
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

void arrayQueueToBitmapDualOrder(struct ArrayQueue *q, struct Bitmap *b, uint32_t *labels)
{

    uint32_t v;
    uint32_t i;
    uint32_t inv_u;
    uint32_t num_threads_max = omp_get_max_threads();

    #pragma omp parallel for default(none) shared(q,b,labels) private(v,i,inv_u) num_threads(num_threads_max)
    for(i = q->head ; i < q->tail; i++)
    {
        v = q->queue[i];
        inv_u = labels[v];
        setBitAtomic(b, inv_u);
    }

    // b->numSetBits = q->q_bitmap->numSetBits;
    q->head = q->tail;
    q->tail_next = q->tail;


}


void bitmapToArrayQueueDualOrder(struct Bitmap *b, struct ArrayQueue *q, struct ArrayQueue **localFrontierQueues, uint32_t *labels)
{

   
    #pragma omp parallel default(none) shared(b,localFrontierQueues,q,labels)
    {
        uint32_t i;
        uint32_t inv_v;
        uint32_t t_id = omp_get_thread_num();
        struct ArrayQueue *localFrontierQueue = localFrontierQueues[t_id];

        #pragma omp for
        for(i = 0 ; i < (b->size); i++)
        {
            if(getBit(b, i))
            {
                inv_v = labels[i];
                // uint32_t shared_q_tail_next = __sync_fetch_and_add(&localFrontierQueue->tail, 1);
                localFrontierQueue->queue[localFrontierQueue->tail] = inv_v;

                // #pragma omp atomic
                localFrontierQueue->tail++;
            }

        }

        flushArrayQueueToShared(localFrontierQueue, q);

    }

    slideWindowArrayQueue(q);

}

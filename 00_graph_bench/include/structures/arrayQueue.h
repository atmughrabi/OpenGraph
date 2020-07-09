#ifndef ARRAYQUEUE_H
#define ARRAYQUEUE_H

#include <stdint.h>
#include "bitmap.h"

struct  ArrayQueue
{
    uint32_t head;
    uint32_t tail;
    uint32_t tail_next;
    uint32_t size;
    uint32_t *queue;
    struct Bitmap *q_bitmap;
    struct Bitmap *q_bitmap_next;
};


struct ArrayQueue *newArrayQueue 	(uint32_t size);
void 	freeArrayQueue				(struct ArrayQueue *q);
void resetArrayQueue(struct ArrayQueue *q);
void softResetArrayQueue(struct ArrayQueue *q);

void enArrayQueue 					(struct ArrayQueue *q, uint32_t k);
void enArrayQueueWithBitmap 		(struct ArrayQueue *q, uint32_t k);
void enArrayQueueAtomic 			(struct ArrayQueue *q, uint32_t k);
void enArrayQueueWithBitmapAtomic 	(struct ArrayQueue *q, uint32_t k);


uint32_t 	deArrayQueue	(struct ArrayQueue *q);
uint32_t 	frontArrayQueue (struct ArrayQueue *q);
uint8_t  isEmptyArrayQueue (struct ArrayQueue *q);
uint8_t  isEnArrayQueued 	(struct ArrayQueue *q, uint32_t k);

void enArrayQueueDelayed 	(struct ArrayQueue *q, uint32_t k);
void enArrayQueueDelayedWithBitmapAtomic (struct ArrayQueue *q, uint32_t k);
void enArrayQueueDelayedWithBitmap (struct ArrayQueue *q, uint32_t k);

void slideWindowArrayQueue (struct ArrayQueue *q);
void slideWindowArrayQueueBitmap (struct ArrayQueue *q);

uint8_t isEmptyArrayQueueNext (struct ArrayQueue *q);
uint8_t isEmptyArrayQueueCurr (struct ArrayQueue *q);

uint32_t sizeArrayQueueCurr(struct ArrayQueue *q);
uint32_t sizeArrayQueueNext(struct ArrayQueue *q);

uint32_t sizeArrayQueue(struct ArrayQueue *q);
uint8_t  isEnArrayQueuedNext 	(struct ArrayQueue *q, uint32_t k);

void arrayQueueGenerateBitmap(struct ArrayQueue *q);
void flushArrayQueueToShared(struct ArrayQueue *local_q, struct ArrayQueue *shared_q);
void arrayQueueToBitmap(struct ArrayQueue *q, struct Bitmap *b);
void bitmapToArrayQueue(struct Bitmap *b, struct ArrayQueue *q, struct ArrayQueue **localFrontierQueues);

#endif



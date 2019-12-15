#ifndef ARRAYQUEUE_H
#define ARRAYQUEUE_H

#include <linux/types.h>
#include "bitmap.h"

struct  ArrayQueue
{
    __u32 head;
    __u32 tail;
    __u32 tail_next;
    __u32 size;
    __u32 *queue;
    struct Bitmap *q_bitmap;
    struct Bitmap *q_bitmap_next;

};


struct ArrayQueue *newArrayQueue 	(__u32 size);
void 	freeArrayQueue				(struct ArrayQueue *q);
void resetArrayQueue(struct ArrayQueue *q);
void softResetArrayQueue(struct ArrayQueue *q);

void enArrayQueue 					(struct ArrayQueue *q, __u32 k);
void enArrayQueueWithBitmap 		(struct ArrayQueue *q, __u32 k);
void enArrayQueueAtomic 			(struct ArrayQueue *q, __u32 k);
void enArrayQueueWithBitmapAtomic 	(struct ArrayQueue *q, __u32 k);


__u32 	deArrayQueue	(struct ArrayQueue *q);
__u32 	frontArrayQueue (struct ArrayQueue *q);
__u8  isEmptyArrayQueue (struct ArrayQueue *q);
__u8  isEnArrayQueued 	(struct ArrayQueue *q, __u32 k);

void enArrayQueueDelayed 	(struct ArrayQueue *q, __u32 k);
void enArrayQueueDelayedWithBitmapAtomic (struct ArrayQueue *q, __u32 k);
void enArrayQueueDelayedWithBitmap (struct ArrayQueue *q, __u32 k);

void slideWindowArrayQueue (struct ArrayQueue *q);
void slideWindowArrayQueueBitmap (struct ArrayQueue *q);

__u8 isEmptyArrayQueueNext (struct ArrayQueue *q);
__u8 isEmptyArrayQueueCurr (struct ArrayQueue *q);

__u32 sizeArrayQueueCurr(struct ArrayQueue *q);
__u32 sizeArrayQueueNext(struct ArrayQueue *q);

__u32 sizeArrayQueue(struct ArrayQueue *q);
__u8  isEnArrayQueuedNext 	(struct ArrayQueue *q, __u32 k);

void arrayQueueGenerateBitmap(struct ArrayQueue *q);
void flushArrayQueueToShared(struct ArrayQueue *local_q, struct ArrayQueue *shared_q);
void arrayQueueToBitmap(struct ArrayQueue *q, struct Bitmap *b);
void bitmapToArrayQueue(struct Bitmap *b, struct ArrayQueue *q, struct ArrayQueue **localFrontierQueues);

#endif



#ifndef ARRAYSTACK_H
#define ARRAYSTACK_H

#include <linux/types.h>
#include "bitmap.h"

struct  ArrayStack
{
    __u32 head;
    __u32 tail;
    __u32 tail_next;
    __u32 size;
    __u32 *Stack;
    struct Bitmap *q_bitmap;
    struct Bitmap *q_bitmap_next;

};


struct ArrayStack *newArrayStack 	(__u32 size);
void 	freeArrayStack				(struct ArrayStack *q);
void resetArrayStack(struct ArrayStack *q);

void pushArrayStack 					(struct ArrayStack *q, __u32 k);
void pushArrayStackWithBitmap 		(struct ArrayStack *q, __u32 k);
void pushArrayStackAtomic 			(struct ArrayStack *q, __u32 k);
void pushArrayStackWithBitmapAtomic 	(struct ArrayStack *q, __u32 k);


__u32 	popArrayStack	(struct ArrayStack *q);
__u32 	frontArrayStack (struct ArrayStack *q);
__u8  isEmptyArrayStack (struct ArrayStack *q);
__u8  ispushArrayStack 	(struct ArrayStack *q, __u32 k);

void pushArrayStackDelayed 	(struct ArrayStack *q, __u32 k);
void pushArrayStackDelayedWithBitmapAtomic (struct ArrayStack *q, __u32 k);
void pushArrayStackDelayedWithBitmap (struct ArrayStack *q, __u32 k);

void slideWindowArrayStack (struct ArrayStack *q);
void slideWindowArrayStackBitmap (struct ArrayStack *q);

__u8 isEmptyArrayStackNext (struct ArrayStack *q);
__u8 isEmptyArrayStackCurr (struct ArrayStack *q);

__u32 sizeArrayStackCurr(struct ArrayStack *q);
__u32 sizeArrayStackNext(struct ArrayStack *q);

__u32 sizeArrayStack(struct ArrayStack *q);
__u8  ispushArrayStackNext 	(struct ArrayStack *q, __u32 k);

void arrayStackGenerateBitmap(struct ArrayStack *q);
void flushArrayStackToShared(struct ArrayStack *local_q, struct ArrayStack *shared_q);
void arrayStackToBitmap(struct ArrayStack *q, struct Bitmap *b);
void bitmapToArrayStack(struct Bitmap *b, struct ArrayStack *q, struct ArrayStack **localFrontierStacks);

#endif



#ifndef ARRAYSTACK_H
#define ARRAYSTACK_H

#include <stdint.h>
#include "bitmap.h"

struct  ArrayStack
{
    uint32_t head;
    uint32_t tail;
    uint32_t tail_next;
    uint32_t size;
    uint32_t *Stack;
    struct Bitmap *q_bitmap;
    struct Bitmap *q_bitmap_next;

};


struct ArrayStack *newArrayStack 	(uint32_t size);
void 	freeArrayStack				(struct ArrayStack *q);
void resetArrayStack(struct ArrayStack *q);

void pushArrayStack 					(struct ArrayStack *q, uint32_t k);
void pushArrayStackWithBitmap 		(struct ArrayStack *q, uint32_t k);
void pushArrayStackAtomic 			(struct ArrayStack *q, uint32_t k);
void pushArrayStackWithBitmapAtomic 	(struct ArrayStack *q, uint32_t k);


uint32_t 	popArrayStack	(struct ArrayStack *q);
uint32_t 	frontArrayStack (struct ArrayStack *q);
uint8_t  isEmptyArrayStack (struct ArrayStack *q);
uint8_t  ispushArrayStack 	(struct ArrayStack *q, uint32_t k);

void pushArrayStackDelayed 	(struct ArrayStack *q, uint32_t k);
void pushArrayStackDelayedWithBitmapAtomic (struct ArrayStack *q, uint32_t k);
void pushArrayStackDelayedWithBitmap (struct ArrayStack *q, uint32_t k);

void slideWindowArrayStack (struct ArrayStack *q);
void slideWindowArrayStackBitmap (struct ArrayStack *q);

uint8_t isEmptyArrayStackNext (struct ArrayStack *q);
uint8_t isEmptyArrayStackCurr (struct ArrayStack *q);

uint32_t sizeArrayStackCurr(struct ArrayStack *q);
uint32_t sizeArrayStackNext(struct ArrayStack *q);

uint32_t sizeArrayStack(struct ArrayStack *q);
uint8_t  ispushArrayStackNext 	(struct ArrayStack *q, uint32_t k);

void arrayStackGenerateBitmap(struct ArrayStack *q);
void flushArrayStackToShared(struct ArrayStack *local_q, struct ArrayStack *shared_q);
void arrayStackToBitmap(struct ArrayStack *q, struct Bitmap *b);
void bitmapToArrayStack(struct Bitmap *b, struct ArrayStack *q, struct ArrayStack **localFrontierStacks);

#endif



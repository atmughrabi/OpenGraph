#ifndef BITMAP_H
#define BITMAP_H

#include <stdint.h>
#include "graphConfig.h"

#define kBitsPerWord  32
#define word_offset(n)  (n / kBitsPerWord)
#define bit_offset(n)  (n & (kBitsPerWord - 1))

#define ba_set(ptr, bit)    ((ptr)[(bit) >> 5] |= (uint32_t)(1 << ((bit) & kBitsPerWord)))
#define ba_clear(ptr, bit)  ((ptr)[(bit) >> 5] &= (uint32_t)(~(1 << ((bit) & kBitsPerWord))))
#define ba_get(ptr, bit)    ((ptr)[(bit) >> 5] & (uint32_t)(1 << ((bit) & kBitsPerWord)) ?  1 : 0 )
#define ba_setbit(ptr, bit, value) { if (value) { ba_set((ptr), (bit)) } else { ba_clear((ptr), (bit)); } }




struct  Bitmap
{
    uint32_t size;
    uint32_t numSetBits;
    uint32_t *bitarray;

};

struct Bitmap *newBitmapSet( uint32_t size);
struct Bitmap *newBitmap( uint32_t size);
void clearBitmap(struct Bitmap *bitmap);
void setBit(struct Bitmap *bitmap, uint32_t pos);
void setBitRange(struct Bitmap *bitmap, uint32_t start, uint32_t end);
void setBitAtomic(struct Bitmap *bitmap, uint32_t pos);
void setBitXOR(struct Bitmap *bitmap, uint32_t pos);
uint32_t getBit(struct Bitmap *bitmap, uint32_t pos);
void clearBit(struct Bitmap *bitmap, uint32_t pos);
void clearBitmap(struct Bitmap *bitmap);
void freeBitmap( struct Bitmap *bitmap);
struct Bitmap  *orBitmap(struct Bitmap *bitmap1, struct Bitmap *bitmap2);
struct Bitmap  *andBitmap(struct Bitmap *bitmap1, struct Bitmap *bitmap2);
uint32_t getNumOfSetBits (struct Bitmap *bitmap);
void swapBitmaps (struct Bitmap **bitmap1, struct Bitmap **bitmap2);
void printSetBits (struct Bitmap *bitmap);

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
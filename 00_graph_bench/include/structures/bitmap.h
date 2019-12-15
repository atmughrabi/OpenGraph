#ifndef BITMAP_H
#define BITMAP_H

#include <linux/types.h>
#include "graphConfig.h"

#define kBitsPerWord  32
#define word_offset(n)  (n / kBitsPerWord)
#define bit_offset(n)  (n & (kBitsPerWord - 1))

#define ba_set(ptr, bit)    ((ptr)[(bit) >> 5] |= (__u32)(1 << ((bit) & kBitsPerWord)))
#define ba_clear(ptr, bit)  ((ptr)[(bit) >> 5] &= (__u32)(~(1 << ((bit) & kBitsPerWord))))
#define ba_get(ptr, bit)    ((ptr)[(bit) >> 5] & (__u32)(1 << ((bit) & kBitsPerWord)) ?  1 : 0 )
#define ba_setbit(ptr, bit, value) { if (value) { ba_set((ptr), (bit)) } else { ba_clear((ptr), (bit)); } }




struct  Bitmap
{
    __u32 size;
    __u32 numSetBits;
    __u32 *bitarray;

};

struct Bitmap *newBitmapSet( __u32 size);
struct Bitmap *newBitmap( __u32 size);
void clearBitmap(struct Bitmap *bitmap);
void setBit(struct Bitmap *bitmap, __u32 pos);
void setBitRange(struct Bitmap *bitmap, __u32 start, __u32 end);
void setBitAtomic(struct Bitmap *bitmap, __u32 pos);
void setBitXOR(struct Bitmap *bitmap, __u32 pos);
__u32 getBit(struct Bitmap *bitmap, __u32 pos);
void clearBit(struct Bitmap *bitmap, __u32 pos);
void clearBitmap(struct Bitmap *bitmap);
void freeBitmap( struct Bitmap *bitmap);
struct Bitmap  *orBitmap(struct Bitmap *bitmap1, struct Bitmap *bitmap2);
struct Bitmap  *andBitmap(struct Bitmap *bitmap1, struct Bitmap *bitmap2);
__u32 getNumOfSetBits (struct Bitmap *bitmap);
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
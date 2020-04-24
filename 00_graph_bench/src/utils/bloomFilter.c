// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : bloomFilter.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>

#include <stdint.h>
#include "bloomFilter.h"
#include "bitmap.h"
#include "hash.h"
#include "myMalloc.h"

#include "graphConfig.h"

struct BloomFilter *newBloomFilter(uint32_t size, uint32_t k)
{



    struct BloomFilter *bloomFilter = (struct BloomFilter *) my_malloc( sizeof(struct BloomFilter));


    bloomFilter->bloom = newBitmap(size);
    bloomFilter->size = ((size + kBitsPerWord - 1) / kBitsPerWord) * kBitsPerWord;
    bloomFilter->k = k;
    bloomFilter->partition = bloomFilter->size / bloomFilter->k;

    return bloomFilter;


}
void freeBloomFilter( struct BloomFilter *bloomFilter)
{

    if(bloomFilter)
    {
        freeBitmap(bloomFilter->bloom);
        free(bloomFilter);
    }


}
void clearBloomFilter( struct BloomFilter *bloomFilter)
{

    clearBitmap(bloomFilter->bloom);

}


void addToBloomFilter(struct BloomFilter *bloomFilter, uint32_t item)
{


    uint64_t z = magicHash64((uint64_t)item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;

    for (i = 0; i < bloomFilter->k; ++i)
    {
        uint64_t k = (h1 + i * h2) % bloomFilter->partition; // bit to set
        uint64_t j = k + (i * bloomFilter->partition);       // in parition 'i'
        setBit(bloomFilter->bloom, (uint32_t)j);
    }

    // bloomFilter->size++;

}
uint32_t findInBloomFilter(struct BloomFilter *bloomFilter, uint32_t item)
{


    // MitzenmacherKirsch optimization
    uint64_t z = magicHash64((uint64_t)item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;

    for (i = 0; i < bloomFilter->k; ++i)
    {
        uint64_t k = (h1 + i * h2) % bloomFilter->partition; // bit to set
        uint64_t j = k + (i * bloomFilter->partition);       // in parition 'i'
        if(!getBit(bloomFilter->bloom, (uint32_t)j))
            return 0;

    }




    return 1;


}


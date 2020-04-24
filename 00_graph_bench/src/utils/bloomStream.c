// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : bloomStream.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <stdint.h>
#include "bloomStream.h"
#include "bitmap.h"
#include "hash.h"

#include "myMalloc.h"
#include "graphConfig.h"

struct BloomStream *newBloomStream(uint32_t size, uint32_t k)
{

    uint32_t i;
    uint32_t alignedSize = ((size + kBitsPerWord - 1) / kBitsPerWord) * kBitsPerWord;
    // uint32_t nextPrimePartition = findNextPrime((alignedSize/k));


    // alignedSize = (((nextPrimePartition*k)+kBitsPerWord - 1)/kBitsPerWord)*kBitsPerWord;



    struct BloomStream *bloomStream = (struct BloomStream *) my_malloc( sizeof(struct BloomStream));
    bloomStream->counter = (uint32_t *) my_malloc(alignedSize * sizeof(uint32_t));
    bloomStream->counterHistory = (uint32_t *) my_malloc(alignedSize * sizeof(uint32_t));


    for(i = 0 ; i < alignedSize; i++)
    {
        bloomStream->counter[i] = 0;
        bloomStream->counterHistory[i] = 0;
    }

    bloomStream->bloom = newBitmap(size);
    bloomStream->bloomPrime = newBitmap(size);
    bloomStream->bloomHistory = newBitmap(size);
    bloomStream->lowestCounter = newBitmap(size);


    bloomStream->size = alignedSize;
    bloomStream->k = k;
    bloomStream->partition = bloomStream->size / bloomStream->k;
    bloomStream->membership = 0;
    bloomStream->temperature  = 0;

    bloomStream->threashold = 0;
    bloomStream->decayPeriod  = 0;
    bloomStream->numIO = 0;



    return bloomStream;

}


void freeBloomStream( struct BloomStream *bloomStream)
{

    if(bloomStream)
    {
        freeBitmap(bloomStream->bloom);
        freeBitmap(bloomStream->bloomPrime);
        freeBitmap(bloomStream->bloomHistory);
        free(bloomStream->counter);
        free(bloomStream->counterHistory);
        free(bloomStream);
    }


}
void clearBloomStream( struct BloomStream *bloomStream)
{

    clearBitmap(bloomStream->bloom);
    clearBitmap(bloomStream->bloomPrime);

}


void addToBloomStream(struct BloomStream *bloomStream, uint64_t item)
{


    // printf("add- %lx %lu \n", item,item);
    uint64_t z = magicHash64(item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;
    uint32_t minCount = UINT_MAX;
    uint32_t freqCount = 0;
    uint32_t index = 0;

    uint32_t found = findInBloomStream(bloomStream, item);


    if(!found)
    {
        for (i = 0; i < bloomStream->k; ++i)
        {
            uint64_t k = (h1 + i * h2) % bloomStream->partition; // bit to set
            uint64_t j = k + (i * bloomStream->partition);       // in parition 'i'
            setBitXOR(bloomStream->bloom, (uint32_t)j);
        }

    }
    else
    {
        for (i = 0; i < bloomStream->k; ++i)
        {
            uint64_t k = (h1 + i * h2) % bloomStream->partition; // bit to set
            uint64_t j = k + (i * bloomStream->partition);       // in parition 'i'
            setBitXOR(bloomStream->bloomPrime, (uint32_t)j);
            bloomStream->counter[(uint32_t)j]++;


            freqCount = bloomStream->counterHistory[(uint32_t)j];
            if(minCount > freqCount)
            {
                index = (uint32_t)j;
                minCount = freqCount;
            }


        }

        swapBitmaps(&bloomStream->bloomPrime, &bloomStream->bloom);
        setBit(bloomStream->lowestCounter, index);

    }

    // BloomStream->size++;

}

uint32_t findInBloomStream(struct BloomStream *bloomStream, uint64_t item)
{


    // MitzenmacherKirsch optimization
    uint64_t z = magicHash64(item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;
    uint32_t index = 0;
    uint32_t found = 0;

    bloomStream->membership = 0;
    bloomStream->temperature  = 0;

    uint64_t k = 0; // bit to set
    uint64_t j = 0;       // in parition 'i'
    uint32_t minCount = UINT_MAX;
    uint32_t freqCount = 0;


    for (i = 0; i < bloomStream->k; ++i)
    {
        k = (h1 + i * h2) % bloomStream->partition; // bit to set
        j = k + (i * bloomStream->partition);       // in parition 'i'

        if(getBit(bloomStream->bloom, j))
        {
            freqCount = bloomStream->counter[(uint32_t)j];
        }
        else
        {
            freqCount = 0;
        }

        if(minCount > freqCount)
        {
            index = (uint32_t)j;
            minCount = freqCount;
        }

    }

    found = getBit(bloomStream->bloomHistory, index) | getBit(bloomStream->bloom, index) | getBit(bloomStream->bloomPrime, index);

    if(found)
    {
        bloomStream->membership = 1;
        bloomStream->temperature = bloomStream->counterHistory[index];


        // printf("FOUND item : %u counter : %u \n", item, bloomStream->counterHistory[index]);
    }
    else
    {
        bloomStream->membership = 0;
        bloomStream->temperature = 0;


        // printf("NOT FOUND\n");
    }

    return bloomStream->membership;
}


void aggregateBloomFilterToHistory(struct BloomStream *bloomStream)
{

    uint32_t i;

    for(i = 0 ; i < bloomStream->size ; i++)
    {

        if(!(bloomStream->counterHistory[i] / 2))
        {
            clearBit(bloomStream->bloomHistory, i);
        }

        // if(getBit(bloomStream->lowestCounter, i)){
        if(getBit(bloomStream->bloom, i) | getBit(bloomStream->bloomPrime, i) | getBit(bloomStream->bloomHistory, i))
        {
            setBit(bloomStream->bloomHistory, i);
        }
        // }

        bloomStream->counterHistory[i] = (bloomStream->counterHistory[i] / 2) + bloomStream->counter[i];
        bloomStream->counter[i] = 0;
    }

    clearBitmap(bloomStream->bloom);
    clearBitmap(bloomStream->bloomPrime);


}

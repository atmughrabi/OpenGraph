// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : bloomMultiHash.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <stdint.h>
#include "bloomMultiHash.h"
#include "bitmap.h"
#include "hash.h"
#include <limits.h>
#include "myMalloc.h"

#include "graphConfig.h"

struct BloomMultiHash *newBloomMultiHash(uint32_t size, double error)
{

    // uint32_t n = ceil(m / (-k / log(1 - exp(log(error) / k))))
    // uint32_t error = pow(1 - exp(-k / (m / n)), k)
    uint32_t m = ceil((size * log(error)) / log(1 / pow(2, log(2))));
    // uint32_t k = round((m / n) * log(2));

    uint32_t i;
    uint32_t alignedSize = ((m + kBitsPerWord - 1) / kBitsPerWord) * kBitsPerWord;





    struct BloomMultiHash *bloomMultiHash = (struct BloomMultiHash *) my_malloc( sizeof(struct BloomMultiHash));
    bloomMultiHash->counter = (uint32_t *) my_malloc(alignedSize * sizeof(uint32_t));
    bloomMultiHash->recency = newBitmap(alignedSize);


    for(i = 0 ; i < alignedSize; i++)
    {
        bloomMultiHash->counter[i] = 0;
    }

    bloomMultiHash->size = alignedSize;

    bloomMultiHash->threashold = 7;
    bloomMultiHash->decayPeriod  = alignedSize;
    bloomMultiHash->numIO = 0;

    bloomMultiHash->error = error;

    double num = log(bloomMultiHash->error);
    double denom = 0.480453013918201; // ln(2)^2
    bloomMultiHash->bpe = -(num / denom);

    bloomMultiHash->k = (uint32_t)ceil(0.693147180559945 * bloomMultiHash->bpe);  // ln(2)
    bloomMultiHash->partition = bloomMultiHash->size / bloomMultiHash->k;

    printf("n: %u \n", size);
    printf("p: %.2f%% \n", error);
    printf("m: %u \n", m);
    printf("k: %u \n", bloomMultiHash->k );

    return bloomMultiHash;

}

void freeBloomMultiHash( struct BloomMultiHash *bloomMultiHash)
{
    if(bloomMultiHash)
    {
        freeBitmap(bloomMultiHash->recency);
        free(bloomMultiHash->counter);
        free(bloomMultiHash);
    }
}

void addToBloomMultiHash(struct BloomMultiHash *bloomMultiHash, uint64_t item)
{

    uint64_t z = magicHash64(item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;

    bloomMultiHash->numIO++;

    for (i = 0; i < bloomMultiHash->k; ++i)
    {
        uint64_t k = (h1 + i * h2) % bloomMultiHash->partition; // bit to set
        uint64_t j = k + (i * bloomMultiHash->partition);       // in parition 'i'

        // if(getBit(bloomMultiHash->recency, j))
        // {
            // bloomMultiHash->counter[(uint32_t)j] += 2;
            // printf("%u %u %u\n",j,bloomMultiHash->counter[(uint32_t)j], item );
        // }
        // else
        // {
            bloomMultiHash->counter[(uint32_t)j]++;
            // printf("%u %u %u\n",j,bloomMultiHash->counter[(uint32_t)j], item );
            // setBit(bloomMultiHash->recency, j);
        // }

    }

    if(bloomMultiHash->numIO > bloomMultiHash->decayPeriod)
    {
        decayBloomMultiHash(bloomMultiHash);
    }

}

uint32_t findInBloomMultiHash(struct BloomMultiHash *bloomMultiHash, uint64_t item)
{


    // MitzenmacherKirsch optimization
    uint64_t z = magicHash64(item);
    uint64_t h1 = z & 0xffffffff;
    uint64_t h2 = z >> 32;
    uint64_t i;

    uint64_t k = 0; // bit to set
    uint64_t j = 0;       // in parition 'i'
    uint32_t freqCount = 0;


    for (i = 0; i < bloomMultiHash->k; ++i)
    {
        k = (h1 + i * h2) % bloomMultiHash->partition; // bit to set
        j = k + (i * bloomMultiHash->partition);       // in parition 'i'


        freqCount = bloomMultiHash->counter[(uint32_t)j];

        if(freqCount < bloomMultiHash->threashold)
            return 0;

    }

    return 1;
}


void decayBloomMultiHash(struct BloomMultiHash *bloomMultiHash)
{


    uint64_t i;
    for(i = 0 ; i < bloomMultiHash->size; i++)
    {
        bloomMultiHash->counter[i] /= 2;
    }

    // clearBitmap(bloomMultiHash->recency);

    bloomMultiHash->numIO = 0;

}
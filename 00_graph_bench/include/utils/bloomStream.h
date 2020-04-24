#ifndef BLOOMSTREAM_H
#define BLOOMSTREAM_H


#include <stdint.h>
#include "bitmap.h"

struct BloomStream
{
    struct Bitmap *bloom;
    struct Bitmap *bloomPrime;
    struct Bitmap *bloomHistory;
    struct Bitmap *lowestCounter;

    uint32_t *counter;
    uint32_t *counterHistory;
    uint32_t size; // size of bloom filter
    uint32_t partition; // partition m/k as a prime number
    uint32_t k; // number of hash function


    //pass these variables after find in bloomfilter
    uint32_t membership;
    uint32_t temperature;


    //pass these variables after find in bloomfilter
    uint32_t threashold;
    uint32_t decayPeriod;
    uint32_t numIO;
};


struct BloomStream *newBloomStream(uint32_t size, uint32_t k);
void freeBloomStream( struct BloomStream *bloomStream);
void clearBloomStream( struct BloomStream *bloomStream);
void addToBloomStream(struct BloomStream *bloomStream, uint64_t item);
uint32_t findInBloomStream(struct BloomStream *bloomStream, uint64_t item);
void aggregateBloomFilterToHistory(struct BloomStream *bloomStream);



#endif
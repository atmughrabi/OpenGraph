#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H


#include <stdint.h>
#include "bitmap.h"

struct BloomFilter
{
    struct Bitmap *bloom;
    uint32_t size; // size of bloom filter
    uint32_t partition; // partition m/k as a prime number
    uint32_t k; // number of hash function
};

struct BloomFilter *newBloomFilter(uint32_t size, uint32_t k);
void freeBloomFilter( struct BloomFilter *bloomFilter);
void clearBloomFilter( struct BloomFilter *bloomFilter);
void addToBloomFilter(struct BloomFilter *bloomFilter, uint32_t item);
uint32_t findInBloomFilter(struct BloomFilter *bloomFilter, uint32_t item);


#endif
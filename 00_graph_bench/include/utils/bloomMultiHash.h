#ifndef BLOOMMULTIHASH_H
#define BLOOMMULTIHASH_H


#include <stdint.h>
#include "bitmap.h"

struct BloomMultiHash
{

    uint32_t *counter;
    struct Bitmap *recency;

    uint32_t size; // size of bloom filter
    uint32_t partition; // partition m/k as a prime number
    uint32_t k; // number of hash function


    //pass these variables after find in bloomfilter
    uint32_t threashold;
    uint32_t decayPeriod;
    uint32_t numIO;

    double bpe;
    double error;
};


struct BloomMultiHash *newBloomMultiHash(uint32_t size, double error);
void freeBloomMultiHash( struct BloomMultiHash *bloomMultiHash);
void addToBloomMultiHash(struct BloomMultiHash *bloomMultiHash, uint64_t item);
uint32_t findInBloomMultiHash(struct BloomMultiHash *bloomMultiHash, uint64_t item);
void decayBloomMultiHash(struct BloomMultiHash *bloomMultiHash);



#endif
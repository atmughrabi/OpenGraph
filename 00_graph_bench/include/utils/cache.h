#ifndef CACHE_H
#define CACHE_H

#include <stdint.h>

//CAPI PSL CACHE default CONFIGS
// #define BLOCKSIZE 128
// #define L1_SIZE 262144
// #define L1_ASSOC 8
// #define POLICY 0

// #define BLOCKSIZE 64
// #define L1_SIZE 32768
// #define L1_ASSOC 8

//GRASP default configs
#define BLOCKSIZE 64
#define L1_SIZE 1048576
#define L1_ASSOC 16
#define POLICY 2

// Policy TYPES
#define LRU_POLICY 0
#define GRASP_POLICY 1
#define LFU_POLICY 2

// Cache states Constants
#define INVALID 0
#define VALID 1
#define DIRTY 2

// LFU Policy Constants
#define FREQ_BITS 8
#define FREQ_MAX ((2 << FREQ_BITS) - 1)

// GRASP Policy Constants RRIP (re-refernece insertion prediction)
#define NUM_BITS_RRIP 3

#define DEFAULT_INSERT_RRPV ((2 << NUM_BITS_RRIP) - 1)
#define COLD_INSERT_RRPV DEFAULT_INSERT_RRPV
#define WARM_INSERT_RRPV DEFAULT_INSERT_RRPV - 1
#define HOT_INSERT_RRPV  1
#define HOT_HIT_RRPV  0
#define RRPV_INIT DEFAULT_INSERT_RRPV


struct CacheLine
{
    uint64_t tag;
    uint8_t Flags;// 0:invalid, 1:valid, 2:dirty
    uint64_t seq; // LRU   POLICY 0
    uint8_t RRPV; // GRASP POLICY 1
    uint8_t freq; // LFU   POLICY 2
};

struct Cache
{
    uint64_t size, lineSize, assoc, sets, log2Sets, log2Blk, tagMask, numLines, evictions;
    uint64_t reads, readMisses, readsPrefetch, readMissesPrefetch, writes, writeMisses, writeBacks;
    uint64_t access_counter;
    uint64_t policy;
    struct CacheLine **cacheLines;

    uint64_t currentCycle;
    uint64_t currentCycle_cache;
    uint64_t currentCycle_preftcher;
    //counters for graph performance on the cache
    uint32_t *verticesMiss;
    uint32_t *verticesHit;
    uint32_t *vertices_base_reuse;
    uint32_t *vertices_total_reuse;
    uint32_t *vertices_accesses;
    uint32_t  numVertices;
};


struct DoubleTaggedCache
{
    struct Cache *cold_cache;// psl_cache
    struct Cache *warm_cache;// warm tag
    struct Cache *hot_cache; // hot_cache
    struct Cache *ref_cache; // psl_cache
};

///cacheline helper functions
void initCacheLine(struct CacheLine *cacheLine);
uint64_t getTag(struct CacheLine *cacheLine);
uint8_t getFlags(struct CacheLine *cacheLine);
uint64_t getSeq(struct CacheLine *cacheLine);
uint8_t getFreq(struct CacheLine *cacheLine);
uint8_t getRRPV(struct CacheLine *cacheLine);
void setFreq(struct CacheLine *cacheLine, uint8_t freq);
void setRRPV(struct CacheLine *cacheLine, uint8_t RRPV);
void setSeq(struct CacheLine *cacheLine, uint64_t Seq);
void setFlags(struct CacheLine *cacheLine, uint8_t flags);
void setTag(struct CacheLine *cacheLine, uint64_t a);

void invalidate(struct CacheLine *cacheLine);
uint32_t isValid(struct CacheLine *cacheLine);


uint64_t calcTag(struct Cache *cache, uint64_t addr);
uint64_t calcIndex(struct Cache *cache, uint64_t addr);
uint64_t calcAddr4Tag(struct Cache *cache, uint64_t tag);


uint64_t getRM(struct Cache *cache);
uint64_t getWM(struct Cache *cache);
uint64_t getReads(struct Cache *cache);
uint64_t getWrites(struct Cache *cache);
uint64_t getWB(struct Cache *cache);
uint64_t getEVC(struct Cache *cache);
uint64_t getRMPrefetch(struct Cache *cache);
uint64_t getReadsPrefetch(struct Cache *cache);
void writeBack(struct Cache *cache, uint64_t addr);

void initCache(struct Cache *cache, int s, int a, int b, int p);
void Access(struct Cache *cache, uint64_t addr, unsigned char op, uint32_t node);
void Prefetch(struct Cache *cache, uint64_t addr, unsigned char op, uint32_t node);

uint32_t checkInCache(struct Cache *cache, uint64_t addr);

struct CacheLine *findLine(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineLRU(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineLFU(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineGRASP(struct Cache *cache, uint64_t addr);

struct CacheLine *fillLine(struct Cache *cache, uint64_t addr);
struct CacheLine *fillLineLRU(struct Cache *cache, uint64_t addr);
struct CacheLine *fillLineLFU(struct Cache *cache, uint64_t addr);
struct CacheLine *fillLineGRASP(struct Cache *cache, uint64_t addr);

struct CacheLine *findLineToReplace(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineToReplaceLRU(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineToReplaceLFU(struct Cache *cache, uint64_t addr);
struct CacheLine *findLineToReplaceGRASP(struct Cache *cache, uint64_t addr);

struct CacheLine *getLRU(struct Cache *cache, uint64_t addr);
struct CacheLine *getLFU(struct Cache *cache, uint64_t addr);
struct CacheLine *getGRASP(struct Cache *cache, uint64_t addr);

void updatePolicy(struct Cache *cache, struct CacheLine *line);
void updateLRU(struct Cache *cache, struct CacheLine *line);
void updateLFU(struct Cache *cache, struct CacheLine *line);
void updateGRASP(struct Cache *cache, struct CacheLine *line);


void printStats(struct Cache *cache);

struct Cache *newCache( uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices);
void freeCache(struct Cache *cache);


struct DoubleTaggedCache *newDoubleTaggedCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices);
void freeDoubleTaggedCache(struct DoubleTaggedCache *cache);
void online_cache_graph_stats(struct Cache *cache, uint32_t node);

#endif
// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : cache.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "quantization.h"
#include "myMalloc.h"
#include "cache.h"


void initCacheLine(struct CacheLine *cacheLine)
{
    cacheLine->tag = 0;
    cacheLine->Flags = 0;
}
uint64_t getAddr(struct CacheLine *cacheLine)
{
    return cacheLine->addr;
}
uint64_t getTag(struct CacheLine *cacheLine)
{
    return cacheLine->tag;
}
uint8_t getFlags(struct CacheLine *cacheLine)
{
    return cacheLine->Flags;
}
uint64_t getSeq(struct CacheLine *cacheLine)
{
    return cacheLine->seq;
}
uint8_t getFreq(struct CacheLine *cacheLine)
{
    return cacheLine->freq;
}
uint8_t getRRPV(struct CacheLine *cacheLine)
{
    return cacheLine->RRPV;
}
uint8_t getSRRPV(struct CacheLine *cacheLine)
{
    return cacheLine->SRRPV;
}
uint8_t getPIN(struct CacheLine *cacheLine)
{
    return cacheLine->PIN;
}
uint8_t getPLRU(struct CacheLine *cacheLine)
{
    return cacheLine->PLRU;
}
uint8_t getXPRRPV(struct CacheLine *cacheLine)
{
    return cacheLine->XPRRPV;
}



void setSeq(struct CacheLine *cacheLine, uint64_t Seq)
{
    cacheLine->seq = Seq;
}
void setFreq(struct CacheLine *cacheLine, uint8_t freq)
{
    cacheLine->freq = freq;
}
void setRRPV(struct CacheLine *cacheLine, uint8_t RRPV)
{
    cacheLine->RRPV = RRPV;
}
void setSRRPV(struct CacheLine *cacheLine, uint8_t SRRPV)
{
    cacheLine->SRRPV = SRRPV;
}
void setPIN(struct CacheLine *cacheLine, uint8_t PIN)
{
    cacheLine->PIN = PIN;
}
void setPLRU(struct CacheLine *cacheLine, uint8_t PLRU)
{
    cacheLine->PLRU = PLRU;
}
void setXPRRPV(struct CacheLine *cacheLine, uint8_t XPRRPV)
{
    cacheLine->XPRRPV = XPRRPV;
}



void setFlags(struct CacheLine *cacheLine, uint8_t flags)
{
    cacheLine->Flags = flags;
}
void setTag(struct CacheLine *cacheLine, uint64_t a)
{
    cacheLine->tag = a;
}
void setAddr(struct CacheLine *cacheLine, uint64_t addr)
{
    cacheLine->addr = addr;
}
void invalidate(struct CacheLine *cacheLine)
{
    cacheLine->idx    = 0;
    cacheLine->seq    = 0;
    cacheLine->tag    = 0;    //useful function
    cacheLine->Flags  = INVALID;
    cacheLine->RRPV   = RRPV_INIT;
    cacheLine->SRRPV  = SRRPV_INIT;
    cacheLine->PIN    = 0;
    cacheLine->PLRU   = 0;
    cacheLine->freq   = 0;
    cacheLine->XPRRPV = XPRRPV_INIT;
}
uint32_t isValid(struct CacheLine *cacheLine)
{
    return ((cacheLine->Flags) != INVALID);
}


//cache helper functions

uint64_t calcTag(struct Cache *cache, uint64_t addr)
{
    return (addr >> (cache->log2Blk) );
}
uint64_t calcIndex(struct Cache *cache, uint64_t addr)
{
    return ((addr >> cache->log2Blk) & cache->tagMask);
}
uint64_t calcAddr4Tag(struct Cache *cache, uint64_t tag)
{
    return (tag << (cache->log2Blk));
}

uint64_t getRM(struct Cache *cache)
{
    return cache->readMisses;
}
uint64_t getWM(struct Cache *cache)
{
    return cache->writeMisses;
}
uint64_t getReads(struct Cache *cache)
{
    return cache->reads;
}
uint64_t getWrites(struct Cache *cache)
{
    return cache->writes;
}
uint64_t getWB(struct Cache *cache)
{
    return cache->writeBacks;
}
uint64_t getEVC(struct Cache *cache)
{
    return cache->evictions;
}
uint64_t getRMPrefetch(struct Cache *cache)
{
    return cache->readMissesPrefetch;
}
uint64_t getReadsPrefetch(struct Cache *cache)
{
    return cache->readsPrefetch;
}
void writeBack(struct Cache *cache, uint64_t addr)
{
    cache->writeBacks++;
}

// ********************************************************************************************
// ***************               Cache comparison                                **************
// ********************************************************************************************

struct DoubleTaggedCache *newDoubleTaggedCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices, uint32_t policy, uint32_t numPropertyRegions)
{
    struct DoubleTaggedCache *cache = (struct DoubleTaggedCache *) my_malloc(sizeof(struct DoubleTaggedCache));

    cache->accel_graph = newAccelGraphCache(PSL_L1_SIZE, PSL_L1_ASSOC, PSL_BLOCKSIZE, num_vertices, PSL_POLICY, numPropertyRegions);
    cache->ref_cache   = newCache( l1_size, l1_assoc, blocksize, num_vertices, policy, numPropertyRegions);

    return cache;
}

void initDoubleTaggedCacheRegion(struct DoubleTaggedCache *cache, struct PropertyMetaData *propertyMetaData)
{
    initAccelGraphCacheRegion     (cache->accel_graph, propertyMetaData);
    initialzeCachePropertyRegions (cache->ref_cache, propertyMetaData, cache->ref_cache->size);
}

void freeDoubleTaggedCache(struct DoubleTaggedCache *cache)
{
    if(cache)
    {
        freeAccelGraphCache(cache->accel_graph);
        freeCache(cache->ref_cache);
        free(cache);
    }
}

// ********************************************************************************************
// ***************              AccelGraph Cache configuration                   **************
// ********************************************************************************************


struct AccelGraphCache *newAccelGraphCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices, uint32_t policy, uint32_t numPropertyRegions)
{
    struct AccelGraphCache *cache = (struct AccelGraphCache *) my_malloc(sizeof(struct AccelGraphCache));

    cache->cold_cache = newCache( l1_size, l1_assoc, blocksize, num_vertices, PSL_POLICY, numPropertyRegions);
    cache->warm_cache = newCache( l1_size, 16, 4, num_vertices, WARM_POLICY, numPropertyRegions);
    cache->hot_cache  = newCache( l1_size, 8, 4, num_vertices, HOT_POLICY, numPropertyRegions);

    return cache;
}

void initAccelGraphCacheRegion(struct AccelGraphCache *cache, struct PropertyMetaData *propertyMetaData)
{
    uint64_t size = cache->cold_cache->size + cache->warm_cache->size + cache->hot_cache->size;
    initialzeCachePropertyRegions (cache->cold_cache, propertyMetaData, size);
    initialzeCachePropertyRegions (cache->warm_cache, propertyMetaData, size);
    initialzeCachePropertyRegions (cache->hot_cache, propertyMetaData, size);
}

void freeAccelGraphCache(struct AccelGraphCache *cache)
{
    if(cache)
    {
        freeCache(cache->cold_cache);
        freeCache(cache->warm_cache);
        freeCache(cache->hot_cache);
        free(cache);
    }
}

// ********************************************************************************************
// ***************              general Cache functions                          **************
// ********************************************************************************************

struct Cache *newCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices, uint32_t policy, uint32_t numPropertyRegions)
{
    uint64_t i;
    uint64_t j;

    struct Cache *cache = ( struct Cache *) my_malloc(sizeof(struct Cache));
    initCache(cache, l1_size, l1_assoc, blocksize, policy);

    cache->num_buckets         = 11;
    cache->numPropertyRegions  = numPropertyRegions;
    cache->propertyRegions     = (struct PropertyRegion *)my_malloc(sizeof(struct PropertyRegion) * numPropertyRegions);

    cache->thresholds              = (uint64_t *)my_malloc(sizeof(uint64_t) * cache->num_buckets );
    cache->thresholds_count        = (uint64_t *)my_malloc(sizeof(uint64_t) * cache->num_buckets );
    cache->thresholds_totalDegrees = (uint64_t *)my_malloc(sizeof(uint64_t) * cache->num_buckets );
    cache->thresholds_avgDegrees   = (uint64_t *)my_malloc(sizeof(uint64_t) * cache->num_buckets );
    cache->regions_avgDegrees      = (uint64_t **)my_malloc(sizeof(uint64_t *) * numPropertyRegions);


    for(i = 0; i < cache->numPropertyRegions; i++)
    {
        cache->regions_avgDegrees[i] = (uint64_t *)my_malloc(sizeof(uint64_t) * (cache->num_buckets + 1) );
        for(j = 0; j < (cache->num_buckets + 1); j++)
        {
            cache->regions_avgDegrees[i][j] = 0;
        }
    }

    for(i = 0; i < cache->num_buckets ; i++)
    {
        cache->thresholds[i]               = 0;
        cache->thresholds_count[i]         = 0;
        cache->thresholds_totalDegrees[i]  = 0;
        cache->thresholds_avgDegrees[i]    = 0;
    }

    for(i = 0; i < numPropertyRegions; i++)
    {
        cache->propertyRegions[i].upper_bound = 0;
        cache->propertyRegions[i].hot_bound   = 0;
        cache->propertyRegions[i].warm_bound  = 0;
        cache->propertyRegions[i].lower_bound = 0;
    }

    cache->numVertices  = num_vertices;
    cache->verticesMiss = (uint64_t *)my_malloc(sizeof(uint64_t) * num_vertices);
    cache->verticesHit  = (uint64_t *)my_malloc(sizeof(uint64_t) * num_vertices);
    cache->vertices_base_reuse  = (uint64_t *)my_malloc(sizeof(uint64_t) * num_vertices);
    cache->vertices_total_reuse = (uint64_t *)my_malloc(sizeof(uint64_t) * num_vertices);
    cache->vertices_accesses    = (uint64_t *)my_malloc(sizeof(uint64_t) * num_vertices);

    for(i = 0; i < num_vertices; i++)
    {
        cache->verticesMiss[i] = 0;
        cache->verticesHit[i] = 0;
        cache->vertices_base_reuse[i] = 0;
        cache->vertices_total_reuse[i] = 0;
        cache->vertices_accesses[i] = 0;
    }

    return cache;
}

void freeCache(struct Cache *cache)
{
    uint64_t i;

    if(cache)
    {
        if(cache->propertyRegions)
            free(cache->propertyRegions);
        if(cache->verticesMiss)
            free(cache->verticesMiss);
        if(cache->verticesHit)
            free(cache->verticesHit);
        if(cache->vertices_base_reuse)
            free(cache->vertices_base_reuse);
        if(cache->vertices_total_reuse)
            free(cache->vertices_total_reuse);
        if(cache->vertices_accesses)
            free(cache->vertices_accesses);
        if(cache->thresholds)
            free(cache->thresholds);
        if(cache->thresholds_count)
            free(cache->thresholds_count);
        if(cache->thresholds_totalDegrees)
            free(cache->thresholds_totalDegrees);
        if(cache->thresholds_avgDegrees)
            free(cache->thresholds_avgDegrees);

        for(i = 0; i < cache->numPropertyRegions; i++)
        {
            if(cache->regions_avgDegrees[i])
                free(cache->regions_avgDegrees[i]);
        }

        if(cache->regions_avgDegrees)
            free(cache->regions_avgDegrees);

        if(cache->cacheLines)
        {
            for(i = 0; i < cache->sets; i++)
            {
                if(cache->cacheLines[i])
                    free(cache->cacheLines[i]);
            }
            free(cache->cacheLines);
        }
        free(cache);
    }
}

void initCache(struct Cache *cache, int s, int a, int b, int p)
{
    uint64_t i, j;
    cache->reads = cache->readMisses = cache->readsPrefetch = cache->readMissesPrefetch = cache->writes = cache->evictions = 0;
    cache->writeMisses = cache->writeBacks = cache->currentCycle_preftcher = cache->currentCycle_cache = cache->currentCycle = 0;

    cache->policy     = (uint32_t)(p);
    cache->size       = (uint64_t)(s);
    cache->lineSize   = (uint64_t)(b);
    cache->assoc      = (uint64_t)(a);
    cache->sets       = (uint64_t)((s / b) / a);
    cache->numLines   = (uint64_t)(s / b);
    cache->log2Sets   = (uint64_t)(log2(cache->sets));
    cache->log2Blk    = (uint64_t)(log2(b));
    cache->access_counter    = 0;

    //*******************//
    //initialize your counters here//
    //*******************//

    cache->tagMask = 0;
    for(i = 0; i < cache->log2Sets; i++)
    {
        cache->tagMask <<= 1;
        cache->tagMask |= 1;
    }

    /**create a two dimentional cache, sized as cache[sets][assoc]**/

    cache->cacheLines = (struct CacheLine **) my_malloc(cache->sets * sizeof(struct CacheLine *));
    for(i = 0; i < cache->sets; i++)
    {
        cache->cacheLines[i] = (struct CacheLine *) my_malloc((cache->assoc + 1) * sizeof(struct CacheLine));
        for(j = 0; j < cache->assoc + 1; j++)
        {
            invalidate(&(cache->cacheLines[i][j]));
        }
    }
}

void online_cache_graph_stats(struct Cache *cache, uint32_t node)
{
    uint32_t first_Access = 0;
    uint32_t i;
    uint32_t v;

    cache->vertices_accesses[node]++;
    cache->access_counter++;

    // if(cache->access_counter % 1000 == 0)
    // {
    //     #pragma omp parallel for
    //     for ( i = 0; i < cache->numVertices; ++i)
    //     {
    //         if(cache->vertices_base_reuse[i] != 0)
    //             cache->vertices_base_reuse[i] = cache->access_counter;
    //     }
    // }

    if(cache->vertices_base_reuse[node] == 0)
        first_Access = 1;

    if(first_Access)
    {
        cache->vertices_total_reuse[node] = 1;
    }
    else
    {
        cache->vertices_total_reuse[node] += (cache->access_counter - cache->vertices_base_reuse[node]);
        // printf("%s\n", );
    }

    cache->vertices_base_reuse[node]   = cache->access_counter;
    v = node;
    for (i = 0; i < 32; ++i)
    {
        cache->vertices_base_reuse[(node / (cache->lineSize / 4)) + i]   = cache->access_counter;
        v++;
        if(v >= cache->numVertices )
            break;
    }
}


uint32_t checkInCache(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *line = findLine(cache, addr);

    if(line == NULL)
        return 1;

    // updatePolicy(cache, findLine(cache, addr));
    // else
    // {
    //     struct CacheLine *lineLRU = getLRU(cache, addr);
    //     if(lineLRU == line)
    //         return 1;
    //     else
    //         return 0;
    // }

    return 0;
}

void Prefetch(struct Cache *cache, uint64_t addr, unsigned char op, uint32_t node)
{
    cache->currentCycle++;/*per cache global counter to maintain LRU order
      among cache ways, updated on every cache access*/
    cache->currentCycle_preftcher++;
    cache->readsPrefetch++;
    struct CacheLine *line = findLine(cache, addr);
    if(line == NULL)/*miss*/
    {
        cache->readMissesPrefetch++;
        fillLine(cache, addr);
    }
    else
    {
        /**since it's a hit, update LRU and update dirty flag**/
        updatePromotionPolicy(cache, line);
    }
}

// ********************************************************************************************
// ***************         AGING POLICIES                                        **************
// ********************************************************************************************

void updateAgingPolicy(struct Cache *cache)
{
    switch(cache->policy)
    {
    case LRU_POLICY:
        updateAgeLRU(cache);
        break;
    case GRASP_POLICY:
        updateAgeGRASP(cache);
        break;
    case LFU_POLICY:
        updateAgeLFU(cache);
        break;
    case GRASPXP_POLICY:
        updateAgeGRASPXP(cache);
        break;
    default :
        updateAgeLRU(cache);
    }
}

/*No aging for LRU*/
void updateAgeLRU(struct Cache *cache)
{

}

void updateAgeLFU(struct Cache *cache)
{
    uint64_t i, j;
    uint8_t freq = 0;
    for(i = 0; i < cache->sets; i++)
    {
        for(j = 0; j < cache->assoc; j++)
        {
            if(isValid(&(cache->cacheLines[i][j])))
            {
                freq = getFreq(&(cache->cacheLines[i][j]));
                if(freq > 0)
                    freq--;
                setFreq(&(cache->cacheLines[i][j]), freq);
            }
        }
    }
}

void updateAgeGRASPXP(struct Cache *cache)
{
    uint64_t i, j;
    uint8_t XPRRPV = 0;
    for(i = 0; i < cache->sets; i++)
    {
        for(j = 0; j < cache->assoc; j++)
        {
            if(isValid(&(cache->cacheLines[i][j])))
            {
                XPRRPV = getXPRRPV(&(cache->cacheLines[i][j]));
                if(XPRRPV > 0)
                    XPRRPV--;
                setXPRRPV(&(cache->cacheLines[i][j]), XPRRPV);
            }
        }
    }
}

/*No aging for GRASP*/
void updateAgeGRASP(struct Cache *cache)
{

}


// ********************************************************************************************
// ***************         INSERTION POLICIES                                    **************
// ********************************************************************************************

void updateInsertionPolicy(struct Cache *cache, struct CacheLine *line)
{
    switch(cache->policy)
    {
    case LRU_POLICY:
        updateInsertLRU(cache, line);
        break;
    case LFU_POLICY:
        updateInsertLFU(cache, line);
        break;
    case GRASP_POLICY:
        updateInsertGRASP(cache, line);
        break;
    case SRRIP_POLICY:
        updateInsertSRRIP(cache, line);
        break;
    case PIN_POLICY:
        updateInsertPIN(cache, line);
        break;
    case PLRU_POLICY:
        updateInsertPLRU(cache, line);
        break;
    case GRASPXP_POLICY:
        updateInsertGRASPXP(cache, line);
        break;
    default :
        updateInsertLRU(cache, line);
    }
}

/*upgrade LRU line to be MRU line*/
void updateInsertLRU(struct Cache *cache, struct CacheLine *line)
{
    setSeq(line, cache->currentCycle);
}

void updateInsertLFU(struct Cache *cache, struct CacheLine *line)
{
    uint8_t freq = 0;
    setFreq(line, freq);
}

void updateInsertGRASP(struct Cache *cache, struct CacheLine *line)
{
    if(inHotRegion(cache, line))
    {
        setRRPV(line, HOT_INSERT_RRPV);
    }
    else if (inWarmRegion(cache, line))
    {
        setRRPV(line, WARM_INSERT_RRPV);
    }
    else
    {
        setRRPV(line, DEFAULT_INSERT_RRPV);
    }
}

void updateInsertSRRIP(struct Cache *cache, struct CacheLine *line)
{
    uint8_t SRRPV = DEFAULT_INSERT_SRRPV;
    setSRRPV(line, SRRPV);
}

void updateInsertPIN(struct Cache *cache, struct CacheLine *line)
{
    if(inHotRegion(cache, line))
    {
        setSeq(line, cache->currentCycle);
        setPIN(line, 1);
    }
    else
    {
        setSeq(line, cache->currentCycle);
        setPIN(line, 0);
    }
}

void updateInsertPLRU(struct Cache *cache, struct CacheLine *line)
{
    uint64_t i, j, tag;
    uint8_t bit_sum = 0;
    uint8_t PLRU = 1;

    i      = calcIndex(cache, line->addr);
    tag = calcTag(cache, line->addr);


    for(j = 0; j < cache->assoc; j++)
    {
        if(getTag(&(cache->cacheLines[i][j])) != tag)
            bit_sum += getPLRU(&(cache->cacheLines[i][j]));
    }

    if(bit_sum == (cache->assoc - 1))
    {
        for(j = 0; j < (cache->assoc); j++)
        {
            setPLRU(&(cache->cacheLines[i][j]), 0);
        }
    }

    setPLRU(line, PLRU);
}

void updateInsertGRASPXP(struct Cache *cache, struct CacheLine *line)
{
    uint8_t XPRRPV = 0;
    XPRRPV = (uint8_t)getCacheRegionGRASPXP(cache, line);
    setXPRRPV(line, XPRRPV);
}

uint32_t inHotRegion(struct Cache *cache, struct CacheLine *line)
{
    uint32_t v;
    uint32_t result = 0;

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        if((line->addr >=  cache->propertyRegions[v].lower_bound) && (line->addr < cache->propertyRegions[v].hot_bound))
        {
            result = 1;
        }
    }

    return result;
}

uint32_t inWarmRegion(struct Cache *cache, struct CacheLine *line)
{
    uint32_t v;
    uint32_t result = 0;

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        if((line->addr >=  cache->propertyRegions[v].hot_bound) && (line->addr < cache->propertyRegions[v].warm_bound))
        {
            result = 1;
        }
    }
    return result;
}

uint32_t inHotRegionAddrGRASP(struct Cache *cache, uint64_t addr)
{
    uint32_t v;
    uint32_t result = 0;

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        if((addr >=  cache->propertyRegions[v].lower_bound) && (addr < cache->propertyRegions[v].hot_bound))
        {
            result = 1;
        }
    }

    return result;
}

uint32_t inWarmRegionAddrGRASP(struct Cache *cache, uint64_t addr)
{
    uint32_t v;
    uint32_t result = 0;

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        if((addr >=  cache->propertyRegions[v].hot_bound) && (addr < cache->propertyRegions[v].warm_bound))
        {
            result = 1;
        }
    }
    return result;
}


// ********************************************************************************************
// ***************         PROMOTION POLICIES                                    **************
// ********************************************************************************************

void updatePromotionPolicy(struct Cache *cache, struct CacheLine *line)
{
    switch(cache->policy)
    {
    case LRU_POLICY:
        updatePromoteLRU(cache, line);
        break;
    case LFU_POLICY:
        updatePromoteLFU(cache, line);
        break;
    case GRASP_POLICY:
        updatePromoteGRASP(cache, line);
        break;
    case SRRIP_POLICY:
        updatePromoteSRRIP(cache, line);
        break;
    case PIN_POLICY:
        updatePromotePIN(cache, line);
        break;
    case PLRU_POLICY:
        updatePromotePLRU(cache, line);
        break;
    case GRASPXP_POLICY:
        updatePromoteGRASPXP(cache, line);
        break;
    default :
        updatePromoteLRU(cache, line);
    }
}

/*upgrade LRU line to be MRU line*/
void updatePromoteLRU(struct Cache *cache, struct CacheLine *line)
{
    setSeq(line, cache->currentCycle);
}

void updatePromoteLFU(struct Cache *cache, struct CacheLine *line)
{
    uint8_t freq = getFreq(line);
    if(freq < FREQ_MAX)
        freq++;
    setFreq(line, freq);
}

void updatePromoteGRASP(struct Cache *cache, struct CacheLine *line)
{
    if(inHotRegion(cache, line))
    {
        setRRPV(line, HOT_HIT_RRPV);
    }
    else
    {
        uint8_t RRPV = getRRPV(line);
        if(RRPV > 0)
            RRPV--;
        setRRPV(line, RRPV);
    }
}

void updatePromoteSRRIP(struct Cache *cache, struct CacheLine *line)
{
    setSRRPV(line, HIT_SRRPV);
}

void updatePromotePIN(struct Cache *cache, struct CacheLine *line)
{
    setSeq(line, cache->currentCycle);
}

void updatePromotePLRU(struct Cache *cache, struct CacheLine *line)
{

    uint64_t i, j, tag;
    uint8_t bit_sum = 0;
    uint8_t PLRU = 1;

    i      = calcIndex(cache, line->addr);
    tag = calcTag(cache, line->addr);


    for(j = 0; j < cache->assoc; j++)
    {
        if(getTag(&(cache->cacheLines[i][j])) != tag)
            bit_sum += getPLRU(&(cache->cacheLines[i][j]));
    }

    if(bit_sum == (cache->assoc - 1))
    {
        for(j = 0; j < (cache->assoc); j++)
        {
            setPLRU(&(cache->cacheLines[i][j]), 0);
        }
    }

    setPLRU(line, PLRU);
}

void updatePromoteGRASPXP(struct Cache *cache, struct CacheLine *line)
{
    uint8_t XPRRPV = getXPRRPV(line);
    uint32_t v;
    uint32_t i;
    uint32_t avg;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        for ( i = 1; i < (cache->num_buckets + 1); ++i)
        {
            if((line->addr >=  cache->regions_avgDegrees[v][i - 1]) && (line->addr < cache->regions_avgDegrees[v][i]))
            {
                avg = cache->thresholds_totalDegrees[i] / cache->thresholds_count[i];
                if(XPRRPV > avg)
                    XPRRPV -= avg;
                else
                    XPRRPV = 0;
                break;
            }
        }
    }

    setXPRRPV(line, XPRRPV);
}

// ********************************************************************************************
// ***************         VICTIM EVICTION POLICIES                              **************
// ********************************************************************************************

struct CacheLine *getVictimPolicy(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *victim = NULL;

    switch(cache->policy)
    {
    case LRU_POLICY:
        victim = getVictimLRU(cache, addr);
        break;
    case LFU_POLICY:
        victim = getVictimLFU(cache, addr);
        break;
    case GRASP_POLICY:
        victim = getVictimGRASP(cache, addr);
        break;
    case SRRIP_POLICY:
        victim = getVictimSRRIP(cache, addr);
        break;
    case PIN_POLICY:
        victim = getVictimPIN(cache, addr);
        break;
    case PLRU_POLICY:
        victim = getVictimPLRU(cache, addr);
        break;
    case GRASPXP_POLICY:
        victim = getVictimGRASPXP(cache, addr);
        break;
    default :
        victim = getVictimLRU(cache, addr);
    }

    return victim;
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *getVictimLRU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = cache->currentCycle;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(getSeq(&(cache->cacheLines[i][j])) <= min)
        {
            victim = j;
            min = getSeq(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr;
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *getVictimLFU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = FREQ_MAX;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(getFreq(&(cache->cacheLines[i][j])) <= min)
        {
            victim = j;
            min = getFreq(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr;
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *getVictimGRASP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }

    victim = 0;
    min = getRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);

    // not in the GRASP paper optimizaiton
    if (min < DEFAULT_INSERT_RRPV)
    {
        int diff = DEFAULT_INSERT_RRPV - min;
        for(j = 0; j < cache->assoc; j++)
        {
            uint8_t RRPV = getRRPV(&(cache->cacheLines[i][j])) + diff;
            setRRPV(&(cache->cacheLines[i][j]), RRPV);
            assert(RRPV <= DEFAULT_INSERT_RRPV);
        }
    }

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr; // update victim with new address so we simulate hot/cold insertion
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *getVictimSRRIP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }

    // do
    // {
    //     for(j = 0; j < cache->assoc; j++)
    //     {
    //         if(getSRRPV(&(cache->cacheLines[i][j])) == SRRPV_INIT)
    //         {
    //             victim = j;
    //             min = getSRRPV(&(cache->cacheLines[i][j]));
    //             break;
    //         }
    //     }

    //     if(!min)
    //     {
    //         for(j = 0; j < cache->assoc; j++)
    //         {
    //             uint8_t SRRPV = getSRRPV(&(cache->cacheLines[i][j])) + 1;
    //             if(SRRPV <= SRRPV_INIT)
    //                 setSRRPV(&(cache->cacheLines[i][j]), SRRPV);
    //         }
    //     }

    // }
    // while(!min);

    victim = 0;
    min = getSRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getSRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getSRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);

    if (min < SRRPV_INIT)
    {
        int diff = SRRPV_INIT - min;
        for(j = 0; j < cache->assoc; j++)
        {
            uint8_t SRRPV = getSRRPV(&(cache->cacheLines[i][j])) + (diff);
            setSRRPV(&(cache->cacheLines[i][j]), SRRPV);
            assert(SRRPV <= SRRPV_INIT);
        }
    }

    assert(min != SRRPV_INIT || min != 0);
    assert(victim != cache->assoc);

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr;
    return &(cache->cacheLines[i][victim]);
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *getVictimPIN(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = cache->currentCycle;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(!getPIN(&(cache->cacheLines[i][j])) && (getSeq(&(cache->cacheLines[i][j])) <= min))
        {
            victim = j;
            min = getSeq(&(cache->cacheLines[i][j]));
        }
    }
    // assert(victim != cache->assoc);

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr; // update victim with new address so we simulate hot/cold insertion
    return &(cache->cacheLines[i][victim]);
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
uint8_t getVictimPINBypass(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, min;
    uint8_t bypass = 1;
    min    = cache->currentCycle;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return 0;
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(!getPIN(&(cache->cacheLines[i][j])) && (getSeq(&(cache->cacheLines[i][j])) <= min))
        {
            min = getSeq(&(cache->cacheLines[i][j]));
            bypass = 0;
        }
    }

    return bypass;
}


struct CacheLine *getVictimPLRU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim;
    victim = cache->assoc;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }

    for(j = 0; j < cache->assoc; j++)
    {
        if(!getPLRU(&(cache->cacheLines[i][j])))
        {
            victim = j;
            break;
        }
    }
    assert(victim != cache->assoc);

    cache->evictions++;

    cache->cacheLines[i][victim].addr = addr;
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *getVictimGRASPXP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }

    // do
    // {
    //     for(j = 0; j < cache->assoc; j++)
    //     {
    //         if(getXPRRPV(&(cache->cacheLines[i][j])) == 0)
    //         {
    //             victim = j;
    //             min = getXPRRPV(&(cache->cacheLines[i][j]));
    //             break;
    //         }
    //     }

    //     if(min)
    //     {
    //         for(j = 0; j < cache->assoc; j++)
    //         {
    //             uint8_t XPRRPV = getXPRRPV(&(cache->cacheLines[i][j])) - 1;
    //             if(XPRRPV > 0)
    //                 setXPRRPV(&(cache->cacheLines[i][j]), XPRRPV);
    //         }
    //     }

    // }
    // while(min);

    victim = 0;
    min = getXPRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getXPRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getXPRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);

    if (min < XPRRPV_INIT)
    {
        int diff = XPRRPV_INIT - min;
        for(j = 0; j < cache->assoc; j++)
        {
            uint8_t XPRRPV = getXPRRPV(&(cache->cacheLines[i][j])) + diff;
            setXPRRPV(&(cache->cacheLines[i][j]), XPRRPV);
            assert(XPRRPV <= XPRRPV_INIT);
        }
    }

    assert(min != XPRRPV_INIT || min != 0);
    assert(victim != cache->assoc);

    // victim = 0;
    // min = getXPRRPV(&(cache->cacheLines[i][0]));

    // for(j = 1; j < cache->assoc; j++)
    // {
    //     if(getXPRRPV(&(cache->cacheLines[i][j])) < min)
    //     {
    //         victim = j;
    //         min = getXPRRPV(&(cache->cacheLines[i][j]));
    //     }
    // }
    // assert(victim != cache->assoc);

    // if (min != XPRRPV_INIT)
    // {
    //     int diff = min;
    //     for(j = 0; j < cache->assoc; j++)
    //     {
    //         uint8_t XPRRPV = getXPRRPV(&(cache->cacheLines[i][j])) - diff;
    //         setXPRRPV(&(cache->cacheLines[i][j]), XPRRPV);
    //     }
    // }

    // min = getXPRRPV(&(cache->cacheLines[i][victim]));
    // assert(min == 0);
    // assert(victim != cache->assoc);

    cache->evictions++;
    cache->cacheLines[i][victim].addr = addr;
    return &(cache->cacheLines[i][victim]);
}

// ********************************************************************************************
// ***************         VICTIM PEEK POLICIES                                  **************
// ********************************************************************************************

struct CacheLine *peekVictimPolicy(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *victim = NULL;

    switch(cache->policy)
    {
    case LRU_POLICY:
        victim = peekVictimLRU(cache, addr);
        break;
    case LFU_POLICY:
        victim = peekVictimLFU(cache, addr);
        break;
    case GRASP_POLICY:
        victim = peekVictimGRASP(cache, addr);
        break;
    case SRRIP_POLICY:
        victim = peekVictimSRRIP(cache, addr);
        break;
    case PIN_POLICY:
        victim = peekVictimPIN(cache, addr);
        break;
    case PLRU_POLICY:
        victim = peekVictimPLRU(cache, addr);
        break;
    case GRASPXP_POLICY:
        victim = peekVictimGRASPXP(cache, addr);
        break;
    default :
        victim = peekVictimLRU(cache, addr);
    }

    return victim;
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *peekVictimLRU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = cache->currentCycle;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            cache->cacheLines[i][j].addr = addr;
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(getSeq(&(cache->cacheLines[i][j])) <= min)
        {
            victim = j;
            min = getSeq(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *peekVictimLFU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = FREQ_MAX;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(getFreq(&(cache->cacheLines[i][j])) <= min)
        {
            victim = j;
            min = getFreq(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *peekVictimGRASP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }

    victim = 0;
    min = getRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}

struct CacheLine *peekVictimSRRIP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }

    victim = 0;
    min = getSRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getSRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getSRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *peekVictimPIN(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = cache->currentCycle;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }
    for(j = 0; j < cache->assoc; j++)
    {
        if(!getPIN(&(cache->cacheLines[i][j])) && (getSeq(&(cache->cacheLines[i][j])) <= min))
        {
            victim = j;
            min = getSeq(&(cache->cacheLines[i][j]));
        }
    }
    // assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}


struct CacheLine *peekVictimPLRU(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim;
    victim = cache->assoc;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }

    for(j = 0; j < cache->assoc; j++)
    {
        if(!getPLRU(&(cache->cacheLines[i][j])))
        {
            victim = j;
            break;
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}


struct CacheLine *peekVictimGRASPXP(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, victim, min;

    victim = cache->assoc;
    min    = 0;
    i      = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
    {
        if(isValid(&(cache->cacheLines[i][j])) == 0)
        {
            return &(cache->cacheLines[i][j]);
        }
    }

    victim = 0;
    min = getXPRRPV(&(cache->cacheLines[i][0]));

    for(j = 1; j < cache->assoc; j++)
    {
        if(getXPRRPV(&(cache->cacheLines[i][j])) > min)
        {
            victim = j;
            min = getXPRRPV(&(cache->cacheLines[i][j]));
        }
    }
    assert(victim != cache->assoc);
    return &(cache->cacheLines[i][victim]);
}


// ********************************************************************************************
// ***************         Cacheline lookups                                     **************
// ********************************************************************************************

/*look up line*/
struct CacheLine *findLine(struct Cache *cache, uint64_t addr)
{
    uint64_t i, j, tag, pos;

    pos = cache->assoc;
    tag = calcTag(cache, addr);
    i   = calcIndex(cache, addr);

    for(j = 0; j < cache->assoc; j++)
        if(isValid((&cache->cacheLines[i][j])))
            if(getTag(&(cache->cacheLines[i][j])) == tag)
            {
                pos = j;
                break;
            }
    if(pos == cache->assoc)
        return NULL;
    else
    {
        return &(cache->cacheLines[i][pos]);
    }
}


/*find a victim*/
struct CacheLine *findLineToReplace(struct Cache *cache, uint64_t addr)
{
    struct CacheLine  *victim = getVictimPolicy(cache, addr);
    updateInsertionPolicy(cache, victim);

    return (victim);
}

/*allocate a new line*/
struct CacheLine *fillLine(struct Cache *cache, uint64_t addr)
{
    uint64_t tag;

    struct CacheLine *victim = findLineToReplace(cache, addr);
    assert(victim != 0);
    if(getFlags(victim) == DIRTY)
    {
        writeBack(cache, addr);
    }

    tag = calcTag(cache, addr);
    setTag(victim, tag);
    setAddr(victim, addr);
    setFlags(victim, VALID);

    /**note that this cache line has been already
       upgraded to MRU in the previous function (findLineToReplace)**/

    return victim;
}

void Access(struct Cache *cache, uint64_t addr, unsigned char op, uint32_t node)
{


    online_cache_graph_stats(cache, node);
    cache->currentCycle++;
    /*per cache global counter to maintain LRU order among cache ways, updated on every cache access*/

    cache->currentCycle_cache++;

    if(op == 'w')
    {
        cache->writes++;
    }
    else if(op == 'r')
    {
        cache->reads++;

    }

    struct CacheLine *line = findLine(cache, addr);
    if(line == NULL)/*miss*/
    {
        if(op == 'w')
        {
            cache->writeMisses++;
        }
        else if(op == 'r')
        {
            cache->readMisses++;
        }

        struct CacheLine *newline = NULL;

        if(cache->policy == PIN_POLICY)
        {
            if(!getVictimPINBypass(cache, addr))
            {
                newline = fillLine(cache, addr);
                newline->idx = node;
                if(op == 'w')
                    setFlags(newline, DIRTY);
            }
        }
        else
        {
            newline = fillLine(cache, addr);
            newline->idx = node;
            if(op == 'w')
                setFlags(newline, DIRTY);
        }


        cache->verticesMiss[node]++;
    }
    else
    {
        /**since it's a hit, update LRU and update dirty flag**/
        updatePromotionPolicy(cache, line);
        if(op == 'w')
            setFlags(line, DIRTY);

        cache->verticesHit[node]++;
    }
}

// ********************************************************************************************
// ***************               ACCElGraph Policy                               **************
// ********************************************************************************************

void AccessDoubleTaggedCacheFloat(struct DoubleTaggedCache *cache, uint64_t addr, unsigned char op, uint32_t node, float value)
{
    // AccessAccelGraphExpressFloat(cache->accel_graph, addr, op, node, value);

    AccessAccelGraphGRASP(cache->accel_graph, addr, op, node);
    Access(cache->ref_cache, addr, op, node);
}

void AccessAccelGraphGRASP(struct AccelGraphCache *accel_graph, uint64_t addr, unsigned char op, uint32_t node)
{
    struct CacheLine *victim = NULL;

    if(checkInCache(accel_graph->warm_cache, addr) && checkInCache(accel_graph->hot_cache, addr))
    {
        if(inHotRegionAddrGRASP(accel_graph->hot_cache, addr))
        {
            Access(accel_graph->cold_cache, addr, op, node);
            Access(accel_graph->hot_cache, addr, op, node);
            victim = peekVictimPolicy(accel_graph->hot_cache, addr);
            if(isValid(victim))
            {
                // Prefetch(accel_graph->warm_cache, victim->addr, 'r', victim_node);
                Access(accel_graph->warm_cache, victim->addr, 'd', victim->idx);
            }
        }
        else if(inWarmRegionAddrGRASP(accel_graph->warm_cache, addr))
        {
            Access(accel_graph->cold_cache, addr, op, node);
            Access(accel_graph->warm_cache, addr, op, node);
        }
        else
        {
            Access(accel_graph->cold_cache, addr, op, node);
        }
    }
    else  if(!checkInCache(accel_graph->warm_cache, addr) && checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->warm_cache, addr, op, node);
    }
    else  if(checkInCache(accel_graph->warm_cache, addr) && !checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->hot_cache, addr, op, node);
    }
    else  if(!checkInCache(accel_graph->warm_cache, addr) && !checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->hot_cache, addr, op, node);
    }
}


void AccessAccelGraphExpressFloat(struct AccelGraphCache *accel_graph, uint64_t addr, unsigned char op, uint32_t node, float value)
{
    struct CacheLine *victim = NULL;

    if(checkInCache(accel_graph->warm_cache, addr) && checkInCache(accel_graph->hot_cache, addr))
    {
        if(value <= 0.0015)
        {
            Access(accel_graph->cold_cache, addr, op, node);
            Access(accel_graph->hot_cache, addr, op, node);
            victim = peekVictimPolicy(accel_graph->hot_cache, addr);
            if(isValid(victim))
            {
                // Prefetch(accel_graph->warm_cache, victim->addr, 'r', victim_node);
                Access(accel_graph->warm_cache, victim->addr, 'd', victim->idx);
            }
        }
        else if(value > 0.0015 && value <= 0.015)
        {
            Access(accel_graph->cold_cache, addr, op, node);
            Access(accel_graph->warm_cache, addr, op, node);
        }
        else
        {
            Access(accel_graph->cold_cache, addr, op, node);
        }
    }
    else  if(!checkInCache(accel_graph->warm_cache, addr) && checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->warm_cache, addr, op, node);
    }
    else  if(checkInCache(accel_graph->warm_cache, addr) && !checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->hot_cache, addr, op, node);
    }
    else  if(!checkInCache(accel_graph->warm_cache, addr) && !checkInCache(accel_graph->hot_cache, addr))
    {
        Access(accel_graph->hot_cache, addr, op, node);
    }
}

// ********************************************************************************************
// ***************               GRASP-XP Policy                                 **************
// ********************************************************************************************

uint64_t getCacheRegionGRASPXP(struct Cache *cache, struct CacheLine *line)
{
    uint32_t v;
    uint32_t i;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        for ( i = 1; i < (cache->num_buckets + 1); ++i)
        {
            if((line->addr >=  cache->regions_avgDegrees[v][i - 1]) && (line->addr < cache->regions_avgDegrees[v][i]))
            {
                return cache->thresholds_avgDegrees[i - 1];
            }
        }
    }

    return 0;
}

uint64_t getCacheRegionAddrGRASPXP(struct Cache *cache, uint64_t addr)
{
    uint32_t v;
    uint32_t i;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        for ( i = 1; i < (cache->num_buckets + 1); ++i)
        {
            if((addr >=  cache->regions_avgDegrees[v][i - 1]) && (addr < cache->regions_avgDegrees[v][i]))
            {
                return cache->thresholds_avgDegrees[i - 1];
            }
        }
    }

    return 0;
}


void setCacheRegionDegreeAvg(struct Cache *cache)
{
    uint32_t v;
    uint32_t i;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        cache->regions_avgDegrees[v][0] = cache->propertyRegions[v].base_address;
        for ( i = 1; i < (cache->num_buckets + 1); ++i)
        {
            cache->regions_avgDegrees[v][i] = cache->regions_avgDegrees[v][i - 1] + (cache->thresholds_count[i - 1] * cache->propertyRegions[v].data_type_size);
        }
    }
}

void setDoubleTaggedCacheThresholdDegreeAvg(struct DoubleTaggedCache *cache, uint32_t  *degrees)
{
    setAccelGraphCacheThresholdDegreeAvg(cache->accel_graph, degrees);
    setCacheThresholdDegreeAvg(cache->ref_cache, degrees);
}


void setAccelGraphCacheThresholdDegreeAvg(struct AccelGraphCache *cache, uint32_t  *degrees)
{
    setCacheThresholdDegreeAvg(cache->cold_cache, degrees);
    setCacheThresholdDegreeAvg(cache->warm_cache, degrees);
    setCacheThresholdDegreeAvg(cache->hot_cache, degrees);
}

void setCacheThresholdDegreeAvg(struct Cache *cache, uint32_t  *degrees)
{
    uint32_t v;
    uint32_t i;
    uint64_t  avgDegrees = 0;
    uint64_t  totalDegrees = 0;
    float *thresholds_avgDegrees;
    thresholds_avgDegrees    = (float *) my_malloc(cache->num_buckets * sizeof(float));

    for (v = 0; v < cache->numVertices; ++v)
    {
        avgDegrees +=  degrees[v];
    }

    avgDegrees /= cache->numVertices;

    // START initialize thresholds
    if(avgDegrees <= 1)
        cache->thresholds[0] = 1;
    else
        cache->thresholds[0] = (avgDegrees / 2);
    for ( i = 1; i < (cache->num_buckets - 1); ++i)
    {
        cache->thresholds[i] = cache->thresholds[i - 1] * 2;
    }
    cache->thresholds[cache->num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds

    // collect stats perbucket
    for (v = 0; v < cache->numVertices; ++v)
    {
        for ( i = 0; i < cache->num_buckets; ++i)
        {
            if(degrees[v] <= cache->thresholds[i])
            {
                cache->thresholds_count[i] += 1;
                cache->thresholds_totalDegrees[i]  += degrees[v];
                break;
            }
        }
    }

    // collect stats perbucket
    for (v = 0; v < cache->numVertices; ++v)
    {
        totalDegrees += degrees[v];
    }


    for ( i = 0; i < cache->num_buckets; ++i)
    {
        if(cache->thresholds_count[i])
        {
            thresholds_avgDegrees[i] = XPRRPV_INIT * ((float)cache->thresholds_totalDegrees[i] / totalDegrees);
        }
        else
        {
            thresholds_avgDegrees[i] = 0;
        }
    }

    struct quant_params_8 rDivD_params;
    getMinMax_8(&rDivD_params, thresholds_avgDegrees, cache->num_buckets);
    rDivD_params.scale = GetScale_8(rDivD_params.min, rDivD_params.max);
    rDivD_params.zero = 0;

    for ( i = 0; i < cache->num_buckets; ++i)
    {
        cache->thresholds_avgDegrees[i]   =  quantize_8(thresholds_avgDegrees[i], rDivD_params.scale, rDivD_params.zero);
    }

    setCacheRegionDegreeAvg(cache);
    free(thresholds_avgDegrees);
}

// ********************************************************************************************
// ***************               GRASP Policy                                    **************
// ********************************************************************************************

// ********************************************************************************************
// ***************               GRASP Initializaiton                            **************
// ********************************************************************************************

void initialzeCachePropertyRegions (struct Cache *cache, struct PropertyMetaData *propertyMetaData, uint64_t size)
{
    uint32_t v;
    uint64_t total_properties_size = 0;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        total_properties_size += (propertyMetaData[v].size * propertyMetaData[v].data_type_size);
        cache->propertyRegions[v].base_address = propertyMetaData[v].base_address;
        cache->propertyRegions[v].size = propertyMetaData[v].size;
        cache->propertyRegions[v].data_type_size = propertyMetaData[v].data_type_size;
    }

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        // cache->propertyRegions[v].fraction    = 100; // classical vs ratio of array size in bytes
        cache->propertyRegions[v].fraction    = (uint64_t)(((uint64_t)(propertyMetaData[v].size * propertyMetaData[v].data_type_size) * 100) / total_properties_size );
        cache->propertyRegions[v].lower_bound = propertyMetaData[v].base_address;
        cache->propertyRegions[v].upper_bound = propertyMetaData[v].base_address + (uint64_t)(propertyMetaData[v].size * propertyMetaData[v].data_type_size);

        cache->propertyRegions[v].hot_bound = cache->propertyRegions[v].lower_bound + ((uint64_t)(size * cache->propertyRegions[v].fraction) / 100);
        if(cache->propertyRegions[v].hot_bound > cache->propertyRegions[v].upper_bound)
        {
            cache->propertyRegions[v].hot_bound = cache->propertyRegions[v].upper_bound;
        }

        cache->propertyRegions[v].warm_bound = cache->propertyRegions[v].hot_bound + ((uint64_t)(size * cache->propertyRegions[v].fraction) / 100);
        if(cache->propertyRegions[v].warm_bound > cache->propertyRegions[v].upper_bound)
        {
            cache->propertyRegions[v].warm_bound = cache->propertyRegions[v].upper_bound;
        }
    }
}

// ********************************************************************************************
// ***************               Stats output                                    **************
// ********************************************************************************************

void printStatsCache(struct Cache *cache)
{
    float missRate = (double)((getWM(cache) + getRM(cache)) * 100) / (getReads(cache) + getWrites(cache)); //calculate miss rate
    missRate       = roundf(missRate * 100) / 100;                                                //rounding miss rate

    float missRateRead = (double)((getRM(cache)) * 100) / (getReads(cache));   //calculate miss rate
    missRateRead       = roundf(missRateRead * 100) / 100;                     //rounding miss rate

    float missRateWrite = (double)((getWM(cache)) * 100) / (getWrites(cache)); //calculate miss rate
    missRateWrite       = roundf(missRateWrite * 100) / 100;                   //rounding miss rate

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Simulation results (Cache)");
    printf(" -----------------------------------------------------\n");

    switch(cache->policy)
    {
    case LRU_POLICY:
        printf("| %-51s | \n", "LRU_POLICY");
        break;
    case LFU_POLICY:
        printf("| %-51s | \n", "LFU_POLICY");
        break;
    case GRASP_POLICY:
        printf("| %-51s | \n", "GRASP_POLICY");
        break;
    case SRRIP_POLICY:
        printf("| %-51s | \n", "SRRIP_POLICY");
        break;
    case PIN_POLICY:
        printf("| %-51s | \n", "PIN_POLICY");
        break;
    case PLRU_POLICY:
        printf("| %-51s | \n", "PLRU_POLICY");
        break;
    case GRASPXP_POLICY:
        printf("| %-51s | \n", "GRASPXP_POLICY");
        break;
    default :
        printf("| %-51s | \n", "LRU_POLICY");
    }

    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Cache Size (KB)", cache->size / 1024 );
    printf("| %-21s | %'-27lu | \n", "Block Size",    cache->lineSize);
    printf("| %-21s | %'-27lu | \n", "Associativity", cache->assoc);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads/Writes", (getReads(cache) + getWrites(cache)) );
    printf("| %-21s | %'-27lu | \n", "Reads/Writes misses", (getWM(cache) + getRM(cache)));
    printf("| %-21s | %-27.2f | \n", "Miss rate(%)", missRate);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads", getReads(cache) );
    printf("| %-21s | %'-27lu | \n", "Read misses", getRM(cache) );
    printf("| %-21s | %-27.2f | \n", "Rd Miss rate(%)", missRateRead);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writes", getWrites(cache) );
    printf("| %-21s | %'-27lu | \n", "Write misses", getWM(cache) );
    printf("| %-21s | %-27.2f | \n", "Wrt Miss rate(%)", missRateWrite);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writebacks", getWB(cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Evictions", getEVC(cache) );
    printf(" -----------------------------------------------------\n");
}

void printStatsGraphReuse(struct Cache *cache, uint32_t *degrees)
{
    uint32_t  i = 0;
    uint32_t  v = 0;
    uint64_t  avgDegrees = 0;
    uint32_t  num_buckets = cache->num_buckets;
    uint32_t  num_vertices = cache->numVertices;

    uint64_t *thresholds;
    uint64_t *thresholds_count;
    uint64_t *thresholds_totalAccesses;
    uint64_t *thresholds_totalDegrees;
    uint64_t *thresholds_totalReuses;
    uint64_t *thresholds_totalMisses;

    float *thresholds_avgAccesses;
    float *thresholds_avgDegrees;
    float *thresholds_avgReuses;
    float *thresholds_avgMisses;
    float  Access_percentage;
    float  Count_percentage;

    uint64_t thresholds_totalAccess  = 0;
    uint64_t thresholds_totalCount   = 0;
    uint64_t thresholds_totalDegree  = 0;
    uint64_t thresholds_totalReuse   = 0;
    uint64_t thresholds_totalMiss    = 0;

    thresholds               = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_count         = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_totalDegrees  = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_totalReuses   = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_totalMisses   = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_totalAccesses = (uint64_t *) my_malloc(num_buckets * sizeof(uint64_t));
    thresholds_avgAccesses   = (float *) my_malloc(num_buckets * sizeof(float));
    thresholds_avgDegrees    = (float *) my_malloc(num_buckets * sizeof(float));
    thresholds_avgReuses     = (float *) my_malloc(num_buckets * sizeof(float));
    thresholds_avgMisses     = (float *) my_malloc(num_buckets * sizeof(float));

    for (i = 0; i < num_buckets; ++i)
    {
        thresholds[i]               = 0;
        thresholds_count[i]         = 0;
        thresholds_totalDegrees[i]  = 0;
        thresholds_totalReuses[i]   = 0;
        thresholds_totalMisses[i]   = 0;
        thresholds_avgDegrees[i]    = 0.0f;
        thresholds_avgReuses[i]     = 0.0f;
        thresholds_avgMisses[i]     = 0.0f;
        thresholds_totalAccesses[i] = 0;
        thresholds_avgAccesses[i]   = 0.0f;
    }

    for (v = 0; v < num_vertices; ++v)
    {
        avgDegrees +=  degrees[v];
    }

    avgDegrees /= num_vertices;

    // START initialize thresholds
    if(avgDegrees <= 1)
        thresholds[0] = 1;
    else
        thresholds[0] = (avgDegrees / 2);
    for ( i = 1; i < (num_buckets - 1); ++i)
    {
        thresholds[i] = thresholds[i - 1] * 2;
    }
    thresholds[num_buckets - 1] = UINT32_MAX;
    // END initialize thresholds


    // collect stats perbucket
    for (v = 0; v < num_vertices; ++v)
    {
        for ( i = 0; i < num_buckets; ++i)
        {
            if(degrees[v] <= thresholds[i])
            {
                if(cache->vertices_accesses[v])
                {
                    thresholds_count[i] += 1;
                    thresholds_totalDegrees[i]  += degrees[v];
                }
                thresholds_totalReuses[i]   += cache->vertices_total_reuse[v];
                thresholds_totalMisses[i]   += cache->verticesMiss[v];
                thresholds_totalAccesses[i] += cache->vertices_accesses[v];
                break;
            }
        }
    }

    // collect stats perbucket
    for ( i = 0; i < num_buckets; ++i)
    {
        if(thresholds_count[i])
        {
            thresholds_avgDegrees[i]   = (float)thresholds_totalDegrees[i]  / thresholds_count[i];
            // thresholds_avgReuses[i]    = (float)thresholds_totalReuses[i]   / thresholds_totalAccesses[i];
            thresholds_avgMisses[i]    = (float)thresholds_totalMisses[i]    / thresholds_count[i];
            thresholds_avgAccesses[i]  = (float)thresholds_totalAccesses[i]  / thresholds_count[i];
        }

        if(thresholds_totalAccesses[i])
        {
            thresholds_avgReuses[i]    = (float)thresholds_totalReuses[i]   / thresholds_totalAccesses[i];
        }

        thresholds_totalAccess += thresholds_totalAccesses[i];
        thresholds_totalCount  += thresholds_count[i];
        thresholds_totalDegree += thresholds_totalDegrees[i];
        thresholds_totalReuse  += thresholds_totalReuses[i];
        thresholds_totalMiss   += thresholds_totalMisses[i];
    }

    printf(" -----------------------------------------------------------------------------------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s | \n", "<= Threshold", "Nodes(%)", "totalAccess(E)", "(%)Accesses", "avgDegrees", "avgReuse CL/Axs", "avgMisses");
    printf(" -----------------------------------------------------------------------------------------------------------------------------\n");
    for ( i = 0; i < num_buckets; ++i)
    {
        Access_percentage = 100 * (float)(thresholds_totalAccesses[i] / (float)thresholds_totalAccess);
        Count_percentage  = 100 * (float)(thresholds_count[i] / (float)num_vertices);
        printf("| %-15lu | %-15.2f | %-15lu | %-15.2f | %-15.2f | %-15.2f | %-15.2f |\n", thresholds[i], Count_percentage, thresholds_totalAccesses[i], Access_percentage, thresholds_avgDegrees[i], thresholds_avgReuses[i], thresholds_avgMisses[i]);
    }
    printf(" -----------------------------------------------------------------------------------------------------------------------------\n");
    printf("| %-15s | %-15s | %-15s | %-15s | %-15s |  %-15s | %-15s | \n", "avgDegrees", "Total Count(V)", "totalAccess", "totalDegrees",  "avgDegrees", "totalReuse", "totalMisses");
    printf(" -----------------------------------------------------------------------------------------------------------------------------\n");
    printf("| %-15lu | %-15lu | %-15lu | %-15lu | %-15lu |  %-15lu | %-15lu |\n", avgDegrees, thresholds_totalCount, thresholds_totalAccess, thresholds_totalDegree, avgDegrees, thresholds_totalReuse, thresholds_totalMiss);
    printf(" -----------------------------------------------------------------------------------------------------------------------------\n");

    free(thresholds);
    free(thresholds_count);
    free(thresholds_totalDegrees);
    free(thresholds_totalReuses);
    free(thresholds_totalMisses);
    free(thresholds_avgDegrees);
    free(thresholds_avgReuses);
    free(thresholds_avgMisses);
    free(thresholds_totalAccesses);
    free(thresholds_avgAccesses);
}

void printStatsGraphCache(struct Cache *cache, uint32_t *in_degree, uint32_t *out_degree)
{
    printStatsCache(cache);
    printf("\n======================  Reuse stats Out Degree =======================\n");
    printStatsGraphReuse(cache, out_degree);
    // printf("\n======================  Reuse stats In Degree  =======================\n");
    // printStatsGraphReuse(cache, in_degree);
    uint32_t v;

    printf("\n=====================      Property Regions          =================\n");
    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        printf(" -----------------------------------------------------\n");
        printf("| %-25s | %-24u| \n", "ID", v);
        printf(" -----------------------------------------------------\n");
        printf("| %-25s | %-24u| \n", "size", cache->propertyRegions[v].size);
        printf("| %-25s | %-24u| \n", "data_type_size", cache->propertyRegions[v].data_type_size);
        printf("| %-25s | 0x%-22lx| \n", "base_address", cache->propertyRegions[v].base_address);
        printf(" -----------------------------------------------------\n");
        printf("| %-25s | %-24u| \n", "fraction", cache->propertyRegions[v].fraction );
        printf("| %-25s | 0x%-22lx| \n", "lower_bound", cache->propertyRegions[v].lower_bound );
        printf("| %-25s | 0x%-22lx| \n", "hot_bound", cache->propertyRegions[v].hot_bound );
        printf("| %-25s | 0x%-22lx| \n", "warm_bound", cache->propertyRegions[v].warm_bound );
        printf("| %-25s | 0x%-22lx| \n", "upper_bound", cache->propertyRegions[v].upper_bound );
        printf(" -----------------------------------------------------\n");
    }
}

void printStatsAccelGraphCache(struct AccelGraphCache *cache, uint32_t *in_degree, uint32_t *out_degree)
{
    //rounding miss rate

    uint64_t readsHits_hot    = getReads(cache->hot_cache)  - getRM(cache->hot_cache);
    uint64_t readsHits_warm   = getReads(cache->warm_cache) - getRM(cache->warm_cache);

    uint64_t readsMisses_cold = getRM(cache->cold_cache);

    uint64_t writesHits_hot   = getWrites(cache->hot_cache)  - getWM(cache->hot_cache);
    uint64_t writesHits_warm  = getWrites(cache->warm_cache) - getWM(cache->warm_cache);

    uint64_t writesMisses_cold = getWM(cache->cold_cache);

    uint64_t ReadWrite_total   = getReads(cache->cold_cache) + getWrites(cache->cold_cache) + readsHits_hot + readsHits_warm + writesHits_hot + writesHits_warm;
    uint64_t ReadWriteMisses_total = readsMisses_cold + writesMisses_cold;

    uint64_t Read_total       = getReads(cache->cold_cache) + readsHits_hot + readsHits_warm;
    uint64_t ReadMisses_total = readsMisses_cold ;

    uint64_t Write_total       = getWrites(cache->cold_cache) + writesHits_hot + writesHits_warm;
    uint64_t WriteMisses_total =  writesMisses_cold;

    float missRate = (double)(ReadWriteMisses_total * 100) / (ReadWrite_total); //calculate miss rate
    missRate       = roundf(missRate * 100) / 100;

    float missRateRead = (double)((ReadMisses_total) * 100) / (Read_total); //calculate miss rate
    missRateRead       = roundf(missRateRead * 100) / 100;

    float missRateWrite = (double)((WriteMisses_total) * 100) / (Write_total); //calculate miss rate
    missRateWrite       = roundf(missRateWrite * 100) / 100;                            //rounding miss rate

    printf("\n====================== cache Stats Accel Graph =======================\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Simulation results (Cache)");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads/Writes", ReadWrite_total);
    printf("| %-21s | %'-27lu | \n", "Reads/Writes misses", ReadWriteMisses_total);
    printf("| %-21s | %-27.2f | \n", "Miss rate(%)", missRate);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads", Read_total);
    printf("| %-21s | %'-27lu | \n", "Read misses", ReadMisses_total );
    printf("| %-21s | %-27.2f | \n", "Rd Miss rate(%)", missRateRead);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writes", Write_total);
    printf("| %-21s | %'-27lu | \n", "Write misses", WriteMisses_total );
    printf("| %-21s | %-27.2f | \n", "Wrt Miss rate(%)", missRateWrite);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writebacks", getWB(cache->cold_cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Evictions", getEVC(cache->cold_cache) );
    printf(" -----------------------------------------------------\n");

    printf("\n================ cache Stats ((PSL) cold_cache Stats) ================\n");
    printStatsGraphCache(cache->cold_cache, in_degree, out_degree);
    printf("\n===================  cache Stats (warm_cache Stats) ==================\n");
    printStatsGraphCache(cache->warm_cache, in_degree, out_degree);
    printf("\n===================  cache Stats (hot_cache Stats)  ==================\n");
    printStatsGraphCache(cache->hot_cache, in_degree, out_degree);
}

void printStatsDoubleTaggedCache(struct DoubleTaggedCache *cache, uint32_t *in_degree, uint32_t *out_degree)
{
    printStatsAccelGraphCache(cache->accel_graph, in_degree, out_degree);
    printf("\n======================================================================\n");
    printf("\n===================== cache Stats (ref_cache Stats)  =================\n");
    printStatsGraphCache(cache->ref_cache, in_degree, out_degree);
    // if(cache->ref_cache->policy == GRASP_POLICY)
    // {
    // uint32_t v;
    // printf("\n=====================      Property Regions          =================\n");

    // for (v = 0; v < cache->ref_cache->numPropertyRegions; ++v)
    // {
    //     printf(" -----------------------------------------------------\n");
    //     printf("| %-25s | %-24u| \n", "ID", v);
    //     printf(" -----------------------------------------------------\n");
    //     printf("| %-25s | %-24u| \n", "size", cache->ref_cache->propertyRegions[v].size);
    //     printf("| %-25s | %-24u| \n", "data_type_size", cache->ref_cache->propertyRegions[v].data_type_size);
    //     printf("| %-25s | 0x%-22lx| \n", "base_address", cache->ref_cache->propertyRegions[v].base_address);
    //     printf(" -----------------------------------------------------\n");
    //     printf("| %-25s | %-24u| \n", "fraction", cache->ref_cache->propertyRegions[v].fraction );
    //     printf("| %-25s | 0x%-22lx| \n", "lower_bound", cache->ref_cache->propertyRegions[v].lower_bound );
    //     printf("| %-25s | 0x%-22lx| \n", "hot_bound", cache->ref_cache->propertyRegions[v].hot_bound );
    //     printf("| %-25s | 0x%-22lx| \n", "warm_bound", cache->ref_cache->propertyRegions[v].warm_bound );
    //     printf("| %-25s | 0x%-22lx| \n", "upper_bound", cache->ref_cache->propertyRegions[v].upper_bound );
    //     printf(" -----------------------------------------------------\n");
    // }
    // }
}

void printStats(struct Cache *cache)
{
    float missRate = (double)((getWM(cache) + getRM(cache)) * 100) / (cache->currentCycle_cache); //calculate miss rate
    missRate = roundf(missRate * 100) / 100;                            //rounding miss rate

    float missRatePrefetch = (double)(( getRMPrefetch(cache)) * 100) / (cache->currentCycle_preftcher); //calculate miss rate
    missRatePrefetch = roundf(missRatePrefetch * 100) / 100;

    uint32_t v;
    uint32_t i;

    printf("\n -----------------------------------------------------\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Simulation results (Cache)");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads", getReads(cache) );
    printf("| %-21s | %'-27lu | \n", "Read misses", getRM(cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writes", getWrites(cache) );
    printf("| %-21s | %'-27lu | \n", "Write misses", getWM(cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27.2f | \n", "Miss rate(%)", missRate);
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Writebacks", getWB(cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Evictions", getEVC(cache) );
    printf(" -----------------------------------------------------\n");


    if(cache->policy == GRASP_POLICY)
    {
        printf("\n=====================      Property Regions          =================\n");

        for (v = 0; v < cache->numPropertyRegions; ++v)
        {
            printf(" -----------------------------------------------------\n");
            printf("| %-25s | %-24u| \n", "ID", v);
            printf(" -----------------------------------------------------\n");
            printf("| %-25s | %-24u| \n", "size", cache->propertyRegions[v].size);
            printf("| %-25s | %-24u| \n", "data_type_size", cache->propertyRegions[v].data_type_size);
            printf("| %-25s | 0x%-22lx| \n", "base_address", cache->propertyRegions[v].base_address);
            printf(" -----------------------------------------------------\n");
            printf("| %-25s | %-24u| \n", "fraction", cache->propertyRegions[v].fraction );
            printf("| %-25s | 0x%-22lx| \n", "lower_bound", cache->propertyRegions[v].lower_bound );
            printf("| %-25s | 0x%-22lx| \n", "hot_bound", cache->propertyRegions[v].hot_bound );
            printf("| %-25s | 0x%-22lx| \n", "warm_bound", cache->propertyRegions[v].warm_bound );
            printf("| %-25s | 0x%-22lx| \n", "upper_bound", cache->propertyRegions[v].upper_bound );
            printf(" -----------------------------------------------------\n");
        }
    }

    double    avgVerticesreuse  = 0;
    uint64_t  numVerticesMiss   = 0;
    uint64_t  totalVerticesMiss = 0;
    uint64_t  accVerticesAccess = 0;
    // uint64_t   minReuse = 0;
    // uint32_t  maxVerticesMiss = 0;
    // uint32_t  maxNode = 0;


    for(i = 0; i < cache->numVertices; i++)
    {
        if(cache->verticesMiss[i] > 5)
        {
            numVerticesMiss++;
            totalVerticesMiss += cache->verticesMiss[i];
        }

        if(cache->vertices_accesses[i])
        {

            // printf("%u. Average reuse:             re %u acc %u\n", i,cache->vertices_total_reuse[i],cache->vertices_accesses[i]);

            avgVerticesreuse += cache->vertices_total_reuse[i] / cache->vertices_accesses[i];
            accVerticesAccess++;
        }

    }

    avgVerticesreuse /= accVerticesAccess;

    float MissNodesRatioReadMisses = (((double)numVerticesMiss / cache->numVertices) * 100.0);
    float ratioReadMissesMissNodes = (((double)totalVerticesMiss / getRM(cache)) * 100.0);

    printf("============ Graph Stats (Nodes cause highest miss stats) ============\n");
    printf("01. number of nodes:                          %lu\n", numVerticesMiss);
    printf("02. number of read misses:                    %lu\n", totalVerticesMiss);
    printf("03. ratio from total nodes :                  %.2f%%\n", MissNodesRatioReadMisses);
    printf("04. ratio from total read misses:             %.2f%%\n", ratioReadMissesMissNodes);
    printf("===================== Graph Stats (Other Stats) ======================\n");
    printf("05. Average reuse:             %.2f%%\n", avgVerticesreuse);



    // char *fname_txt = (char *) malloc((strlen(fname) + 20) * sizeof(char));
    // char *fname_topMisses = (char *) malloc((strlen(fname) + 20) * sizeof(char));


    // fname_txt = strcpy (fname_txt, fname);
    // fname_topMisses  = strcat (fname_txt, ".top");// out-degree

    // FILE *fptr;
    // fptr = fopen(fname_topMisses, "w");

    // for(i = 0; i < numVertices; i++)
    // {
    //     if(verticesMiss[i] > 100)
    //         fprintf(fptr, "%u %u\n", i, verticesMiss[i]);
    // }

    // fclose(fptr);


}

// void printStatsDoubleTaggedCache(struct DoubleTaggedCache *cache)
// {
//     for (v = 0; v < numPropertyRegions; ++v)
//     {
//         printf(" -----------------------------------------------------\n");
//         printf("| %-25s | %-24u| \n", "size", propertyMetaData[v].size);
//         printf("| %-25s | %-24u| \n", "data_type_size", propertyMetaData[v].data_type_size);
//         printf("| %-25s | 0x%-22lx| \n", "base_address", propertyMetaData[v].base_address);
//         printf(" -----------------------------------------------------\n");
//     }
//     printf("\n===================== cache Stats (cold_cache Stats) =================\n");
//     printStats(cache->cold_cache);
//     printf("\n===================== cache Stats (warm_cache Stats) =================\n");
//     printStats(cache->warm_cache);
//     printf("\n===================== cache Stats (hot_cache Stats)  =================\n");
//     printStats(cache->hot_cache);
//     printf("\n===================== cache Stats (ref_cache Stats)  =================\n");
//     printStats(cache->ref_cache);
// }
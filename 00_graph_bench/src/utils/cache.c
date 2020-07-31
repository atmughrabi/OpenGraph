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
    cacheLine->seq = 0;
    cacheLine->tag = 0;    //useful function
    cacheLine->Flags = INVALID;
    cacheLine->RRPV = RRPV_INIT;
    cacheLine->freq = 0;
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


struct DoubleTaggedCache *newDoubleTaggedCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices, uint32_t policy, uint32_t numPropertyRegions)
{
    struct DoubleTaggedCache *cache = (struct DoubleTaggedCache *) my_malloc(sizeof(struct DoubleTaggedCache));

    cache->cold_cache = newCache( l1_size, l1_assoc, blocksize, num_vertices, LRU_POLICY, numPropertyRegions);
    cache->warm_cache = newCache( l1_size / 2, 8, 4, num_vertices, policy, numPropertyRegions);
    cache->hot_cache  = newCache( l1_size / 2, 8, 4, num_vertices, policy, numPropertyRegions);
    cache->ref_cache  = newCache( l1_size, l1_assoc, blocksize, num_vertices, policy, numPropertyRegions);

    return cache;
}

void initDoubleTaggedCacheRegion(struct DoubleTaggedCache *cache, struct PropertyMetaData *propertyMetaData)
{
    initialzeCachePropertyRegions (cache->cold_cache, propertyMetaData);
    initialzeCachePropertyRegions (cache->warm_cache, propertyMetaData);
    initialzeCachePropertyRegions (cache->hot_cache, propertyMetaData);
    initialzeCachePropertyRegions (cache->ref_cache, propertyMetaData);
}

void freeDoubleTaggedCache(struct DoubleTaggedCache *cache)
{
    if(cache)
    {
        freeCache(cache->cold_cache);
        freeCache(cache->warm_cache);
        freeCache(cache->hot_cache);
        freeCache(cache->ref_cache);
        free(cache);
    }
}


struct Cache *newCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices, uint32_t policy, uint32_t numPropertyRegions)
{

    uint64_t i;

    struct Cache *cache = ( struct Cache *) my_malloc(sizeof(struct Cache));
    initCache(cache, l1_size, l1_assoc, blocksize, policy);

    cache->numPropertyRegions  = numPropertyRegions;
    cache->propertyRegions     = (struct PropertyRegion *)my_malloc(sizeof(struct PropertyRegion) * numPropertyRegions);

    for(i = 0; i < numPropertyRegions; i++)
    {
        cache->propertyRegions[i].upper_bound = 0;
        cache->propertyRegions[i].hot_bound   = 0;
        cache->propertyRegions[i].warm_bound  = 0;
        cache->propertyRegions[i].lower_bound = 0;
    }

    cache->numVertices  = num_vertices;
    cache->verticesMiss = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->verticesHit  = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_base_reuse  = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_total_reuse = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_accesses    = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);

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
        cache->cacheLines[i] = (struct CacheLine *) my_malloc(cache->assoc * sizeof(struct CacheLine));
        for(j = 0; j < cache->assoc; j++)
        {
            invalidate(&(cache->cacheLines[i][j]));
        }
    }
}

void online_cache_graph_stats(struct Cache *cache, uint32_t node)
{
    uint32_t first_Access = 0;

    cache->vertices_accesses[node]++;
    cache->access_counter++;

    if(cache->vertices_base_reuse[node] == 0)
        first_Access = 1;

    if(first_Access)
    {
        cache->vertices_total_reuse[node] = 1;
        cache->vertices_base_reuse[node]  = cache->access_counter;
    }
    else
    {
        cache->vertices_total_reuse[node] += (cache->access_counter - cache->vertices_base_reuse[node]);
        cache->vertices_base_reuse[node]   = cache->access_counter;
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
    case GRASP_POLICY:
        updateInsertGRASP(cache, line);
        break;
    case LFU_POLICY:
        updateInsertLFU(cache, line);
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
    case GRASP_POLICY:
        updatePromoteGRASP(cache, line);
        break;
    case LFU_POLICY:
        updatePromoteLFU(cache, line);
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
    case GRASP_POLICY:
        victim = getVictimGRASP(cache, addr);
        break;
    case LFU_POLICY:
        victim = getVictimLFU(cache, addr);
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
        if(isValid(&(cache->cacheLines[i][j])) == 0) return &(cache->cacheLines[i][j]);
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
        if(isValid(&(cache->cacheLines[i][j])) == 0) return &(cache->cacheLines[i][j]);
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
        if(isValid(&(cache->cacheLines[i][j])) == 0) return &(cache->cacheLines[i][j]);
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

    // not in the paper optimizaiton
    if (min < DEFAULT_INSERT_RRPV)
    {
        int diff = DEFAULT_INSERT_RRPV - min;
        for(j = 1; j < cache->assoc; j++)
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
        else
        {
            cache->readMisses++;
        }

        struct CacheLine *newline = fillLine(cache, addr);
        if(op == 'w')
            setFlags(newline, DIRTY);

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
// ***************               GRASP Policy                                    **************
// ********************************************************************************************


// ********************************************************************************************
// ***************               GRASP Initializaiton                            **************
// ********************************************************************************************

void initialzeCachePropertyRegions (struct Cache *cache, struct PropertyMetaData *propertyMetaData)
{
    uint32_t v;
    uint64_t total_properties_size = 0;
    // uint32_t property_fraction = 100 / cache->numPropertyRegions; //classical vs ratio of array size in bytes

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        total_properties_size += (propertyMetaData[v].size * propertyMetaData[v].data_type_size);
    }

    for (v = 0; v < cache->numPropertyRegions; ++v)
    {
        // cache->propertyRegions[v].fraction    = property_fraction; // classical vs ratio of array size in bytes
        cache->propertyRegions[v].fraction    = ((propertyMetaData[v].size * propertyMetaData[v].data_type_size) * 100) / total_properties_size;

        cache->propertyRegions[v].lower_bound = propertyMetaData[v].base_address;
        cache->propertyRegions[v].upper_bound = propertyMetaData[v].base_address + (propertyMetaData[v].size * propertyMetaData[v].data_type_size);

        cache->propertyRegions[v].hot_bound = cache->propertyRegions[v].lower_bound + ((cache->size * cache->propertyRegions[v].fraction) / 100);
        if(cache->propertyRegions[v].hot_bound > cache->propertyRegions[v].upper_bound)
        {
            cache->propertyRegions[v].hot_bound = cache->propertyRegions[v].upper_bound;
        }

        cache->propertyRegions[v].warm_bound = cache->propertyRegions[v].hot_bound + ((cache->size * cache->propertyRegions[v].fraction) / 100);
        if(cache->propertyRegions[v].warm_bound > cache->propertyRegions[v].upper_bound)
        {
            cache->propertyRegions[v].warm_bound = cache->propertyRegions[v].upper_bound;
        }
    }

}

// ********************************************************************************************
// ***************               Stats output                                    **************
// ********************************************************************************************


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
        for (v = 0; v < cache->numPropertyRegions; ++v)
        {
            printf(" -----------------------------------------------------\n");
            printf("| %-25s | %-24u| \n", "ID", v);
            printf(" -----------------------------------------------------\n");
            printf("| %-25s | %-24u| \n", "fraction", cache->propertyRegions[v].fraction );
            printf("| %-25s | 0x%-22lx| \n", "lower_bound", cache->propertyRegions[v].lower_bound );
            printf("| %-25s | 0x%-22lx| \n", "hot_bound", cache->propertyRegions[v].hot_bound );
            printf("| %-25s | 0x%-22lx| \n", "warm_bound", cache->propertyRegions[v].warm_bound );
            printf("| %-25s | 0x%-22lx| \n", "upper_bound", cache->propertyRegions[v].upper_bound );
            printf(" -----------------------------------------------------\n");
        }
    }
    uint64_t  numVerticesMiss = 0;
    uint64_t  totalVerticesMiss = 0;
    double  avgVerticesreuse = 0;
    uint64_t   accVerticesAccess = 0;
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
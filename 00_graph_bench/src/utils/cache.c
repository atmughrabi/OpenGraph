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


struct DoubleTaggedCache *newDoubleTaggedCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices)
{

    struct DoubleTaggedCache *cache = (struct DoubleTaggedCache *) my_malloc(sizeof(struct DoubleTaggedCache));

    cache->cold_cache = newCache( l1_size, l1_assoc, blocksize, num_vertices);
    cache->warm_cache = newCache( l1_size / 2, 8, 4, num_vertices);
    cache->hot_cache  = newCache( l1_size / 2, 8, 4, num_vertices);
    cache->ref_cache  = newCache( l1_size, l1_assoc, blocksize, num_vertices);

    return cache;

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


struct Cache *newCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices)
{

    uint64_t i;

    struct Cache *cache = ( struct Cache *) my_malloc(sizeof(struct Cache));

    cache->numVertices = num_vertices;
    initCache(cache, l1_size, l1_assoc, blocksize, POLICY);

    cache->verticesMiss = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->verticesHit = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_base_reuse = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_total_reuse = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);
    cache->vertices_accesses = (uint32_t *)my_malloc(sizeof(uint32_t) * num_vertices);


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

    cache->policy     = (uint64_t)(p);
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

void Access(struct Cache *cache, uint64_t addr, unsigned char op, uint32_t node)
{

    online_cache_graph_stats(cache, node);
    cache->currentCycle++;/*per cache global counter to maintain LRU order

    among cache ways, updated on every cache access*/

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
            cache->verticesMiss[node]++;
        }
        else
        {
            cache->readMisses++;
            cache->verticesMiss[node]++;
        }

        struct CacheLine *newline = fillLine(cache, addr);
        if(op == 'w')
            setFlags(newline, DIRTY);

    }
    else
    {
        /**since it's a hit, update LRU and update dirty flag**/
        updatePolicy(cache, line);
        if(op == 'w')
            setFlags(line, DIRTY);

        cache->verticesHit[node]++;
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
        updatePolicy(cache, line);
    }
}

/*look up line*/
struct CacheLine *findLine(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *line = NULL;

    switch(cache->policy)
    {
    case LRU_POLICY:
        line = findLineLRU(cache, addr);
        break;
    case GRASP_POLICY:
        line = findLineGRASP(cache, addr);
        break;
    case LFU_POLICY:
        line = findLineLFU(cache, addr);
        break;
    default :
        line = findLineLRU(cache, addr);
    }

    return line;
}

struct CacheLine *findLineLRU(struct Cache *cache, uint64_t addr)
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

struct CacheLine *findLineLFU(struct Cache *cache, uint64_t addr)
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

struct CacheLine *findLineGRASP(struct Cache *cache, uint64_t addr)
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

void updatePolicy(struct Cache *cache, struct CacheLine *line)
{
    switch(cache->policy)
    {
    case LRU_POLICY:
        updateLRU(cache, line);
        break;
    case GRASP_POLICY:
        updateGRASP(cache, line);
        break;
    case LFU_POLICY:
        updateLFU(cache, line);
        break;
    default :
        updateLRU(cache, line);
    }
}

/*upgrade LRU line to be MRU line*/
void updateLRU(struct Cache *cache, struct CacheLine *line)
{
    setSeq(line, cache->currentCycle);
}

void updateLFU(struct Cache *cache, struct CacheLine *line)
{
    uint8_t freq = getFreq(line);
    if(freq < FREQ_MAX)
        freq++;
    setFreq(line, freq);
}

void updateGRASP(struct Cache *cache, struct CacheLine *line)
{
    setRRPV(line, cache->currentCycle);
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *getLRU(struct Cache *cache, uint64_t addr)
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

struct CacheLine *getLFU(struct Cache *cache, uint64_t addr)
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

struct CacheLine *getGRASP(struct Cache *cache, uint64_t addr)
{

}

/*find a victim, move it to MRU position*/
struct CacheLine *findLineToReplace(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *victim = NULL;

    switch(cache->policy)
    {
    case LRU_POLICY:
        victim = findLineToReplaceLRU(cache, addr);
        break;
    case GRASP_POLICY:
        victim = findLineToReplaceGRASP(cache, addr);
        break;
    case LFU_POLICY:
        victim = findLineToReplaceLFU(cache, addr);
        break;
    default :
        victim = findLineToReplaceLRU(cache, addr);
    }

    return victim;
}

struct CacheLine *findLineToReplaceLRU(struct Cache *cache, uint64_t addr)
{
    struct CacheLine  *victim = getLRU(cache, addr);
    updateLRU(cache, victim);

    return (victim);
}

struct CacheLine *findLineToReplaceLFU(struct Cache *cache, uint64_t addr)
{
    struct CacheLine  *victim = getLFU(cache, addr);
    updateLFU(cache, victim);

    return (victim);
}

struct CacheLine *findLineToReplaceGRASP(struct Cache *cache, uint64_t addr)
{

}

/*allocate a new line*/
struct CacheLine *fillLine(struct Cache *cache, uint64_t addr)
{
    struct CacheLine *victim = NULL;

    switch(cache->policy)
    {
    case LRU_POLICY:
        victim = fillLineLRU(cache, addr);
        break;
    case GRASP_POLICY:
        victim = fillLineGRASP(cache, addr);
        break;
    case LFU_POLICY:
        victim = fillLineLFU(cache, addr);
        break;
    default :
        victim = fillLineLRU(cache, addr);
    }

    return victim;
}

struct CacheLine *fillLineLRU(struct Cache *cache, uint64_t addr)
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
    setFlags(victim, VALID);


    /**note that this cache line has been already
       upgraded to MRU in the previous function (findLineToReplace)**/

    return victim;
}

struct CacheLine *fillLineLFU(struct Cache *cache, uint64_t addr)
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
    setFlags(victim, VALID);


    /**note that this cache line has been already
       upgraded to MRU in the previous function (findLineToReplace)**/

    return victim;

}

struct CacheLine *fillLineGRASP(struct Cache *cache, uint64_t addr)
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
    setFlags(victim, VALID);


    /**note that this cache line has been already
       upgraded to MRU in the previous function (findLineToReplace)**/

    return victim;

}

void printStats(struct Cache *cache)
{



    float missRate = (double)((getWM(cache) + getRM(cache)) * 100) / (cache->currentCycle_cache); //calculate miss rate
    missRate = roundf(missRate * 100) / 100;                            //rounding miss rate

    float missRatePrefetch = (double)(( getRMPrefetch(cache)) * 100) / (cache->currentCycle_preftcher); //calculate miss rate
    missRatePrefetch = roundf(missRatePrefetch * 100) / 100;


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

    printf(" -----------------------------------------------------\n");

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Prefetcher Stats");
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %'-27lu | \n", "Reads", getReadsPrefetch(cache) );
    printf("| %-21s | %'-27lu | \n", "Read misses", getRMPrefetch(cache) );
    printf(" -----------------------------------------------------\n");
    printf("| %-21s | %-27.2f | \n", "Effeciency(%)", missRatePrefetch);
    printf(" -----------------------------------------------------\n\n");


    uint64_t  numVerticesMiss = 0;
    uint64_t  totalVerticesMiss = 0;
    double  avgVerticesreuse = 0;
    uint64_t   accVerticesAccess = 0;
    // uint64_t   minReuse = 0;
    // uint32_t  maxVerticesMiss = 0;
    // uint32_t  maxNode = 0;

    uint32_t i;
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
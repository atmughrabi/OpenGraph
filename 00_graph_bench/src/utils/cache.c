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
ulong getTag(struct CacheLine *cacheLine)
{
    return cacheLine->tag;
}

ulong getFlags(struct CacheLine *cacheLine)
{
    return cacheLine->Flags;
}
ulong getSeq(struct CacheLine *cacheLine)
{
    return cacheLine->seq;
}
void setSeq(struct CacheLine *cacheLine, ulong Seq)
{
    cacheLine->seq = Seq;
}
ulong getFreq(struct CacheLine *cacheLine)
{
    return cacheLine->freq;
}
void setFreq(struct CacheLine *cacheLine, ulong freq)
{
    cacheLine->freq = freq;
}
void setFlags(struct CacheLine *cacheLine, ulong flags)
{
    cacheLine->Flags = flags;
}
void setTag(struct CacheLine *cacheLine, ulong a)
{
    cacheLine->tag = a;
}
void invalidate(struct CacheLine *cacheLine)
{
    cacheLine->tag = 0;    //useful function
    cacheLine->Flags = INVALID;
}
uint32_t isValid(struct CacheLine *cacheLine)
{
    return ((cacheLine->Flags) != INVALID);
}


//cache helper functions

ulong calcTag(struct Cache *cache, ulong addr)
{
    return (addr >> (cache->log2Blk) );
}
ulong calcIndex(struct Cache *cache, ulong addr)
{
    return ((addr >> cache->log2Blk) & cache->tagMask);
}
ulong calcAddr4Tag(struct Cache *cache, ulong tag)
{
    return (tag << (cache->log2Blk));
}

ulong getRM(struct Cache *cache)
{
    return cache->readMisses;
}
ulong getWM(struct Cache *cache)
{
    return cache->writeMisses;
}
ulong getReads(struct Cache *cache)
{
    return cache->reads;
}
ulong getWrites(struct Cache *cache)
{
    return cache->writes;
}
ulong getWB(struct Cache *cache)
{
    return cache->writeBacks;
}
ulong getEVC(struct Cache *cache)
{
    return cache->evictions;
}
ulong getRMPrefetch(struct Cache *cache)
{
    return cache->readMissesPrefetch;
}

ulong getReadsPrefetch(struct Cache *cache)
{
    return cache->readsPrefetch;
}
void writeBack(struct Cache *cache, ulong addr)
{
    cache->writeBacks++;
}


struct DoubleTaggedCache *newDoubleTaggedCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices)
{

    struct DoubleTaggedCache *cache = (struct DoubleTaggedCache *) my_malloc(sizeof(struct DoubleTaggedCache));

    cache->cold_cache = newCache( l1_size, l1_assoc, blocksize, num_vertices);
    cache->warm_cache = newCache( l1_size / 2, 8, 4, num_vertices);
    cache->hot_cache = newCache( l1_size / 2, 8, 4, num_vertices);
    cache->ref_cache = newCache( l1_size, l1_assoc, blocksize, num_vertices);

    return cache;

}

void freeDoubleTaggedCache(struct DoubleTaggedCache *cache)
{

    if(cache)
    {
        freeCache(cache->cold_cache);
        freeCache(cache->ref_cache);
        freeCache(cache->warm_cache);
        freeCache(cache->hot_cache);
        free(cache);

    }

}


struct Cache *newCache(uint32_t l1_size, uint32_t l1_assoc, uint32_t blocksize, uint32_t num_vertices)
{

    ulong i;

    struct Cache *cache = ( struct Cache *) my_malloc(sizeof(struct Cache));

    cache->numVertices = num_vertices;
    initCache(cache, l1_size, l1_assoc, blocksize);

    cache->verticesMiss = (uint *)my_malloc(sizeof(uint) * num_vertices);
    cache->verticesHit = (uint *)my_malloc(sizeof(uint) * num_vertices);
    cache->vertices_base_reuse = (uint *)my_malloc(sizeof(uint) * num_vertices);
    cache->vertices_total_reuse = (uint *)my_malloc(sizeof(uint) * num_vertices);
    cache->vertices_accesses = (uint *)my_malloc(sizeof(uint) * num_vertices);


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
    ulong i;

    if(cache)
    {

        if(cache->verticesMiss)
            free(cache->verticesMiss);
        if(cache->verticesHit)
            free(cache->verticesHit);


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

void initCache(struct Cache *cache, int s, int a, int b )
{
    ulong i, j;
    cache->reads = cache->readMisses = cache->readsPrefetch = cache->readMissesPrefetch = cache->writes = cache->evictions = 0;
    cache->writeMisses = cache->writeBacks = cache->currentCycle_preftcher = cache->currentCycle_cache = cache->currentCycle = 0;

    cache->size       = (ulong)(s);
    cache->lineSize   = (ulong)(b);
    cache->assoc      = (ulong)(a);
    cache->sets       = (ulong)((s / b) / a);
    cache->numLines   = (ulong)(s / b);
    cache->log2Sets   = (ulong)(log2(cache->sets));
    cache->log2Blk    = (ulong)(log2(b));
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

void cache_graph_stats(struct Cache *cache, uint node)
{
    uint first_Access = 0;

    cache->vertices_accesses[node]++;
    cache->access_counter++;

    if(cache->vertices_base_reuse[node] == 0)
        first_Access = 1;

    if(first_Access)
    {
        cache->vertices_total_reuse[node] = 1;
        cache->vertices_base_reuse[node] = cache->access_counter;
    }
    else
    {
        cache->vertices_total_reuse[node] += cache->access_counter - cache->vertices_base_reuse[node];
        cache->vertices_base_reuse[node] = cache->access_counter;
    }
}

void Access(struct Cache *cache, ulong addr, uchar op, uint node)
{

    cache_graph_stats(cache, node);
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
        updateLRU(cache, line);
        if(op == 'w')
            setFlags(line, DIRTY);

        cache->verticesHit[node]++;
    }
}

uint32_t checkInCache(struct Cache *cache, ulong addr)
{
    struct CacheLine *line = findLine(cache, addr);

    if(line == NULL)
        return 1;

    // updateLRU(cache, findLine(cache, addr));
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

void Prefetch(struct Cache *cache, ulong addr, uchar op, uint node)
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
        updateLRU(cache, line);
    }
}


/*look up line*/
struct CacheLine *findLine(struct Cache *cache, ulong addr)
{
    ulong i, j, tag, pos;

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

/*upgrade LRU line to be MRU line*/
void updateLRU(struct Cache *cache, struct CacheLine *line)
{
    setSeq(line, cache->currentCycle);
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
struct CacheLine *getLRU(struct Cache *cache, ulong addr)
{
    ulong i, j, victim, min;

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

/*find a victim, move it to MRU position*/
struct CacheLine *findLineToReplace(struct Cache *cache, ulong addr)
{
    struct CacheLine  *victim = getLRU(cache, addr);
    updateLRU(cache, victim);

    return (victim);
}

/*allocate a new line*/
struct CacheLine *fillLine(struct Cache *cache, ulong addr)
{
    ulong tag;

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


    ulong  numVerticesMiss = 0;
    ulong  totalVerticesMiss = 0;
    double  avgVerticesreuse = 0;
    ulong   accVerticesAccess = 0;
    // ulong   minReuse = 0;
    // uint  maxVerticesMiss = 0;
    // uint  maxNode = 0;

    uint i;
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
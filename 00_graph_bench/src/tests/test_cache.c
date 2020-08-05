#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>
#include <stdbool.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "myMalloc.h"
#include <assert.h>


#include "cache.h"

int main(int argc, char *argv[])
{

    uint32_t num_vertices = 0;
    uint32_t numPropertyRegions = 2;
    struct PropertyMetaData *propertyMetaData = (struct PropertyMetaData *) my_malloc(numPropertyRegions * sizeof(struct PropertyMetaData));
    struct DoubleTaggedCache *cache = newDoubleTaggedCache(L1_SIZE,  L1_ASSOC,  BLOCKSIZE, num_vertices, POLICY, 1);

    propertyMetaData[0].base_address = (uint64_t)&riDividedOnDiClause[0];
    propertyMetaData[0].size = graph->num_vertices;
    propertyMetaData[0].data_type_size = sizeof(float);

    propertyMetaData[1].base_address = (uint64_t)&pageRanksNext[0];
    propertyMetaData[1].size = graph->num_vertices;
    propertyMetaData[1].data_type_size = sizeof(float);

    initDoubleTaggedCacheRegion(cache, propertyMetaData);
    setDoubleTaggedCacheThresholdDegreeAvg(cache, graph->vertices->out_degree);



    printStatsDoubleTaggedCache(cache, graph->vertices->in_degree, graph->vertices->out_degree);
    freeDoubleTaggedCache(cache);
    if(propertyMetaData)
        free(propertyMetaData);



    FILE *fp = fopen(argv[1], "rb");
    if ( !fp )
    {
        printf("Could not open file: %s\n", argv[1]);
    }
    // https://github.com/ease-lab/grasp/blob/master/trace-based-simulators/cache.h
    // File format
    // #regions (8 bytes) --> number of regions
    // Followed by #regions records, each record occupying 25+8=33 bytes.
    //      Each record contains name of the region(25 bytes) and value(8 bytes).
    // Followed by #address (8 bytes) --> number of L2 misses
    // Followed by a trace of addresses, 8 bytes each.

    // Region name convention:
    //  Region-1 --> lower bound address for the region
    //  Region-n --> upper bound address for the region
    //  Region-f --> percentage of LLC capacity to be allocated to this region
    //  Any other region can be ignored.
    assert(fp != NULL);
    uint64_t res = 0;
    uint64_t len = 0;
    res = fread(&len, 8, 1, fp);
    char key[25];
    uint64_t val;
    printf("%35s %20s\n", "Region", "Value");
    for ( int i = 0 ; i < len ; i++ )
    {
        res = fread(&key, 1, 25, fp);
        res = fread(&val, 8, 1, fp);
        magic_map[string(key)] = val;
        printf("%35s %20lu\n", key, val);
    }
    printf("\n");

    add_all_regions();

    fread(&total_addresses, 8, 1, fp);

    printf("%35s %20lu\n", "Total Addresses:", total_addresses);

    trace = (uint64_t *) malloc(total_addresses * 8);
    uint64_t bytes = fread(trace, 8, total_addresses, fp);
    printf("%35s %20lu\n", "Elements read:", bytes);
    if ( bytes != total_addresses )
    {
        printf("Could not read all addresses.\n");
    }
    assert(trace);

}
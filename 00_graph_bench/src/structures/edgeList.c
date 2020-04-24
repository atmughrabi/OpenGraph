// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : edgeList.c
// Create : 2019-06-29 12:31:24
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <err.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>


#include "mt19937.h"
#include "myMalloc.h"
#include "graphConfig.h"
#include "edgeList.h"



uint32_t maxTwoIntegers(uint32_t num1, uint32_t num2)
{

    if(num1 >= num2)
        return num1;
    else
        return num2;

}


void writeEdgeListToTXTFile(struct EdgeList *edgeList, const char *fname)
{

    FILE *fp;
    uint32_t i;


    char *fname_txt = (char *) malloc((strlen(fname) + 10) * sizeof(char));

    fname_txt = strcpy (fname_txt, fname);
    fname_txt = strcat (fname_txt, ".txt");

    fp = fopen (fname_txt, "w");

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Vertices (V)");
    printf("| %-51u | \n", edgeList->num_vertices);
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Number of Edges (E)");
    printf("| %-51u | \n", edgeList->num_edges);
    printf(" -----------------------------------------------------\n");


    for(i = 0; i < edgeList->num_edges; i++)
    {

        fprintf(fp, "%u %u\n", edgeList->edges_array_src[i], edgeList->edges_array_dest[i]);

    }

    fclose (fp);
}

// read edge file to edge_array in memory
struct EdgeList *newEdgeList( uint32_t num_edges)
{


    struct EdgeList *newEdgeList = (struct EdgeList *) my_malloc(sizeof(struct EdgeList));
    newEdgeList->edges_array_src = (uint32_t *) my_malloc(num_edges * sizeof(uint32_t));
    newEdgeList->edges_array_dest = (uint32_t *) my_malloc(num_edges * sizeof(uint32_t));

#if WEIGHTED
    newEdgeList->edges_array_weight = (uint32_t *) my_malloc(num_edges * sizeof(uint32_t));
#endif

    uint32_t i;
    #pragma omp parallel for
    for(i = 0; i < num_edges; i++)
    {
        newEdgeList->edges_array_dest[i] = 0;
        newEdgeList->edges_array_src[i] = 0;
#if WEIGHTED
        newEdgeList->edges_array_weight[i] = 0;
#endif
    }

    newEdgeList->num_edges = num_edges;
    newEdgeList->num_vertices = 0;
    // newEdgeList->edges_array = newEdgeList(num_edges);

#if WEIGHTED
    newEdgeList->max_weight = 0;
#endif

    return newEdgeList;

}


struct EdgeList *removeDulpicatesSelfLoopEdges( struct EdgeList *edgeList)
{

    struct EdgeList *tempEdgeList = newEdgeList(edgeList->num_edges);
    uint32_t tempSrc = 0;
    uint32_t tempDest = 0;
    uint32_t tempWeight = 0;
    uint32_t j = 0;
    uint32_t i = 0;


    do
    {
        tempSrc = edgeList->edges_array_src[i];
        tempDest = edgeList->edges_array_dest[i];
#if WEIGHTED
        tempWeight = edgeList->edges_array_weight[i];
#endif
        i++;
    }
    while(tempSrc == tempDest);

        tempEdgeList->edges_array_src[j] = tempSrc;
    tempEdgeList->edges_array_dest[j] = tempDest;
#if WEIGHTED
    tempEdgeList->edges_array_weight[j] = tempWeight;
#endif
    j++;

    for(; i < tempEdgeList->num_edges; i++)
    {
        tempSrc = edgeList->edges_array_src[i];
        tempDest = edgeList->edges_array_dest[i];
#if WEIGHTED
        tempWeight = edgeList->edges_array_weight[i];
#endif
        if(tempSrc != tempDest)
        {
            if(tempEdgeList->edges_array_src[j - 1] != tempSrc || tempEdgeList->edges_array_dest[j - 1] != tempDest )
            {
                tempEdgeList->edges_array_src[j] = tempSrc;
                tempEdgeList->edges_array_dest[j] = tempDest;
#if WEIGHTED
                tempEdgeList->edges_array_weight[j] = tempWeight;
#endif
                j++;
            }
        }
    }

    tempEdgeList->num_edges = j;
    tempEdgeList->num_vertices = edgeList->num_vertices ;
    freeEdgeList(edgeList);
    return tempEdgeList;

}

void freeEdgeList( struct EdgeList *edgeList)
{

    if(edgeList)
    {
        // freeEdgeArray(edgeList->edges_array);
        if(edgeList->edges_array_src)
            free(edgeList->edges_array_src);
        if(edgeList->edges_array_dest)
            free(edgeList->edges_array_dest);

#if WEIGHTED
        if(edgeList->edges_array_weight)
            free(edgeList->edges_array_weight);
#endif

        free(edgeList);
    }


}


char *readEdgeListstxt(const char *fname, uint32_t weighted)
{

    FILE *pText, *pBinary;
    uint32_t size = 0, i;
    uint32_t src = 0, dest = 0;

#if WEIGHTED
    uint32_t weight = 1;
#endif

    char *fname_txt = (char *) malloc((strlen(fname) + 10) * sizeof(char));
    char *fname_bin = (char *) malloc((strlen(fname) + 10) * sizeof(char));

    fname_txt = strcpy (fname_txt, fname);

#if WEIGHTED
    fname_bin = strcat (fname_txt, ".wbin");
#else
    fname_bin = strcat (fname_txt, ".bin");
#endif



    // printf("Filename : %s \n",fname);
    // printf("Filename : %s \n",fname_bin);


    pText = fopen(fname, "r");
    pBinary = fopen(fname_bin, "wb");



    if (pText == NULL)
    {
        err(1, "open: %s", fname);
        return NULL;
    }
    if (pBinary == NULL)
    {
        err(1, "open: %s", fname_bin);
        return NULL;
    }

    while (1)
    {

        // if(size > 48){
#if WEIGHTED
        if(weighted)
        {
            i = fscanf(pText, "%u\t%u\n", &src, &dest);
            // weight = (generateRandInt(mt19937var) % 256) + 1;
            weight = 1;
        }
        else
        {
            i = fscanf(pText, "%u\t%u\t%u\n", &src, &dest, &weight);
        }
#else
        i = fscanf(pText, "%u\t%u\n", &src, &dest);
#endif

        if( i == EOF )
            break;


        fwrite(&src, sizeof (src), 1, pBinary);
        fwrite(&dest, sizeof (dest), 1, pBinary);

#if WEIGHTED
        fwrite(&weight, sizeof (weight), 1, pBinary);
#endif

        size++;
    }


    fclose(pText);
    fclose(pBinary);


    return fname_bin;


}




struct EdgeList *readEdgeListsbin(const char *fname, uint8_t inverse, uint32_t symmetric, uint32_t weighted)
{


    int fd = open(fname, O_RDONLY);
    struct stat fs;
    char *buf_addr;
    uint32_t  *buf_pointer;
    uint32_t  src = 0, dest = 0;
    uint32_t offset;

    if (fd == -1)
    {
        err(1, "open: %s", fname);
        return 0;
    }

    if (fstat(fd, &fs) == -1)
    {
        err(1, "stat: %s", fname);
        return 0;
    }

    /* fs.st_size could have been 0 actually */
    buf_addr = mmap(0, fs.st_size, PROT_WRITE, MAP_PRIVATE, fd, 0);

    if (buf_addr == (void *) -1)
    {
        err(1, "mmap: %s", fname);
        close(fd);
        return 0;
    }


    buf_pointer = (uint32_t *) buf_addr;

#if WEIGHTED
    if(weighted)
        offset = 2;
    else
        offset = 3;

#else
    offset = 2;
#endif

    uint32_t num_edges = (uint64_t)fs.st_size / ((offset) * sizeof(uint32_t));
    // uint32_t num_edges = 32;
    struct EdgeList *edgeList;

#if DIRECTED
    if(symmetric)
    {
        edgeList = newEdgeList((num_edges) * 2);
    }
    else
    {
        edgeList = newEdgeList(num_edges);
    }
#else
    if(symmetric)
    {
        edgeList = newEdgeList((num_edges) * 2);
    }
    else
    {
        edgeList = newEdgeList(num_edges);
    }
#endif

    uint32_t i;
    uint32_t num_vertices = 0;

#if WEIGHTED
    uint32_t max_weight = 0;
#endif

    // #pragma omp parallel for reduction(max:num_vertices)
    for(i = 0; i < num_edges; i++)
    {
        src = buf_pointer[((offset) * i) + 0];
        dest = buf_pointer[((offset) * i) + 1];
        // printf(" %u %lu -> %lu \n",i,src,dest);
#if DIRECTED
        if(!inverse)
        {
            if(symmetric)
            {
                edgeList->edges_array_src[i] = src;
                edgeList->edges_array_dest[i] = dest;
                edgeList->edges_array_src[i + (num_edges)] = dest;
                edgeList->edges_array_dest[i + (num_edges)] = src;

#if WEIGHTED
                if(weighted)
                {
                    edgeList->edges_array_weight[i] = 1;
                    edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
                }
                else
                {
                    edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
                    edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
                }
#endif

            }
            else
            {
                edgeList->edges_array_src[i] = src;
                edgeList->edges_array_dest[i] = dest;

#if WEIGHTED
                if(weighted)
                {
                    edgeList->edges_array_weight[i] = 1;
                }
                else
                {
                    edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
                }
#endif
            } // symmetric
        } // inverse
        else
        {
            if(symmetric)
            {
                edgeList->edges_array_src[i] = dest;
                edgeList->edges_array_dest[i] = src;
                edgeList->edges_array_src[i + (num_edges)] = src;
                edgeList->edges_array_dest[i + (num_edges)] = dest;
#if WEIGHTED
                if(weighted)
                {
                    edgeList->edges_array_weight[i] = 1;
                    edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
                }
                else
                {
                    edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
                    edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
                }
#endif
            }
            else
            {
                edgeList->edges_array_src[i] = dest;
                edgeList->edges_array_dest[i] = src;
#if WEIGHTED
                if(weighted)
                {
                    edgeList->edges_array_weight[i] = 1;
                }
                else
                {
                    edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
                }
#endif
            }// symmetric
        }// inverse
#else
        if(symmetric)
        {
            edgeList->edges_array_src[i] = src;
            edgeList->edges_array_dest[i] = dest;
            edgeList->edges_array_src[i + (num_edges)] = dest;
            edgeList->edges_array_dest[i + (num_edges)] = src;
#if WEIGHTED
            if(weighted)
            {
                edgeList->edges_array_weight[i] = 1;
                edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
            }
            else
            {
                edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
                edgeList->edges_array_weight[i + (num_edges)] = edgeList->edges_array_weight[i];
            }
#endif
        }
        else
        {
            edgeList->edges_array_src[i] = src;
            edgeList->edges_array_dest[i] = dest;
#if WEIGHTED
            if(weighted)
            {
                edgeList->edges_array_weight[i] = 1;
            }
            else
            {
                edgeList->edges_array_weight[i] = buf_pointer[((offset) * i) + 2];
            }
#endif
        }
#endif

        num_vertices = maxTwoIntegers(num_vertices, maxTwoIntegers(edgeList->edges_array_src[i], edgeList->edges_array_dest[i]));

#if WEIGHTED
        max_weight = maxTwoIntegers(max_weight, edgeList->edges_array_weight[i]);
#endif

    }

    edgeList->num_vertices = num_vertices + 1; // max number of veritices Array[0-max]

#if WEIGHTED
    edgeList->max_weight = max_weight;
#endif
    // printf("DONE Reading EdgeList from file %s \n", fname);
    // edgeListPrint(edgeList);

    munmap(buf_addr, fs.st_size);
    close(fd);

    return edgeList;
}


struct EdgeList *readEdgeListsMem( struct EdgeList *edgeListmem,  uint8_t inverse, uint32_t symmetric, uint32_t weighted)
{


    uint32_t num_edges = edgeListmem->num_edges;
    uint32_t  src = 0, dest = 0, weight = 1;;
    struct EdgeList *edgeList;

    edgeList = newEdgeList((num_edges));

    uint32_t i;

    // #pragma omp parallel for
    for(i = 0; i < num_edges; i++)
    {
        src = edgeListmem->edges_array_src[i];
        dest = edgeListmem->edges_array_dest[i];
#if WEIGHTED
        weight = edgeListmem->edges_array_weight[i];
#endif
        // printf(" %u %lu -> %lu \n",src,dest);
#if DIRECTED
        if(!inverse)
        {
            if(symmetric)
            {
                edgeList->edges_array_src[i] = src;
                edgeList->edges_array_dest[i] = dest;


#if WEIGHTED

                edgeList->edges_array_weight[i] = weight;

#endif

            }
            else
            {
                edgeList->edges_array_src[i] = src;
                edgeList->edges_array_dest[i] = dest;

#if WEIGHTED
                if(weighted)
                {
                    edgeList->edges_array_weight[i] = weight;
                }
                else
                {
                    edgeList->edges_array_weight[i] = weight;
                }
#endif
            } // symmetric
        } // inverse
        else
        {
            if(symmetric)
            {
                edgeList->edges_array_src[i] = dest;
                edgeList->edges_array_dest[i] = src;

#if WEIGHTED

                edgeList->edges_array_weight[i] = weight;

#endif
            }
            else
            {
                edgeList->edges_array_src[i] = dest;
                edgeList->edges_array_dest[i] = src;
#if WEIGHTED

                edgeList->edges_array_weight[i] = weight;

#endif
            }// symmetric
        }// inverse
#else
        if(symmetric)
        {
            edgeList->edges_array_src[i] = src;
            edgeList->edges_array_dest[i] = dest;

#if WEIGHTED

            edgeList->edges_array_weight[i] = weight;

#endif
        }
        else
        {
            edgeList->edges_array_src[i] = src;
            edgeList->edges_array_dest[i] = dest;
#if WEIGHTED

            edgeList->edges_array_weight[i] = weight;

#endif
        }
#endif
    }

    edgeList->num_vertices = edgeListmem->num_vertices; // max number of veritices Array[0-max]

#if WEIGHTED
    edgeList->max_weight =  edgeListmem->max_weight;
#endif

    return edgeList;
}

void edgeListPrint(struct EdgeList *edgeList)
{


    printf("number of vertices (V) : %u \n", edgeList->num_vertices);
    printf("number of edges    (E) : %u \n", edgeList->num_edges);

    uint32_t i;
    for(i = 0; i < edgeList->num_edges; i++)
    {
#if WEIGHTED
        printf("%u -> %u w: %d \n", edgeList->edges_array_src[i], edgeList->edges_array_dest[i], edgeList->edges_array_weight[i]);
#else
        printf("%u -> %u \n", edgeList->edges_array_src[i], edgeList->edges_array_dest[i]);
#endif
    }

}
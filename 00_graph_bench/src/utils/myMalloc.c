// -----------------------------------------------------------------------------
//
//      "OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : myMalloc.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#include "myMalloc.h"

extern int errno ;

void *aligned_malloc( size_t size )
{

    // void *dataOut = aligned_alloc(CACHELINE_BYTES, size);
    void *dataOut = NULL ;
    int err = posix_memalign(
      (void**)&dataOut, CACHELINE_BYTES, size);

    assert(err == 0 && "Error, aligned_alloc!");
    // if ( !dataOut )
    // {
    //     fprintf ( stderr, "Error, aligned_alloc: %s\n", strerror ( errno )) ;
    //     exit ( EXIT_FAILURE ) ;
    // }
    return dataOut ;
}

void *regular_malloc( size_t size )
{

    // void *dataOut = aligned_alloc(CACHELINE_BYTES, size);
    void *dataOut = malloc ( size ) ;
    if ( !dataOut )
    {
        fprintf ( stderr, "Error, malloc: %s\n", strerror ( errno )) ;
        exit ( EXIT_FAILURE ) ;
    }
    return dataOut ;
}

void *my_malloc( size_t size )
{

    // void *dataOut = aligned_alloc(CACHELINE_BYTES, size);
    void *dataOut = NULL;
#if ALIGNED
    dataOut =   aligned_malloc(size);
#else
    dataOut =   regular_malloc(size);
#endif

    return dataOut;
}



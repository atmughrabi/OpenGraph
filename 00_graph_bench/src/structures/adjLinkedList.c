// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : adjLinkedList.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:36:13
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#include "adjLinkedList.h"
#include "myMalloc.h"
#include "graphConfig.h"

// A utility function to create a new adjacency list node
struct AdjLinkedListNode *newAdjLinkedListOutNode(uint32_t dest)
{

    struct AdjLinkedListNode *newNode = (struct AdjLinkedListNode *) my_malloc(sizeof(struct AdjLinkedListNode));


    newNode->dest = dest;
#if WEIGHTED
    newNode->weight = 0;
#endif

    newNode->next = NULL;

    return newNode;

}


struct AdjLinkedListNode *newAdjLinkedListInNode( uint32_t src)
{

    struct AdjLinkedListNode *newNode = (struct AdjLinkedListNode *) my_malloc(sizeof(struct AdjLinkedListNode));


    newNode->dest = src;
#if WEIGHTED
    newNode->weight = 0;
#endif

    newNode->next = NULL;

    return newNode;

}
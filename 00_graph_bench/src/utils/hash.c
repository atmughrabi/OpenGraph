// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : hash.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------
#include <stdint.h>
#include "hash.h"

uint32_t magicHash32(uint32_t x)
{
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

uint32_t magicHash32Reverse(uint32_t x)
{
    x = ((x >> 16) ^ x) * 0x119de1f3;
    x = ((x >> 16) ^ x) * 0x119de1f3;
    x = (x >> 16) ^ x;
    return x;
}


uint64_t magicHash64(uint64_t x)
{
    x = (x ^ (x >> 30)) * (uint64_t)0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * (uint64_t)0x94d049bb133111eb;
    x = x ^ (x >> 31);
    return x;
}

uint64_t magicHash64Reverse(uint64_t x)
{
    x = (x ^ (x >> 31) ^ (x >> 62)) * (uint64_t)0x319642b2d24d8ec3;
    x = (x ^ (x >> 27) ^ (x >> 54)) * (uint64_t)0x96de1b173f119089;
    x = x ^ (x >> 30) ^ (x >> 60);
    return x;
}


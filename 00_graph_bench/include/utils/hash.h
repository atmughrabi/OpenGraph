#ifndef HASH_H
#define HASH_H


#include <stdint.h>

uint32_t magicHash32(uint32_t x);
uint32_t magicHash32Reverse(uint32_t x);
uint64_t magicHash64(uint64_t x);
uint64_t magicHash64Reverse(uint64_t x);


#endif
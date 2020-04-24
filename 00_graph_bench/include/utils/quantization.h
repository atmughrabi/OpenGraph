// -----------------------------------------------------------------------------
//
//      "00_AccelGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Mohannad Ibrahim/Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : bellmanFord.c
// Create : 2019-09-28 15:19:50
// Revise : 2019-09-28 15:34:29
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#ifndef QUANTIZATION_H
#define QUANTIZATION_H

#include <stdint.h>
#include <limits.h>

//ranges for unisigned 8_bit quantization
#define RANGE_MAX_8 UCHAR_MAX
#define RANGE_MIN_8 (uint8_t)0
//ranges for unisigned 16-bit quantization
#define RANGE_MAX_16 USHRT_MAX
#define RANGE_MIN_16 (uint16_t)0
//ranges for unisigned 32-bit quantization
#define RANGE_MAX_32 UINT_MAX
#define RANGE_MIN_32 (uint32_t)0

//ranges for unisigned 32-bit quantization
#define RANGE_MAX UINT_MAX
#define RANGE_MIN (uint32_t)0

struct quant_params
{
    double scale;
    uint32_t zero;  //zero point or  zero-offset
    float min, max; //range
};

struct quant_params_8
{
    float scale;
    uint8_t zero;  //zero point or  zero-offset
    float min, max; //range
};

struct quant_params_16
{
    float scale;
    uint16_t zero;  //zero point or  zero-offset
    float min, max; //range
};

struct quant_params_32
{
    double scale;
    uint32_t zero;  //zero point or  zero-offset
    float min, max, non_zero_min; //range
};

/* function to find min and max values simultanuously amongst the ranks (array)
 it has an O(N) complexity*/
void getMinMax(struct quant_params *q_params, float *ranks, uint32_t size);
/* function to find min and max values simultanuously amongst the ranks (array)
    it has an O(N) complexity*/
void getMinMax_8(struct quant_params_8 *q_params, float *ranks, uint32_t size);
/* function to find min and max values simultanuously amongst the ranks (array)
 it has an O(N) complexity*/
void getMinMax_16(struct quant_params_16 *q_params, float *ranks, uint32_t size);
/* function to find min and max values simultanuously amongst the ranks (array)
    it has an O(N) complexity*/
void getMinMax_32(struct quant_params_32 *q_params, float *ranks, uint32_t size);



#define ABS(num) (double)((num<0)?(-num):(num))
#define ROUND(num) (uint32_t)((num)>=0?((num)+0.5):((num)-0.5))

//to keep the number in the range 0 - 255
#define CLAMP(num,min,max) (num < min ? min : (num > max ? max : num))

//quantization parameters
#define GetScale(min,max) (double)(min == max ? 1.0f : ABS(((double)max - (double)min))/(double)RANGE_MAX)
#define GetZeroPoint(max,scale) (uint32_t)CLAMP((RANGE_MAX - ROUND((double)max/(double)scale)), RANGE_MIN, RANGE_MAX)

//quantize
#define quantize(num,scale,zero) (uint32_t)(CLAMP(ROUND((double)num/(double)scale) + zero, RANGE_MIN, RANGE_MAX))
#define dequantize_f(num,scale,zero) (float)(scale*(CLAMP(num, RANGE_MIN, RANGE_MAX)))
#define dequantize_d(num,scale,zero) (double)(scale*(CLAMP(num, RANGE_MIN, RANGE_MAX)))



#define ABS_8(num) (double)((num<0)?(-num):(num))
#define ROUND_8(num) (uint32_t)((num)>=0?((num)+0.5):((num)-0.5))

//to keep the number in the range 0 - 255
#define CLAMP_8(num,min,max) (num < min ? min : (num > max ? max : num))

//quantization parameters
#define GetScale_8(min,max) (double)(min == max ? 1.0f : ABS_8(((double)max - (double)min))/(double)RANGE_MAX_8)
#define GetZeroPoint_8(max,scale) (uint8_t)CLAMP_8((RANGE_MAX_8 - ROUND_8((double)max/(double)scale)), RANGE_MIN_8, RANGE_MAX_8)

//quantize
#define quantize_8(num,scale,zero) (uint8_t)(CLAMP_8(ROUND_8((double)num/(double)scale) + zero, RANGE_MIN_8, RANGE_MAX_8))
#define dequantize_8_f(num,scale,zero) (float)(scale*(CLAMP_8(num, RANGE_MIN_8, RANGE_MAX_8)))
#define dequantize_8_d(num,scale,zero) (double)(scale*(CLAMP_8(num, RANGE_MIN_8, RANGE_MAX_8)))




#define ABS_16(num) ((double)(num<0)?(-num):(num))
#define ROUND_16(num) (uint32_t)((num)>=0?((num)+0.5):((num)-0.5))

//to keep the number in the range 0 - 255
#define CLAMP_16(num,min,max) (num < min ? min : (num > max ? max : num))

//quantization parameters
#define GetScale_16(min,max) (double)(min == max ? 1.0f : ABS_16(((double)max - (double)min))/(double)RANGE_MAX_16)
#define GetZeroPoint_16(max,scale) (uint16_t)CLAMP_16(RANGE_MAX_16 - ROUND_16(max/scale), RANGE_MIN_16, RANGE_MAX_16)

//quantize_16
#define quantize_16(num,scale,zero) (uint16_t)(CLAMP_16(ROUND_16((double)num/(double)scale) + zero, RANGE_MIN_16, RANGE_MAX_16))
#define dequantize_16_f(num,scale,zero) (float)(scale*(CLAMP_16(num, RANGE_MIN_16, RANGE_MAX_16)))
#define dequantize_16_d(num,scale,zero) (double)(scale*(CLAMP_16(num, RANGE_MIN_16, RANGE_MAX_16)))




#define ABS_32(num) (double)((num<0)?(-num):(num))
#define ROUND_32(num) (uint32_t)((num)>=0?((num)+0.5):((num)-0.5))

//to keep the number in the range 0 - 255
#define CLAMP_32(num,min,max) (num < min ? min : (num > max ? max : num))

//quantization parameters
#define GetScale_32(min,max) (double)(min == max ? 1.0f : ABS_32(((double)max - (double)min))/(double)RANGE_MAX_32)
#define GetZeroPoint_32(max,scale) (uint32_t)CLAMP_32((RANGE_MAX_32 - ROUND_32((double)max/scale)), RANGE_MIN_32, RANGE_MAX_32)

//quantize_32
#define quantize_32(num,scale,zero) (uint32_t)(CLAMP_32(ROUND_32((double)num/(double)scale) + zero, RANGE_MIN_32, RANGE_MAX_32))
#define dequantize_32_f(num,scale,zero) (float)(scale*(CLAMP_32(num, RANGE_MIN_32, RANGE_MAX_32)))
#define dequantize_32_d(num,scale,zero) (double)(scale*(CLAMP_32(num, RANGE_MIN_32, RANGE_MAX_32)))

#endif /* QUANTIZATION_H */
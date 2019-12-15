// -----------------------------------------------------------------------------
//
//		"OpenGraph"
//
// -----------------------------------------------------------------------------
// Copyright (c) 2014-2019 All rights reserved
// -----------------------------------------------------------------------------
// Author : Abdullah Mughrabi
// Email  : atmughra@ncsu.edu||atmughrabi@gmail.com
// File   : timer.c
// Create : 2019-06-21 17:15:17
// Revise : 2019-09-28 15:37:12
// Editor : Abdullah Mughrabi
// -----------------------------------------------------------------------------

#include "timer.h"

void Start(struct Timer *timer)
{
    gettimeofday(&(timer->start_time), NULL);
}

void Stop(struct Timer *timer)
{
    gettimeofday(&(timer->elapsed_time), NULL);
    timer->elapsed_time.tv_sec  -= timer->start_time.tv_sec;
    timer->elapsed_time.tv_usec -= timer->start_time.tv_usec;
}

double Seconds(struct Timer *timer)
{
    return timer->elapsed_time.tv_sec + timer->elapsed_time.tv_usec / 1e6;
}

double Millisecs(struct Timer *timer)
{
    return 1000 * timer->elapsed_time.tv_sec + timer->elapsed_time.tv_usec / 1000;
}

double Microsecs(struct Timer *timer)
{
    return 1e6 * timer->elapsed_time.tv_sec + timer->elapsed_time.tv_usec;
}

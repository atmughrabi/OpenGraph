#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

struct Timer{
  struct timeval start_time;
  struct timeval elapsed_time;
 };

  void Start(struct Timer* timer);
  void Stop(struct Timer* timer);
  double Seconds(struct Timer* timer);
  double Millisecs(struct Timer* timer);
  double Microsecs(struct Timer* timer);


#endif  // TIMER_H
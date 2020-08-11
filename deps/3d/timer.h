#ifndef timer_h
#define timer_h
#include <map>
#include "print.h"
#include <sys/time.h>

namespace exafmm {
  timeval t;
  std::map<std::string,timeval> timer;

  void start(std::string event) {
    gettimeofday(&t, NULL);
    timer[event] = t;
  }

  real_t stop(std::string event) {
    gettimeofday(&t, NULL);
    real_t eventTime = t.tv_sec - timer[event].tv_sec +
      (t.tv_usec - timer[event].tv_usec) * 1e-6;
    print(event, eventTime);
    return eventTime;
  }
}
#endif

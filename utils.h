#ifndef __UTILS_H__
#define __UTILS_H__
#include <ctime>
#include <iomanip>
#include <iostream>

class CTimer {
  clock_t _start;
 public:
  void start() { _start = std::clock(); }
  void stop(std::string msg = "") {
    std::cerr << msg +  " Elapsed: "
	      << (std::clock() - _start) / (double) CLOCKS_PER_SEC
	      << " secs" << std::endl;	      
  }
};


#endif
















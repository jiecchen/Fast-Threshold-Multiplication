#ifndef __UTILS_H__
#define __UTILS_H__
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stack>

class CTimer {
  std::stack<std::string> _msg;
  std::stack<clock_t> _start;
 public:
  void start(std::string msg = "") {
    this->_msg.push(msg);
    _start.push(std::clock());
    for (unsigned int i = 0; i < _start.size(); ++i)
      std::cerr << ">> ";
    std::cerr << msg << std::endl;


  }
  void stop(std::string msg = "") {
    for (unsigned int i = 0; i < _start.size(); ++i)
      std::cerr << ">> ";
    std::cerr << _msg.top() +  " Elapsed: "
	      << (std::clock() - _start.top()) / (double) CLOCKS_PER_SEC
	      << " secs" << std::endl;	     
    _start.pop();
    _msg.pop();
  }
};


#endif
















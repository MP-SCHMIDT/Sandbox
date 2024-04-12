#pragma once

// Standard lib
#include <chrono>
#include <vector>


static std::vector<std::chrono::high_resolution_clock::time_point> TimerStack;

namespace Timer {
  inline int PushTimer() {
    TimerStack.push_back(std::chrono::high_resolution_clock::now());
    return (int)TimerStack.size();
  }

  inline double CheckTimer() {
    if (TimerStack.empty()) return -1.0;
    std::chrono::high_resolution_clock::time_point timeNew= std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point timeOld= TimerStack[TimerStack.size() - 1];
    return std::chrono::duration<double>(timeNew - timeOld).count();
  }

  inline double PopTimer() {
    if (TimerStack.empty()) return -1.0;
    std::chrono::high_resolution_clock::time_point timeNew= std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point timeOld= TimerStack[TimerStack.size() - 1];
    double elapsed= std::chrono::duration<double>(timeNew - timeOld).count();
    TimerStack.pop_back();
    return elapsed;
  }
}  // namespace Timer

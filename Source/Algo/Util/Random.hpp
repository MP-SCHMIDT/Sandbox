#pragma once

// Standard lib
#include <cstdlib>

namespace Random {
  inline int Val(int iMin, int iMax) {
    if (iMax <= iMin) return iMin;
    return iMin + rand() % (iMax - iMin + 1);
  }
  inline float Val(float iMin, float iMax) {
    if (iMax <= iMin) return iMin;
    return iMin + (iMax - iMin) * ((float)rand() / (float)RAND_MAX);
  }
  inline double Val(double iMin, double iMax) {
    if (iMax <= iMin) return iMin;
    return iMin + (iMax - iMin) * ((double)rand() / (double)RAND_MAX);
  }
};

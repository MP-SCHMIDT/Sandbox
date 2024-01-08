#include "Random.hpp"

// Standard lib
#include <cstdlib>


int Random::Val(int const iMin, int const iMax) {
  if (iMax <= iMin) return iMin;
  return iMin + rand() % (iMax - iMin + 1);
}

float Random::Val(float const iMin, float const iMax) {
  if (iMax <= iMin) return iMin;
  return iMin + (iMax - iMin) * ((float)rand() / (float)RAND_MAX);
}

double Random::Val(double const iMin, double const iMax) {
  if (iMax <= iMin) return iMin;
  return iMin + (iMax - iMin) * ((double)rand() / (double)RAND_MAX);
}

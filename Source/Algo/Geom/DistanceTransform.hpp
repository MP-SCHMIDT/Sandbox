#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"


class DistanceTransform
{
  public:
  static void ApplyDistanceTransform(
      const int iMaxIterConvergence,
      const float iSizeX,
      const float iSizeY,
      const float iSizeZ,
      const Field::Field3<char> &iVoxType,
      Field::Field3<float> &ioDist,
      const int iVerbose);
};

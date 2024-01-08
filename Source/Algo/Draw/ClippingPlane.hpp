#pragma once

// Standard lib
#include <array>


class ClippingPlane
{
  public:
  static void DrawClipping(const bool iUse, const int iDim, const double iPos, const bool iSide,
                           const std::array<double, 3> iBoxMin, const std::array<double, 3> iBoxMax);
};

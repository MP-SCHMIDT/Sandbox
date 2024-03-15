#pragma once

// Standard lib
#include <array>
#include <vector>


class MergeVertices
{
  public:
  static void QuadraticMerge(
      double const iTolerance,
      std::vector<std::array<double, 3>>& ioVertices,
      std::vector<std::array<int, 3>>& ioTriangles);
};

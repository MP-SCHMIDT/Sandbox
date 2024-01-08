#pragma once

// Standard lib
#include <array>
#include <vector>


class Bresenham
{
  public:
  static std::vector<std::array<int, 2>> Line2D(int x0, int y0, int x1, int y1);
  static std::vector<std::array<int, 3>> Line3D(int x0, int y0, int z0, int x1, int y1, int z1);
};

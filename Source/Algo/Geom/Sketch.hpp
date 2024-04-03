#pragma once

// Standard lib
#include <array>
#include <vector>


class Sketch
{
  public:
  static void PolylineSubdivideAndSmooth(
      bool const iIsOpenPolyline,
      int const iNbNewNodesPerSegment,
      int const iNbSmoothingSteps,
      std::vector<std::array<double, 3>>& ioPolyline);
};

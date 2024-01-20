#pragma once

// Standard lib
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


class Sketch
{
  public:
  static void PolylineSubdivideAndSmooth(
      bool const iIsOpenPolyline,
      int const iNbNewNodesPerSegment,
      int const iNbSmoothingSteps,
      std::vector<Vec::Vec3<double>>& ioPolyline);
};

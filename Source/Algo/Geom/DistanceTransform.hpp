#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Type/Field.hpp"


class DistanceTransform
{
  public:
  // Distance transform algorithm to spread a narrow band distance field to the full grid
  // Can iterate multiple times to handle inactive voxels acting as obstacle to the distance spread
  // iVoxelType defines if the voxel has a fixed distance (>0) or is inactive (<0)
  static void ApplyDistanceTransform(
      const int iMaxIterConvergence,
      const float iElemSizeX,
      const float iElemSizeY,
      const float iElemSizeZ,
      const Field::Field3<char> &iVoxType,
      Field::Field3<float> &ioDist,
      const int iVerbose);
};

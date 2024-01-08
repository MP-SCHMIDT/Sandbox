#pragma once

// Standard lib
#include <array>
#include <vector>


class MarchingCubes
{
  private:
  static void Interpolate(
      double const iIsoval,
      double const iEpsilon,
      double const iVal1,
      double const iVal2,
      std::array<double, 3> const& iPos1,
      std::array<double, 3> const& iPos2,
      std::array<double, 3>& oPos);

  public:
  static void ComputeMarchingCubes(
      double const iIsoval,
      std::array<double, 3> const& iBBoxMin,
      std::array<double, 3> const& iBBoxMax,
      std::vector<std::vector<std::vector<double>>> const& iField,
      std::vector<std::array<double, 3>>& oVertices,
      std::vector<std::array<int, 3>>& oTriangles);
};

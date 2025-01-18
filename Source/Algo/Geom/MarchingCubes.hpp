#pragma once

// Standard lib
#include <array>
#include <vector>


class MarchingCubes
{
  public:
  static void BuildMesh(const int nX,
                        const int nY,
                        const int nZ,
                        const double iIsoval,
                        const bool isCentered,
                        const bool isPositiveInside,
                        const std::array<double, 3>& iBBoxMin,
                        const std::array<double, 3>& iBBoxMax,
                        const std::vector<double>& iField,
                        std::vector<std::array<double, 3>>& oVertices,
                        std::vector<std::array<int, 3>>& oTriangles);

  private:
  static void Interpolate(const double iIsoval,
                          const double iEpsilon,
                          const double iVal1,
                          const double iVal2,
                          const std::array<double, 3>& iPos1,
                          const std::array<double, 3>& iPos2,
                          std::array<double, 3>& oPos);
};

#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include <Math/FieldArray.hpp>

class DrawField
{
  public:
  static void DrawColored3DField(const std::vector<std::vector<std::vector<bool>>>& iShow,
                                 const std::vector<std::vector<std::vector<std::array<float, 4>>>>& iColor,
                                 const std::array<double, 3>& iBoxMin,
                                 const std::array<double, 3>& iBoxMax,
                                 const std::array<double, 3>& iCamDir,
                                 const bool iWireframe,
                                 const bool iShading,
                                 const bool iAlphaPanels);

  static void DrawColored3DField(const Field3D<char>& iShow,
                                 const Field3D<std::array<float, 4>>& iColor,
                                 const std::array<double, 3>& iBoxMin,
                                 const std::array<double, 3>& iBoxMax,
                                 const std::array<double, 3>& iCamDir,
                                 const bool iWireframe,
                                 const bool iShading,
                                 const bool iAlphaPanels);
};

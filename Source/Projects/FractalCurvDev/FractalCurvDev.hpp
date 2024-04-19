#pragma once

// Standard lib
#include <array>
#include <vector>

// Algo headers
#include "Math/Vec.hpp"


// Generate 3D surfaces from fractal curves
// - Extrude fractals in the out-of-plane rirection to create surfaces
// - Progressively increase fractal depth along the extrusion direction
// - Written to generate 3D printable math shapes
//
// Example for reference
// https://naml.us/paper/fractal/
// https://naml.us/paper/irving_segerman2012_fractal.pdf
class FractalCurvDev
{
  private:
  // List of UI parameters for this project
  enum ParamType
  {
    FractalMode_____,
    MaxDepth________,
    StepZVal________,
    StepZExpo_______,
    SpreadCoeff_____,
    VerboseLevel____,
  };

  std::vector<std::vector<Vec::Vec3<float>>> Nodes;
  std::vector<std::array<Vec::Vec3<float>, 3>> Faces;

  public:
  bool isActivProj;
  bool isAllocated;
  bool isRefreshed;

  FractalCurvDev();

  void SetActiveProject();
  bool CheckAlloc();
  bool CheckRefresh();
  void Allocate();
  void Refresh();
  void KeyPress(const unsigned char key);
  void MousePress(const unsigned char mouse);
  void Animate();
  void Draw();
};
